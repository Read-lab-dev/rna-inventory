library(dplyr)
library(glmnet)
library(survival)
library(survminer)
rm(list = ls())
load("~/rna/sc/IMP3/TCGA.Rdata")
surv <- pd[,c(2,4)]
cox.pcutoff <- 0.05 # cox的p阈值
Coxoutput.OS <- NULL

DEG_RIOK2 <- read.csv("./bulk/DEseq2_GBM301_shRIOK2.csv")
DEG_IMP3 <- read.csv("./bulk/DEseq2_GBM301_sh77.csv")

select_RIOK2<- DEG_RIOK2 %>% filter(change!="not")
select_IMP3<- DEG_IMP3 %>% filter(change!="not")
select <- inner_join(select_IMP3,select_RIOK2,by="X")
select <- select %>% filter(log2FoldChange.x*log2FoldChange.y>0)
expr<-exprset[rownames(exprset)%in%select$X,]
display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  
for (i in 1:nrow(expr)) {
  display.progress(index = i,totalN = nrow(expr)) # 显示进度
  # 产生临时变量存储生存以及变量表达值
  tmp <- data.frame(gene = as.numeric(expr[i,]),
                    OS.time = surv[,"OS.time"],
                    OS = surv[,"OS"],
                    stringsAsFactors = F)
  
  # 单变量cox比例风险模型
  cox <- coxph(Surv(OS.time, OS) ~ gene, data = tmp)
  coxSummary = summary(cox)
  
  # 生成cox结果数据框，包括基因名，风险比，z值，waldtest p值，以及HR置信区间
  Coxoutput.OS=rbind.data.frame(Coxoutput.OS,data.frame(gene=rownames(expr)[i],
                                                        HR=as.numeric(coxSummary$coefficients[,"exp(coef)"]),
                                                        z=as.numeric(coxSummary$coefficients[,"z"]),
                                                        pvalue=as.numeric(coxSummary$coefficients[,"Pr(>|z|)"]),
                                                        lower=as.numeric(coxSummary$conf.int[,3]),
                                                        upper=as.numeric(coxSummary$conf.int[,4]),
                                                        stringsAsFactors = F),
                                stringsAsFactors = F)
}
# gene.sel <- c("RIOK2", "EGFR","IGF2BP3","CRTC2")
gene.sel <- Coxoutput.OS[which(Coxoutput.OS$pvalue < cox.pcutoff),"gene"]
expr <- t(exprset[gene.sel,])
colnames(surv) <- c("status","time")
surv <- Surv(surv$time,surv$status)
fit = glmnet(expr,surv,alpha = 1,family = "cox")
plot(fit,xvar = "lambda",label = T)
cvfit <- cv.glmnet(expr,surv,alpha = 1,family = "cox")
plot(cvfit)
coefficient <- coef(cvfit,s=cvfit$lambda.min)
selected_index <- which(as.numeric(coefficient)!=0)
selected_feature <- names(coefficient[selected_index,])

risk_score <- apply(expr[,selected_feature], 1, function(x)sum(x*coefficient[selected_index,]))
pd$riskscore <- as.numeric(risk_score)
pd$risk <- ifelse(pd$riskscore<=median(pd$riskscore),"low","high")
surv.fit <- survfit(Surv(OS.time,OS)~risk,data = pd)
# res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "IMP3",minprop = 0.4)
# plot(res.cut, "IMP3", palette = "npg")
# res.cat <- surv_categorize(res.cut)
# fit <- survfit(Surv(OS.time, OS) ~IMP3, data = res.cat)
p.tcga <- ggsurvplot(surv.fit,risk.table = TRUE,conf.int = TRUE,
           pval=T,pval.coord=c(1000,1),palette = c("#F76D5E","#26C3FF"))+
  ggtitle("Training Set:TCGA GBM")
save(coefficient,selected_index,selected_feature,file = "model.rdata")

cox <- coxph(Surv(OS.time, OS) ~ riskscore, data = pd)


p.tcga
######CGGA#########
rm(list = ls())
load("cgga.Rdata")
load("model2.rdata")
pd <- pd %>% filter(IDH_status=="Wildtype")
# selected_feature[28] <- "FAM35A"
cgga <- scale(t(cgga[selected_feature,pd$Sample]))
risk_score <- apply(cgga, 1, function(x)sum(x*coefficient[selected_index,]))
pd$riskscore <- as.numeric(risk_score)
pd <- pd[!is.na(pd$OS.time),]
pd$risk <- ifelse(pd$riskscore<=median(pd$riskscore),"low","high")
surv.fit <- survfit(Surv(OS.time,OS)~risk,data = pd)
p.cgga <-ggsurvplot(surv.fit,risk.table = TRUE,conf.int = TRUE,
                    pval=T,palette= c("#F76D5E","#26C3FF"),
                    pval.coord=c(max(pd$OS.time)/2,1))+
  ggtitle("CGGA GBM")
p.cgga
######GRA#########
rm(list = ls())
load("Gravendeel.Rdata")
load("model2.rdata")
# selected_feature[28] <- "FAM35A"
pd[is.na(pd$IDH1_status),11] <- "Wild_type"
pd <- pd %>% filter(IDH1_status!="Mut")
idx <- selected_feature%in%rownames(expr)
exprdat <- scale(t(expr[selected_feature[idx],rownames(pd)]))
risk_score <- apply(exprdat, 1, function(x)sum(x*coefficient[selected_index[idx],]))
pd$riskscore <- as.numeric(risk_score)
pd <- pd[!is.na(pd$OS.time),]
pd$risk <- ifelse(pd$riskscore<=median(pd$riskscore),"low","high")
surv.fit <- survfit(Surv(OS.time,OS)~risk,data = pd)
p.gra <-ggsurvplot(surv.fit,risk.table = TRUE,conf.int = TRUE,xlim = c(0,2200),
                   pval=T,palette= c("#F76D5E","#26C3FF"),
                   pval.coord=c(max(pd$OS.time)/5,1))+
                   ggtitle("Gravendeel et.al 2009 GBM")
p.gra

#####REMBRANDT###
# expr <- fread("2024-06-03_Rembrandt_exp.txt")
# expr <- as.data.frame(expr)
# expr <- tibble::column_to_rownames(expr,var = "Sample")
# expr <- t(expr)
# pd <- read.table("2024-06-03_Rembrandt_pheno.txt",header = T)
# pd <- tibble::column_to_rownames(pd,var = "Sample")
# save(pd,expr,file="Rembrandt.Rdata")

load("Rembrandt.Rdata")
load("model2.rdata")
pd <- pd %>% filter(Histology=="GBM")
idx <- selected_feature%in%rownames(expr)
exprdat <- scale(t(expr[selected_feature[idx],rownames(pd)]))
risk_score <- apply(exprdat, 1, function(x)sum(x*coefficient[selected_index[idx],]))
pd$riskscore <- as.numeric(risk_score)
pd$risk <- ifelse(pd$riskscore<=median(pd$riskscore),"low","high")
surv.fit <- survfit(Surv(OS.time,OS)~risk,data = pd)
p.remb <-ggsurvplot(surv.fit,risk.table = TRUE,conf.int = TRUE,xlim = c(0,3000),
                    pval=T,palette = c("#F76D5E","#26C3FF"),
                    pval.coord=c(1000,1))+
           ggtitle("Rembrandt et.al 2009 GBM")
p.remb

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

ggsave("M3-TCGA-survival.pdf", plot = p.tcga,width = 6, height = 6)
ggsave("M3-CGGA-survival.pdf", plot = p.cgga,width = 6, height = 6)
ggsave("M3-Gravendeel-survival.pdf", plot = p.gra,width = 6, height = 6)
ggsave("M3-Rembrandt-survival.pdf", plot = p.remb,width = 6, height = 6)

sig <- data.frame(row.names = NULL,feature=selected_feature,
                  coef=coefficient[selected_index,])
DT::datatable(sig)
write.csv(sig,file = "sig-M3.csv")

#############ROC############
library(timeROC)
with(pd,
     ROC_riskscore <<- timeROC(T = OS.time,
                               delta = OS,
                               marker = riskscore,
                               cause = 1,
                               weighting = "marginal",
                               times = c(365,730,1095),
                               ROC = TRUE,
                               iid = TRUE)
)
plot(ROC_riskscore, time = 365,lwd=2, col = "#F76D5E", add = F,title = "")
plot(ROC_riskscore, time = 730,lwd=2,col = "#72D8FF", add = T)
plot(ROC_riskscore, time = 1095,lwd=2 ,col = "#BA55D3", add = T)
legend("bottom",bty='n',bg = NA,
       c(paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3))
                       ,paste("2-Year AUC = ",round(ROC_riskscore$AUC[2],3))
                        ,paste("3-Year AUC = ",round(ROC_riskscore$AUC[3],3))),col=c("#F76D5E","#72D8FF","#BA55D3"),lty=1,lwd=2)

4.5*4.5

#########SUBTYPE#########
my_comparisons <- list(c("Mesenchymal", "Classical"), 
                       c("Classical", "Proneural"),
                       c("Mesenchymal", "Proneural"))
ggplot(pd,aes(Subtype_Verhaak_2010,riskscore,fill=Subtype_Verhaak_2010))+
  stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA,alpha=0.9)+
  geom_jitter(width = 0.2,color="black",alpha=0.5)+
  ylab("Riskscore")+
  scale_fill_manual(values = c("#F76D5E", "#FFFFBF", "#72D8FF","purple"))+
  stat_compare_means(method="wilcox",comparisons=my_comparisons,
                     label="p.signif")+
  theme_few()+
  theme(legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text  = element_text(face = "bold"),
        axis.title.x = element_blank())
ggsave(file="Gravendeel-SUBTYPE-Boxplot.pdf",width = 4,height = 4)
