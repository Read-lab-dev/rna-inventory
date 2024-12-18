##########SVM#######
#multicox
library("survival")
library("survminer")
library(e1071)
library(dplyr)
setwd("~/rna/sc/IMP3/")
selected_feature <-c("RIOK2", "EGFR","IGF2BP3","WWTR1")
selected_feature <-c("RIOK2", "EGFR","IGF2BP3","G3BP1","G3BP2","ILF3","ATXN2L","PABPC4")
selected_feature <-c("RIOK2","IGF2BP3")
selected_feature <-c("RIOK2", "EGFR","IGF2BP3")
load("~/rna/sc/IMP3/TCGA.Rdata")
expr <- as.data.frame(scale(t(exprset[selected_feature,])))
expr$time <- pd$OS.time
expr$status <- pd$OS
res.cox1 <- coxph(as.formula(paste("Surv(time, status) ~",
                            paste(selected_feature,collapse = "+"))),
                 data = expr)
summary(res.cox1)

# Plot the baseline survival function
ggsurvplot(survfit(res.cox1), color = "#2E9FDF",data = expr,
           ggtheme = theme_minimal())

# coefficient <- res.cox1$coefficients
pd$score <- apply(expr[,1:length(selected_feature)], 1, function(x)sum(x*coefficient[1:length(selected_feature)]))


res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "score",minprop = 0.01)
plot(res.cut, "score", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
summary(coxph(Surv(OS.time, OS) ~score, data = res.cat))
p.tcga <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,pval.method	=T,pval.method.coord=c(max(res.cat$OS.time)*0.65-700,1),
                     ggtheme = theme_bw(),
                     risk.table.pos="in",surv.median.line="hv",
                     legend.title="TCGA-IDHwt-GBM",
                     pval=T,pval.coord=c(max(res.cat$OS.time)*0.65,1),palette = c("#ff4d40","#1263ff"))
p.tcga

pd$group <- ifelse(expr$RIOK2>=quantile(expr$RIOK2)[4]&expr$IGF2BP3>=quantile(expr$IGF2BP3)[4],"Coexpression","Non-coexpression")
pd$group <- ifelse(expr$RIOK2>=quantile(expr$RIOK2)[3]&expr$IGF2BP3>=quantile(expr$IGF2BP3)[3],"Coexpression","Non-coexpression")
fit <- survfit(Surv(OS.time, OS) ~group, data = pd)
p.tcga <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,pval.method	=T,pval.method.coord=c(max(pd$OS.time)*0.65-700,1),
                     ggtheme = theme_bw(),
                     risk.table.pos="in",surv.median.line="hv",
                     legend.title="TCGA-IDHwt-GBM",
                     pval=T,pval.coord=c(max(pd$OS.time)*0.65,1),palette = c("#ff4d40","#1263ff"))
# exp$gene2score <- pd$score
# exp$gene2group <- res.cat$score
# exp$gene3score <- pd$score
# exp$gene3group <- res.cat$score
# exp$gene8score <- pd$score
# exp$gene8group <- res.cat$score
# write.table(exp,file = "Remb_rawdata.txt")
##############
load("cgga.Rdata")
expr <- cgga
pd <- pd[!is.na(pd$OS.time),]
pd <- pd %>% filter(IDH_status=="Wildtype")
expr <- as.data.frame(t(cgga[selected_feature,pd$Sample]))
expr$time <- pd$OS.time
expr$status <- pd$OS
res.cox2<- coxph(as.formula(paste("Surv(time, status) ~",paste(selected_feature,collapse = "+"))), data = expr)

# coefficient <- res.cox2$coefficients

pd$score <- apply(expr[,1:length(selected_feature)], 1, function(x)sum(x*coefficient[1:length(selected_feature)]))

res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "score")
plot(res.cut, "score", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
p.cgga <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,pval.method	=T,pval.method.coord=c(max(res.cat$OS.time)*0.75-1000,1),
                     ggtheme = theme_bw(),
                     risk.table.pos="in",surv.median.line="hv",
                     legend.title="CGGA-IDHwt-GBM",
                     pval=T,pval.coord=c(max(res.cat$OS.time)*0.75,1),palette = c("#ff4d40","#1263ff"))
p.cgga
summary(coxph(Surv(OS.time, OS) ~score, data = res.cat))
summary(coxph(Surv(OS.time, OS) ~score, data = pd))
##############
load("Gravendeel.Rdata")
pd <- pd[!is.na(pd$OS.time),]
# pd[is.na(pd$IDH1_status),11] <- "Wild_type"
pd <- pd %>% filter(IDH1_status!="Mut")
expr <-as.data.frame(t(expr[selected_feature,pd$Sample]))
expr$time <- pd$OS.time
expr$status <- pd$OS
res.cox3 <- coxph(as.formula(paste("Surv(time, status) ~",paste(selected_feature,collapse = "+"))), data = expr)


# coefficient <- res.cox3$coefficients
pd$score <- apply(expr[,1:length(selected_feature)], 1, function(x)sum(x*coefficient[1:length(selected_feature)]))

res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "score")
plot(res.cut, "score", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
p.gra <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,pval.method	=T,pval.method.coord=c(max(res.cat$OS.time)*0.75-1000,1),
                     ggtheme = theme_bw(),
                     risk.table.pos="in",surv.median.line="hv",
                     legend.title="Gravendeel-IDHwt-GBM",
                     pval=T,pval.coord=c(max(res.cat$OS.time)*0.75,1),palette = c("#ff4d40","#1263ff"))
p.gra
summary(coxph(Surv(OS.time, OS) ~score, data = res.cat))
summary(coxph(Surv(OS.time, OS) ~score, data = pd))
##############
load("Rembrandt.Rdata")
pd <- pd[!is.na(pd$OS.time),]
pd <- pd %>% filter(Histology=="GBM")
expr <- as.data.frame(t(expr[selected_feature,rownames(pd)]))
expr$time <- pd$OS.time
expr$status <- pd$OS
res.cox4 <- coxph(as.formula(paste("Surv(time, status) ~",paste(selected_feature,collapse = "+"))), data = expr)

coefficient <- res.cox4$coefficients
pd$score <- apply(expr[,1:length(selected_feature)], 1, function(x)sum(x*coefficient[1:length(selected_feature)]))

res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "score")
plot(res.cut, "score", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
p.remb <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,pval.method	=T,pval.method.coord=c(max(res.cat$OS.time)*0.75-1000,1),
                    ggtheme = theme_bw(),
                    risk.table.pos="in",surv.median.line="hv",
                    legend.title="Rembrandt-IDHwt-GBM",
                    pval=T,pval.coord=c(max(res.cat$OS.time)*0.75,1),palette = c("#ff4d40","#1263ff"))
p.remb
summary(coxph(Surv(OS.time, OS) ~score, data = res.cat))
summary(coxph(Surv(OS.time, OS) ~score, data = pd))
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

ggsave("8gene-TCGA-survival.pdf", plot = p.tcga,width = 4.4, height = 3.9)
ggsave("8gene-CGGA-survival.pdf", plot = p.cgga,width = 4.4, height = 3.9)
ggsave("8gene-Gravendeel-survival.pdf", plot = p.gra,width = 4.4, height = 3.9)
ggsave("8gene-Rembrandt-survival.pdf", plot = p.remb,width = 4.4, height = 3.9)


library(timeROC)
with(pd,
     ROC_riskscore <<- timeROC(T = OS.time,
                               delta = OS,
                               marker = score,
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

sink(file = "3gene_output.txt")
print("TCGA")
summary(res.cox1)
print("CGGA")
summary(res.cox2)
print("Gravendeel")
summary(res.cox3)
print("Rembrandt")
summary(res.cox4)
sink(file = NULL)

