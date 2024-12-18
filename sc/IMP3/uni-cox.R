display.progress = function (index, totalN, breakN=20) {
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
}  
Coxoutput.OS <- NULL
load("~/rna/sc/IMP3/TCGA.Rdata")
expr <- as.data.frame(scale(exprset[selected_feature,]))


load("cgga.Rdata")
expr <- cgga
pd <- pd[!is.na(pd$OS.time),]
pd <- pd %>% filter(IDH_status=="Wildtype")
expr <- as.data.frame(cgga[selected_feature,pd$Sample])

load("Gravendeel.Rdata")
pd <- pd[!is.na(pd$OS.time),]
# pd[is.na(pd$IDH1_status),11] <- "Wild_type"
pd <- pd %>% filter(IDH1_status!="Mut")
expr <-as.data.frame(expr[selected_feature,pd$Sample])

load("Rembrandt.Rdata")
pd <- pd[!is.na(pd$OS.time),]
pd <- pd[!is.na(pd$OS),]
pd <- pd %>% filter(Histology=="GBM")
expr <- as.data.frame(expr[selected_feature,rownames(pd)])
for (i in 1:nrow(expr)) {
  display.progress(index = i,totalN = nrow(expr)) # 显示进度
  # 产生临时变量存储生存以及变量表达值
  tmp <- data.frame(gene = as.numeric(expr[i,]),
                    OS.time = pd[,"OS.time"],
                    OS = pd[,"OS"],
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
pd$IMP3 <- as.numeric(t(expr["IGF2BP3",]))
pd$group <- ifelse(pd$IMP3<=median(pd$IMP3),"low-IMP3","high-IMP3")
fit <- survfit(Surv(OS.time, OS) ~group, data = pd)
p.tcga <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,pval.method	=T,pval.method.coord=c(max(pd$OS.time)*0.65-700,1),
                     ggtheme = theme_bw(),
                     risk.table.pos="in",surv.median.line="hv",
                     legend.title="Rembrandt-IDHwt-GBM     HR=1.55",
                     pval=T,pval.coord=c(max(pd$OS.time)*0.65,1),
                     palette = c("#ff4d40","#1263ff"))
p.tcga
