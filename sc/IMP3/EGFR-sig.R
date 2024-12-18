##########EGFR#######
#EGFR
library("survival")
library("survminer")
library(e1071)
library(dplyr)
setwd("~/rna/sc/IMP3/")
load("~/rna/sc/IMP3/TCGA.Rdata")
expr <- as.data.frame(scale(t(exprset["EGFR",])))
expr$time <- pd$OS.time
expr$status <- pd$OS
res.cox <- coxph(Surv(time, status) ~ EGFR, data = expr)
res.cox
coefficient <- res.cox$coefficients
pd$score <- expr$EGFR*coefficient
res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "score")
plot(res.cut, "score", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
p.tcga <- ggsurvplot(fit,risk.table = TRUE,conf.int = TRUE,
                     pval=T,pval.coord=c(1000,1),palette = c("#F76D5E","#26C3FF"))

##############
rm(list = ls())
load("cgga.Rdata")
load("model2.rdata")
expr <- cgga
pd <- pd %>% filter(IDH_status=="Wildtype")
expr <- as.data.frame(scale(t(expr["EGFR",pd$Sample])))
expr$time <- pd$OS.time
expr$status <- pd$OS
res.cox <- coxph(Surv(time, status) ~ EGFR, data = expr)
res.cox
coefficient <- res.cox$coefficients
pd$score <- expr$EGFR*coefficient
res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "score")
plot(res.cut, "score", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
p.tcga <- ggsurvplot(fit,risk.table = TRUE,conf.int = TRUE,
                     pval=T,pval.coord=c(1000,1),palette = c("#F76D5E","#26C3FF"))
##############
load("Gravendeel.Rdata")
pd[is.na(pd$IDH1_status),11] <- "Wild_type"
pd <- pd %>% filter(IDH1_status!="Mut")
expr <-as.data.frame(expr["EGFR",pd$Sample])
colnames(expr)[1]<- "EGFR"
expr$time <- pd$OS.time
expr$status <- pd$OS
res.cox <- coxph(Surv(time, status) ~ EGFR, data = expr)
res.cox
coefficient <- res.cox$coefficients
pd$score <- expr$EGFR*coefficient
res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "score")
plot(res.cut, "score", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
p.gra <- ggsurvplot(fit,risk.table = TRUE,conf.int = TRUE,
                    pval=T,pval.coord=c(1000,1),palette = c("#F76D5E","#26C3FF"))
##############
load("Rembrandt.Rdata")
pd <- pd %>% filter(Histology=="GBM")
expr <- as.data.frame(expr["EGFR",rownames(pd)])
expr$time <- pd$OS.time
expr$status <- pd$OS
colnames(expr)[1]<- "EGFR"
res.cox <- coxph(Surv(time, status) ~ EGFR, data = expr)
res.cox
coefficient <- res.cox$coefficients
pd$score <- expr$EGFR*coefficient
res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "score")
plot(res.cut, "score", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~score, data = res.cat)
p.remb <- ggsurvplot(fit,risk.table = TRUE,conf.int = TRUE,
                     pval=T,pval.coord=c(1000,1),palette = c("#F76D5E","#26C3FF"))
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

ggsave("5Gene-TCGA-survival.pdf", plot = p.tcga,width = 6, height = 6)
ggsave("5Gene-CGGA-survival.pdf", plot = p.cgga,width = 6, height = 6)
ggsave("5Gene-Gravendeel-survival.pdf", plot = p.gra,width = 6, height = 6)
ggsave("5Gene-Rembrandt-survival.pdf", plot = p.remb,width = 6, height = 6)

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

4.5*4.5