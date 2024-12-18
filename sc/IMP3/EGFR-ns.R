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
pd$EGFR <- expr$EGFR
res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "EGFR")
plot(res.cut, "EGFR", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~EGFR, data = res.cat)

p.tcga <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,
                     ggtheme = theme_bw(),
                     pval=T,pval.coord=c(max(res.cat$OS.time)*0.75,1),palette = c("#ff4d40","#1263ff"))
p.tcga

##############
load("cgga.Rdata")
expr <- cgga
pd <- pd[!is.na(pd$OS.time),]
pd <- pd %>% filter(IDH_status=="Wildtype")
expr <- as.data.frame(scale(t(expr["EGFR",pd$Sample])))
pd$EGFR <- expr$EGFR
res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "EGFR")
plot(res.cut, "EGFR", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~EGFR, data = res.cat)
p.cgga <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,
                     ggtheme = theme_bw(),
                     pval=T,pval.coord=c(max(res.cat$OS.time)*0.8,1),palette = c("#ff4d40","#1263ff"))
##############
load("Gravendeel.Rdata")
pd <- pd[!is.na(pd$OS.time),]
# pd[is.na(pd$IDH1_status),11] <- "Wild_type"
pd <- pd %>% filter(IDH1_status!="Mut")
expr <-as.data.frame(expr["EGFR",pd$Sample])
colnames(expr)[1]<- "EGFR"
pd$EGFR <- expr$EGFR
res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "EGFR")
plot(res.cut, "EGFR", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~EGFR, data = res.cat)
p.gra <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,
                    ggtheme = theme_bw(),
                    pval=T,pval.coord=c(max(res.cat$OS.time)*0.7,1),palette = c("#ff4d40","#1263ff"))
##############
load("Rembrandt.Rdata")
pd <- pd[!is.na(pd$OS.time),]
pd <- pd %>% filter(Histology=="GBM")
expr <- as.data.frame(expr["EGFR",rownames(pd)])
expr$time <- pd$OS.time
expr$status <- pd$OS
colnames(expr)[1]<- "EGFR"
pd$EGFR <- expr$EGFR
res.cut <-surv_cutpoint(pd,time = "OS.time",event = "OS",variables = "EGFR")
plot(res.cut, "EGFR", palette = "npg")
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(OS.time, OS) ~EGFR, data = res.cat)
p.remb <- ggsurvplot(fit,risk.table = T,conf.int = TRUE,
                     ggtheme = theme_bw(),
                     pval=T,pval.coord=c(max(res.cat$OS.time)*0.75,1),palette = c("#ff4d40","#1263ff"))
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

ggsave("EGFR-TCGA-survival.pdf", plot = p.tcga,width = 6, height = 6)
ggsave("EGFR-CGGA-survival.pdf", plot = p.cgga,width = 6, height = 6)
ggsave("EGFR-Gravendeel-survival.pdf", plot = p.gra,width = 6, height = 6)
ggsave("EGFR-Rembrandt-survival.pdf", plot = p.remb,width = 6, height = 6)

library(timeROC)
with(pd,
     ROC_riskEGFR <<- timeROC(T = OS.time,
                               delta = OS,
                               marker = EGFR,
                               cause = 1,
                               weighting = "marginal",
                               times = c(365,730,1095),
                               ROC = TRUE,
                               iid = TRUE)
)
plot(ROC_riskEGFR, time = 365,lwd=2, col = "#F76D5E", add = F,title = "")
plot(ROC_riskEGFR, time = 730,lwd=2,col = "#72D8FF", add = T)
plot(ROC_riskEGFR, time = 1095,lwd=2 ,col = "#BA55D3", add = T)
legend("bottom",bty='n',bg = NA,
       c(paste("1-Year AUC = ",round(ROC_riskEGFR$AUC[1],3))
         ,paste("2-Year AUC = ",round(ROC_riskEGFR$AUC[2],3))
         ,paste("3-Year AUC = ",round(ROC_riskEGFR$AUC[3],3))),col=c("#F76D5E","#72D8FF","#BA55D3"),lty=1,lwd=2)

4.5*4.5