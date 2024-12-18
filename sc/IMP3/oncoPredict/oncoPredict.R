install.packages("oncoPredict")
rm(list = ls())  ## 魔幻操作，一键清空~
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
setwd('./oncoPredict/')
GDSC2_Expr = readRDS(file="GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res = readRDS(file ="GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 

#######
#CTRP2
GDSC2_Expr = readRDS(file = "CTRP2_Expr (TPM, not log transformed).rds")
GDSC2_Res = readRDS(file = "CTRP2_Res.rds")
#######

gene_length <- read.csv("~/rna/human_genome/gene_length.csv")
length <- clusterProfiler::bitr(gene_length$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
length <- inner_join(gene_length,length)
exprSet$SYMBOL <- rownames(exprSet)
length <- left_join(exprSet,length)[,c(7,9)]
length <- length[!duplicated(length$SYMBOL),]
exprSet <- exprSet[,-7]
kb <- length$Length/1000
rpk <- na.omit(exprSet / kb)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)

testExpr<- log2(tpm+1)
# testExpr<- tpm
testExpr[1:4,1:4]  
dim(testExpr)  

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'standardize',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

resultPtype <- read.csv('./calcPhenotype_Output/DrugPredictions-BT109-GDSC.csv', header = T ,stringsAsFactors = F ,row.names = 1)

km_drug <- cbind(metadata,resultPtype) 
km_drug[1:4,1:4]
# km_drug <- km_drug[-4,]
km_drug_melt <- melt(km_drug)

idx <- as.character(unique(km_drug_melt$variable))
outtab <- NULL

for (i in idx) {
  aa <-km_drug_melt %>% filter(variable==i)
  wilcox <- t.test(value~Treat,data=aa,paired=T)
  wilcox.data <- data.frame(pval=wilcox$p.value,
                            W=wilcox$statistic,
                            drug=i)
  outtab <- rbind(outtab,wilcox.data)
}
outtab <- outtab %>% filter(pval<0.05) %>% arrange(desc(W))

outtab_sen <- outtab %>% filter(W>0&pval<0.05)
outtab_res <- outtab %>% filter(W<0&pval<0.05)


my_comparisons <- list(c("DMSO", "VP"))

ggplot(km_drug_melt %>% filter(variable%in%outtab$drug),aes(Treat,value,fill=Treat))+
  stat_boxplot(geom = "errorbar",width=0.6)+
  geom_boxplot()+
  facet_wrap(vars(variable),scales="free",ncol = 5)+
  stat_compare_means(comparisons = my_comparisons,method = "t.test",label="p.signif",paired = T)+
  theme(strip.background = element_rect(
    color="black", fill="#FC4E07AA", linetype="solid",linewidth = 0.1
  ))+scale_fill_manual(values = c("#CF5A79","#544799"))
ggsave("Predicted_Drug_BT109_VP_CTPC.pdf",height = 10,width = 10)

########

