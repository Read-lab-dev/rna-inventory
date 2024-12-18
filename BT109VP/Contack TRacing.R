##Contact Tracing
##Zhengang Hu 12-07-2023 15.11
library(Seurat)
library(SeuratDisk)
int.seu <- qs::qread("int.qs")
Idents(int.seu) <- int.seu$orig.ident
bb <- subset(int.seu,idents=c("Eng_DMSO","Eng_VP_S","Eng_VP_L"))
bb$condition <- ifelse(bb$orig.ident=="Eng_DMSO","DMSO","VP")
bb$sample <- bb$orig.ident
cc <- GetAssayData(bb,slot = "counts")
output.seu <- CreateSeuratObject(counts = cc,meta.data = bb@meta.data,assay = "RNA")
output.seu$celltype <- as.character(output.seu$celltype)
SaveH5Seurat(output.seu, filename = "int.h5Seurat")
Convert("int.h5Seurat",dest = "h5ad")

output.seu <- standard10X(output.seu)
DimPlot(output.seu)
