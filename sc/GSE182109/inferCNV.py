
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import os as os
sc.settings.set_figure_params(dpi=80, figsize=(5, 5))
os.chdir("/home/hzg/rna/sc/GSE182109/")
adata = sc.read_h5ad("int.seu.h5ad")
locfiles = pd.read_csv("gene.csv")
adata
adata.var['ensg'] = locfiles['SYMBOL'].to_list()
adata.var['chromosome'] = locfiles['chr'].to_list()
adata.var['start'] = locfiles['start'].to_list()
adata.var['end'] = locfiles['end'].to_list()

metadata = pd.read_csv("metadata.csv")
adata.obs["celltype"] = metadata["celltype"].values
adata.obs["celltype"] = adata.obs["celltype"].astype('category')
adata.obs["celltype"]
adata.obs['celltype'] = adata.obs['celltype'].astype('string')
adata.var.loc[:, ["ensg", "chromosome", "start", "end"]].head()

cnv.tl.infercnv(
  adata,
  window_size=350,
  reference_key="celltype",
  reference_cat=["T cell","Oligodendrocyte","Endothelial"],
  )

cnv.pl.chromosome_heatmap(adata, groupby="celltype",log=True,dendrogram=True,save="fig2.png")

cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)

cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True)

cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
cnv.pl.umap(adata,color="cnv_leiden",legend_loc="on data",legend_fontoutline=2,ax=ax1,show=False)
cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata, color="celltype", ax=ax3,show=False,legend_loc='on data',legend_fontsize=10)
cnv.pl.umap(adata, color="cnv_status", ax=ax4,legend_loc='on data')

cnv.pl.umap(adata,color="cnv_leiden",legend_loc="on data",legend_fontoutline=2)
cnv.pl.umap(adata, color="cnv_score")

adata.obs["cnv_status"] = "normal"
adata.obs.loc[
    adata.obs["cnv_leiden"].isin(["12","8","3","11"]), "cnv_status"
] = "tumor"

cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "tumor", :])

cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "normal", :])

cnv.pl.chromosome_heatmap(adata, groupby="cnv_status")

cord = pd.DataFrame(data=adata.obsm['X_cnv_umap'], columns=['x', 'y'])
cord.to_csv('cnv_umap.csv', sep=",")

cord = pd.DataFrame(data=adata.obs['cnv_score'])
cord.to_csv('cnv_score.csv', sep=",")

cord = pd.DataFrame(data=adata.obs['cnv_leiden'])
cord.to_csv('cnv_leiden.csv', sep=",")
