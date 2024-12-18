import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
sc.settings.set_figure_params(dpi=80, figsize=(5, 5))
os.chdir('/home/hzg/backup/lab/BT109VP/Rcnv/')

adata = sc.read_h5ad("/home/hzg/backup/lab/BT109VP/Rcnv/int.seu.h5ad")
locfiles = pd.read_csv("/home/hzg/backup/lab/BT109VP/Rcnv/gene.csv")

adata.var['ensg'] =locfiles['ensg'] .to_list()
adata.var['chromosome'] = locfiles['chr'].to_list()
adata.var['start'] = locfiles['start'].to_list()
adata.var['end'] = locfiles['end'].to_list()

adata.var.loc[:, ["ensg", "chromosome", "start", "end"]].head()

adata.obs['seurat_clusters'] = adata.obs['seurat_clusters'].astype('string')
adata.obs['celltype'] = adata.obs['celltype'].astype('string')

sc.pl.umap(adata, color="seurat_clusters")

cnv.tl.infercnv(
  adata,
  window_size=250, 
  reference_key="celltype",
  reference_cat=["8"],
  )

cnv.pl.chromosome_heatmap(adata, groupby="celltype")

cnv.tl.pca(adata)
cnv.pp.neighbors(adata)
cnv.tl.leiden(adata)

cnv.pl.chromosome_heatmap(adata, groupby="cnv_leiden", dendrogram=True)

cnv.tl.umap(adata)
cnv.tl.cnv_score(adata)
cnv.pl.umap(adata,color="seurat_clusters")

metadata = pd.read_csv("metadata.csv")

adata.obs["cell_type"] = metadata["seurat_clusters"].values
adata.obs["cell_type"] = adata.obs["cell_type"].astype('category')

adata.obs["cell_type"]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
cnv.pl.umap(adata,color="cnv_leiden",legend_loc="on data",legend_fontoutline=2,ax=ax1,show=False)
cnv.pl.umap(adata, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata, color="celltype", ax=ax3)
cnv.pl.umap(adata, color="cnv_status", ax=ax4)

adata.obs["cnv_status"] = "normal"
adata.obs.loc[
adata.obs["cnv_leiden"].isin(["1", "7","16","17","2","15","8","0"]), "cnv_status"
] = "tumor"

cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "tumor", :])

cnv.pl.chromosome_heatmap(adata[adata.obs["cnv_status"] == "normal", :])

cnv.pl.chromosome_heatmap(adata)

cord = pd.DataFrame(data=adata.obsm['X_cnv_umap'], columns=['x', 'y'])
cord.to_csv('cnv_umap.csv', sep=",")

cord = pd.DataFrame(data=adata.obs['cnv_score'])
cord.to_csv('cnv_score.csv', sep=",")

cord = pd.DataFrame(data=adata.obs['cnv_leiden'])
cord.to_csv('cnv_leiden.csv', sep=",")
