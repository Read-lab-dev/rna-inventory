# -*- coding: utf-8 -*-
import scanpy as sc
import scvelo as scv
import pandas as pd
import os
import anndata
import numpy as np
import matplotlib.pyplot as plt
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo')  # for beautified visualization
scv.settings.set_figure_params('scvelo')
os.chdir('/home/hzg/rna/BT109VP/ac_analysis/')
adata = sc.read("/home/hzg/rna/BT109VP/velocyto/velo.h5ad")

sample_obs = pd.read_csv("cellID_obs_ac.csv")
umap = pd.read_csv("cell_embeddings_ac.csv")
cell_clusters = pd.read_csv("clusters_ac.csv")

adata = adata[np.isin(adata.obs.index,sample_obs["x"])]
adata_index = pd.DataFrame(adata.obs.index)
umap = umap.rename(columns = {'Unnamed: 0':'CellID'})
adata_index = adata_index.rename(columns = {0:'CellID'})
umap_ordered = adata_index.merge(umap, on = "CellID")
umap_ordered = umap_ordered.iloc[:,1:]
adata.obsm['X_umap'] = umap_ordered.values
adata.var_names_make_unique()
adata.var.index.is_unique

cell_clusters_ordered = adata_index.merge(cell_clusters, on = "CellID")
adata.obs['clusters'] = pd.Categorical(cell_clusters_ordered['x'])

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)
scv.tl.velocity(adata, mode = "stochastic")
scv.tl.velocity_graph(adata,n_jobs=10)
scv.pl.velocity_embedding_stream(adata, basis = 'umap',color="clusters",palette="Set3",save='velo_umap.pdf',density=1.5)
