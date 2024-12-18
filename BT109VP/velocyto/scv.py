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
os.chdir('/home/hzg/rna/BT109VP/velocyto/')
#############################################
adata = anndata.read_loom('21047FL-89-01-02.loom')
sample_obs = pd.read_csv("cellID_obs.csv")
umap = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters.csv")

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
scv.pl.velocity_embedding(adata, basis = 'umap',color="clusters")
scv.pl.velocity_embedding_stream(adata, basis = 'X_pca',color="clusters")

scv.pl.velocity_embedding_stream(adata, basis = 'umap',color="clusters",palette="Set2",save='velo_umap.pdf',density=3)


adata.write("velo.h5ad", compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading

scv.pl.velocity(adata, ['CNN3',  'SLC1A2', 'TET1', 'GLI3'], ncols=2)

scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

sc.pl.umap(adata,color="clusters",palette="Set1",legend_loc="on data")
scv.tl.score_genes_cell_cycle(adata)
scv.pl.scatter(adata, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95])

s_genes, g2m_genes = scv.utils.get_phase_marker_genes(adata)
s_genes = scv.get_df(adata[:, s_genes], 'spearmans_score', sort_values=True).index
g2m_genes = scv.get_df(adata[:, g2m_genes], 'spearmans_score', sort_values=True).index

kwargs = dict(frameon=False, ylabel='cell cycle genes')
scv.pl.scatter(adata, list(s_genes[:2]) + list(g2m_genes[:3]), **kwargs)

scv.pl.velocity(adata,['MET','CENPF'],ncols=2)

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])

import seaborn as sns
pd.set_option('display.latex.repr', True)
df = adata.obs.groupby('clusters')[keys].mean().T

sns.clustermap(df,annot=True,cmap='coolwarm')
sc.pl.umap(adata,color="clusters",palette="Set1",legend_loc="on data")

scv.pl.velocity_graph(adata, threshold=.1,palette="Set2")

x, y = scv.utils.get_cell_transitions(adata, basis='umap', starting_cell='21047-120-03-01-01:TGTTGAGTCTGAGAAAx')
ax = scv.pl.velocity_graph(adata, c='lightgrey', edge_width=.05, show=False)
ax = scv.pl.scatter(adata, x=x, y=y, s=120, c='ascending', cmap='gnuplot', ax=ax)

scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='clusters')
df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
sns.heatmap(df,annot=True,cmap='Blues')

scv.pl.paga(adata, basis='umap', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5)

scv.tl.recover_dynamics(adata,n_jobs=10)
scv.tl.latent_time(adata)
scv.pl.umap(adata, color="latent_time", color_map="gnuplot")

top_genes = adata.var["fit_likelihood"].sort_values(ascending=False).index[:300]
top_genes = ["EGFR","CD44","PROM1","SOX2","NEU1","NES","GFAP","MET","COL1A1","PDGFRA","DCX","MAP2"]
scv.pl.heatmap(adata, var_names=top_genes, sortby="velocity_pseudotime", col_color="clusters", n_convolve=100)
scv.pl.umap(adata, color=['DCX','velocity_pseudotime'], color_map="coolwarm")
scv.pl.scatter(adata,x='DCX',y='velocity_pseudotime')
adata.write("velo.h5ad", compression='gzip')  # `compression='gzip'` saves disk space, but slows down writing and subsequent reading