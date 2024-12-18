import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os as os
import matplotlib; matplotlib.use('TkAgg')
os.chdir("/home/hzg/rna/xenium/liu/")
sample='MUV1'
path='/home/hzg/rna/xenium/liu/reads/DECODED_'+sample+'.csv'


re=pd.read_csv(path)
pciseq=pd.read_csv('reads/'+sample+'_pciseq_results.csv')
re=re.rename(columns={'x':'X','y':'Y'})
re.shape
data=re.merge(pciseq,on=['X','Y'],how='outer')
data.shape

metadata=data.iloc[:, list(range(1, 29)) + list(range(136, 142))]
data.iloc[:,136:143]
expdata=data.iloc[:,range(28,135)]

adata = sc.AnnData(expdata)
adata.obs=metadata

maximum_numbers_of_reads=50
minimum_number_of_reads=5

adata.obs['total_counts']=list(np.sum(expdata.loc[:,:],axis=1))
sns.displot(adata.obs,x='total_counts')
plt.axvline(x=maximum_numbers_of_reads,color='green',linestyle='--')
plt.axvline(x=minimum_number_of_reads,color='red',linestyle='--')
plt.show()
adata=adata[adata.obs['total_counts']>minimum_number_of_reads] #adjust the numbers to keep good cells according to appropiate ones
sc.pp.filter_cells(adata,min_genes=3)

adata.raw=adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

plt.rcParams['figure.facecolor'] = 'white'
sc.pp.highly_variable_genes(adata, min_mean=0.3, max_mean=7, min_disp=-0.5)
sc.pl.highly_variable_genes(adata)
adata=adata[:,adata.var.highly_variable==True]

ames=adata[adata.obs['name']=='MES_like']
aoc=adata[adata.obs['name']=='OC_like']
aopc=adata[adata.obs['name']=='OPC_like']
aac=adata[adata.obs['name']=='AC_like']
acycling=adata[adata.obs['name']=='Cycling']

plt.figure(figsize=(15,10))
sns.kdeplot(ames.obs['AC_like'], bw=0.05,color='green')
sns.kdeplot(aoc.obs['AC_like'], bw=0.05,color='blue')
sns.kdeplot(aopc.obs['AC_like'], bw=0.05,color='red')
sns.kdeplot(aac.obs['AC_like'], bw=0.05,color='orange')
sns.kdeplot(acycling.obs['AC_like'], bw=0.05,color='pink')
plt.show()

adata=sc.read("adata_16samples_DMG.h5ad")

adata=adata[~adata.obs['name'].isin(['nan'])]

adata.obsm["spatial"]=np.array([adata.obs.X*0.325,adata.obs.Y*0.325]).transpose().astype('float64')

for samp in adata.obs['sample'].unique():
    asub=adata[adata.obs['sample']==samp]
    print(samp)
    sc.pl.spatial(
    asub,
    color="name",
    neighbors_key="spatial_neighbors",
    spot_size=40,
    edges=False,
    edges_width=2,
    img_key=None,
    )