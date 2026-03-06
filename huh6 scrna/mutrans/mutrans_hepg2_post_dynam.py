## MuTrans HepG2 post-treatment: landscapes
## August 2025

## [ Load dependencies ] ----

import sys
import os

print(f"Current directory: {os.getcwd()}")
os.chdir("/data/rds/DMP/DUDMP/PAEDCANC/echen/MuTrans-release-main/Example")
print(f"Changed directory: {os.getcwd()}")

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import seaborn as sns
import hdf5plugin
import gc

import matplotlib.pyplot as plt
import pyMuTrans as pm

import matlab
import matlab.engine

sc.settings.set_figure_params(dpi=100, frameon=False)

datadir = "data_hepg2/" 
plotdir = "../plots/hepg2_post/" 

filename = "post_scanpy_nn.h5ad"
file_path = os.path.join(datadir, filename)
print(file_path)
sample_name = filename.replace("_scanpy_nn.h5ad", "")
adata = sc.read_h5ad(file_path)

## [ Dyanamical analysis ] ----

k = 6
par = {"choice_distance":"cosine", "K_cluster":k, "trials":50, "weight_scale":True, "initial":"pca", "reduce_large_scale":True, "reduce_num_meta_cell":1000.0, "fig_name": f"{plotdir}{sample_name}_dynam_out"} 
adata_mu = pm.dynamical_analysis(adata, par)
out = adata_mu.uns['da_out']
ind = np.argsort(np.asarray(out['perm_class']).ravel().astype(int)-1)
labels_in_meta = np.asarray(out['reduce_class']).ravel().astype(int)-1
attractor_meta = np.asarray(out['class_order']).ravel()[ind]-1
adata_mu.obs['attractor'] = [ int(attractor_meta[labels_in_meta[i]]) for i in range(len(labels_in_meta))]
entropy_meta = np.asarray(out['H']).ravel()[ind]
adata_mu.obs['entropy'] = [entropy_meta[labels_in_meta[i]] for i in range(len(labels_in_meta))]

print(adata_mu)

print(np.asarray(adata_mu.uns['da_out']['mu_hat']))

## [ MPFT and MPPT ] ----

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, method = "MPFT", size_point = 40, alpha_point = 0.5, size_text = 15)
plt.savefig(f"{plotdir}{sample_name}_mpft.pdf")
plt.savefig(f"{plotdir}{sample_name}_mpft.png")

## 0 to 4
si = 0
sf = 4
fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si,sf = sf, flux_fraction = 1, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.png")

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si, sf = sf, flux_fraction = 0.9, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.png")

## 0 to 5
si = 0
sf = 5
fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si,sf = sf, flux_fraction = 1, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.png")

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si, sf = sf, flux_fraction = 0.9, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.png")

## 1 to 4
si = 1
sf = 4
fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si,sf = sf, flux_fraction = 1, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.png")

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si, sf = sf, flux_fraction = 0.9, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.png")

## 1 to 5
si = 1
sf = 5
fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si,sf = sf, flux_fraction = 1, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.png")

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si, sf = sf, flux_fraction = 0.9, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.png")

## 2 to 4
si = 2
sf = 4
fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si,sf = sf, flux_fraction = 1, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.png")

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si, sf = sf, flux_fraction = 0.9, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.png")

## 2 to 5
si = 2
sf = 5
fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si,sf = sf, flux_fraction = 1, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.png")

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si, sf = sf, flux_fraction = 0.9, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.png")

## 3 to 4
si = 3
sf = 4
fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si,sf = sf, flux_fraction = 1, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.png")

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si, sf = sf, flux_fraction = 0.9, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.png")

## 3 to 5
si = 3
sf = 5
fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si,sf = sf, flux_fraction = 1, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_100.png")

fig = plt.figure(figsize = (10, 6))
pm.infer_lineage(adata_mu, si = si, sf = sf, flux_fraction = 0.9, method = "MPPT", size_state = 0.2, size_point = 40, size_text = 20, alpha_point = 0.7)
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.pdf")
plt.savefig(f"{plotdir}{sample_name}_mppt_{si}_to_{sf}_90.png")

## [ Attractor state plots ] ----

adata_mu.obs['attractor'] = adata_mu.obs.attractor.astype('category')

## Entropy UMAPs
color_palette = sns.color_palette('Set1', k)
sc.pl.umap(adata_mu, color = ['attractor', 'entropy', 'ICAM1','EPCAM', 'CLDN3', 'FHIT', 'CORIN'], vmax = 'p95', palette = color_palette, show = False)
plt.savefig(f"{plotdir}{sample_name}_umap_attractor_entropy.pdf")
plt.savefig(f"{plotdir}{sample_name}_umap_attractor_entropy.png")

sc.pl.violin(adata_mu, keys = 'entropy', groupby = 'attractor', palette = color_palette, size = 0, show = False)
plt.savefig(f"{plotdir}{sample_name}_violin_attractor_entropy.pdf")
plt.savefig(f"{plotdir}{sample_name}_violin_attractor_entropy.png")

sc.pl.violin(adata_mu, keys = ['ICAM1','EPCAM', 'CLDN3', 'FHIT', 'CORIN'], groupby = 'attractor', palette = color_palette, size = 0, show = False)
plt.savefig(f"{plotdir}{sample_name}_violin_markers.pdf")
plt.savefig(f"{plotdir}{sample_name}_violin_markers.png")

## Save Anndata object
for keys in list(adata_mu.uns['da_out'].keys()):
    if type(adata_mu.uns['da_out'][keys]).__name__ == 'double':
        adata_mu.uns['da_out'][keys] = np.asarray(adata_mu.uns['da_out'][keys])    

for keys in list(adata_mu.uns['land'].keys()):
    if type(adata_mu.uns['land'][keys]).__name__ == 'double':
        adata_mu.uns['land'][keys] = np.asarray(adata_mu.uns['land'][keys])    

del adata_mu.uns['land']['model']

adata_mu.write(f"{datadir}{sample_name}_mutrans.h5ad")