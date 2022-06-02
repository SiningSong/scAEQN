import scvi
import pandas as pd
import scanpy as sc
import anndata as ad
import numpy as np
from numpy import unique
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans
from sklearn import metrics
import time

adata = sc.read_h5ad('/data/brain/brain.h5ad')
print(adata)

# BBKNN方法
start1 = time.time()
sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key='batch')
data_bbknn = pd.DataFrame(adata.X.T, index=adata.var['genes'], columns=adata.obs_names)
data_bbknn.to_csv('E:/mymodel/methodcompare/brain/brain_bbknn.csv')

sc.tl.umap(adata)
sc.pl.umap(adata, color=['batch', 'cell type'])
end1 = time.time()
tbbknn = end1-start1
# kmeans聚类，并计算ARI
ARI_bbknn = 0
for i in range(10):
    cluster = KMeans(n_clusters=6)
    cluster.fit(adata.X)
    cluster_result = cluster.labels_
    unique(cluster_result)
    ARI_bbknn = ARI_bbknn+metrics.adjusted_rand_score(adata.obs['cell type'],cluster_result)

ARI_bbknn = ARI_bbknn/10
print(ARI_bbknn)
print(tbbknn)

###scvi
start2 = time.time()
scvi.data.setup_anndata(adata, batch_key="batch")
vae = scvi.model.SCVI(adata)
vae.train()
adata.obsm["X_scVI"] = vae.get_latent_representation()
adata.obsm["X_normalized_scVI"] = vae.get_normalized_expression()
scvi.data.view_anndata_setup(adata)
end2 = time.time()
tscvi = end2 - start2

ARI_scvi = 0
for i in range(10):
    cluster_scvi = KMeans(n_clusters=6)
    cluster_scvi.fit(adata.obsm['X_scVI'])
    cluster_result_scvi = cluster_scvi.labels_
    ARI_scvi = ARI_scvi + metrics.adjusted_rand_score(adata.obs['cell type'], cluster_result_scvi)

ARI_scvi = ARI_scvi/10    #scvi的聚类指标
print(ARI_scvi)
print(tscvi)
data_scvi = pd.DataFrame(adata.obsm['X_scVI'], index=adata.obs_names)
data_scvi.to_csv('E:/mymodel/methodcompare/brain/brain_scvi.csv')

# scanorama
start3 = time.time()
adata_sca = sc.read_h5ad('/data/brain/brain.h5ad')
sc.tl.pca(adata_sca)
sc.external.pp.scanorama_integrate(adata_sca, key='batch')
sca_adata = adata_sca.obsm['X_scanorama']
end3 = time.time()
time_sca = end3-start3

ARI_sca = 0
for i in range(10):
    cluster_sca = KMeans(n_clusters=6)
    cluster_sca.fit(sca_adata)
    cluster_result_sca = cluster_sca.labels_
    unique(cluster_result_sca)
    ARI_sca = ARI_sca + metrics.adjusted_rand_score(adata_sca.obs['cell type'], cluster_result_sca)

ARI_sca = ARI_sca/10
print(ARI_sca)
print(time_sca)

data_scanorama = pd.DataFrame(sca_adata, index=adata.obs_names)
data_scanorama.to_csv('E:/mymodel/methodcompare/brain/brain_scanorama.csv')

# harmony
start4 = time.time()
sc.tl.pca(adata)
sc.external.pp.harmony_integrate(adata, key='batch')
harmonymod = adata.obsm['X_pca_harmony']
end4 = time.time()
time_har = end4-start4

ARI_har = 0
for i in range(10):
    cluster_har = KMeans(n_clusters=6)
    cluster_har.fit(harmonymod)
    cluster_result_har = cluster_har.labels_
    unique(cluster_result_har)
    ARI_har = ARI_har + metrics.adjusted_rand_score(adata.obs['cell type'], cluster_result_har)

ARI_har = ARI_har/10
print(ARI_har)
print(time_har)
data_harmony = pd.DataFrame(harmonymod, index=adata.obs_names)
data_harmony.to_csv('E:/mymodel/methodcompare/brain/brain_harmony.csv')
