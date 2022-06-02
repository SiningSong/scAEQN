import pandas as pd
import anndata as ad
import scanpy as sc

##############读取数据#####################
data_path= 'E:/mymodel/methodcompare/datasets/panc2tech/'
save_dir = 'E:/mymodel/methodcompare/datasets/panc2tech_pca/'

batch = pd.read_csv(data_path+'batch.csv', index_col=0)
celltype = pd.read_csv(data_path+'celltype.csv', index_col=0)

mnndata = pd.read_csv(data_path+'panc2tech_mnn.csv', index_col=0)
limmadata = pd.read_csv(data_path+'panc2tech_limma.csv', index_col=0)
combatdata = pd.read_csv(data_path+'panc2tech_combat.csv', index_col=0)
quantdata = pd.read_csv(data_path+'panc2tech_quant.csv', index_col=0)
#scbatchdata = pd.read_csv(data_path+'panc2tech_scbatch.csv', index_col=0)
bbknndata = pd.read_csv(data_path+'panc2tech_bbknn.csv', index_col=0)
scaeqndata = pd.read_csv(data_path+'panc2tech_scaeqn_latent.csv', index_col=0)
harmonydata = pd.read_csv(data_path+'panc2tech_harmony.csv', index_col=0)
scanoramadata = pd.read_csv(data_path+'panc2tech_scanorama.csv', index_col=0)
scvidata = pd.read_csv(data_path+'panc2tech_scvi.csv', index_col=0)

adata = ad.read_h5ad('E:/mymodel/data/panc 2 tech/panc2tech.h5ad')

######################将数据转换为adata格式#############################
mnndata.columns = adata.obs.index
adata_mnn = ad.AnnData(mnndata.T,obs=adata.obs)
limmadata.columns = adata.obs.index
adata_limma = ad.AnnData(limmadata.T,obs=adata.obs)
combatdata.columns = adata.obs.index
adata_combat = ad.AnnData(combatdata.T,obs=adata.obs)
quantdata.columns = adata.obs.index
adata_quant = ad.AnnData(quantdata.T,obs=adata.obs)
#scbatchdata.columns = adata.obs.index
#adata_scbatch = ad.AnnData(scbatchdata.T,obs=adata.obs)
bbknndata.columns = adata.obs.index
adata_bbknn = ad.AnnData(bbknndata.T,obs=adata.obs)
adata_scaeqn = ad.AnnData(scaeqndata.T,obs=adata.obs)
harmonydata.index = adata.obs.index
adata_harmony = ad.AnnData(harmonydata, obs=adata.obs)
scanoramadata.index = adata.obs.index
adata_scanorama = ad.AnnData(scanoramadata, obs=adata.obs)
scvidata.index = adata.obs.index
adata_scvi = ad.AnnData(scvidata, obs=adata.obs)

batch.index = adata.obs.index
celltype.index = adata.obs.index

# 将校正后的数据降维到10PCs保存降维后的数据
sc.tl.pca(adata, 10)  # uncorrected
data_pca = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs.index)
sc.tl.pca(adata_bbknn, 10)   # bbknn
data_bbknn_pca = pd.DataFrame(adata_bbknn.obsm['X_pca'], index=adata.obs.index)
sc.tl.pca(adata_combat, 10)   # combat
data_combat_pca = pd.DataFrame(adata_combat.obsm['X_pca'], index=adata.obs.index)
sc.tl.pca(adata_harmony, 10)   # harmony
data_harmony_pca = pd.DataFrame(adata_harmony.obsm['X_pca'], index=adata.obs.index)
sc.tl.pca(adata_limma, 10)     # limma
data_limma_pca = pd.DataFrame(adata_limma.obsm['X_pca'], index=adata.obs.index)
sc.tl.pca(adata_mnn, 10)    # mnn
data_mnn_pca = pd.DataFrame(adata_mnn.obsm['X_pca'], index=adata.obs.index)
sc.tl.pca(adata_quant, 10)   # quantnorm
data_quant_pca = pd.DataFrame(adata_quant.obsm['X_pca'], index=adata.obs.index)
sc.tl.pca(adata_scaeqn, 10)  # scaeqn
data_scAEQN_pca = pd.DataFrame(adata_scaeqn.obsm['X_pca'], index=adata.obs.index)
sc.tl.pca(adata_scanorama, 10)  # scanorama
data_scanorama_pca = pd.DataFrame(adata_scanorama.obsm['X_pca'], index=adata.obs.index)
#sc.tl.pca(adata_scbatch, 10)  # scbatch
#data_scbatch_pca = pd.DataFrame(adata_scbatch.obsm['X_pca'], index=adata.obs.index)
data_scvi_pca = pd.DataFrame(adata_scvi.X, index=adata.obs.index)

data_pca['batch'] = batch['batch']
data_pca['celltype'] = celltype['celltype']
data_bbknn_pca['batch'] = batch['batch']
data_bbknn_pca['celltype'] = celltype['celltype']
data_combat_pca['batch'] = batch['batch']
data_combat_pca['celltype'] = celltype['celltype']
data_harmony_pca['batch'] = batch['batch']
data_harmony_pca['celltype'] = celltype['celltype']
data_limma_pca['batch'] = batch['batch']
data_limma_pca['celltype'] = celltype['celltype']
data_mnn_pca['batch'] = batch['batch']
data_mnn_pca['celltype'] = celltype['celltype']
data_quant_pca['batch'] = batch['batch']
data_quant_pca['celltype'] = celltype['celltype']
data_scAEQN_pca['batch'] = batch['batch']
data_scAEQN_pca['celltype'] = celltype['celltype']
data_scanorama_pca['batch'] = batch['batch']
data_scanorama_pca['celltype'] = celltype['celltype']
#data_scbatch_pca['batch'] = batch['batch']
#data_scbatch_pca['celltype'] = celltype['celltype']
data_scvi_pca['batch'] = batch['batch']
data_scvi_pca['celltype'] = celltype['celltype']

data_pca.to_csv(save_dir+'panc2tech_raw_pca.csv')
data_bbknn_pca.to_csv(save_dir+'panc2tech_bbknn_pca.csv')
data_combat_pca.to_csv(save_dir+'panc2tech_combat_pca.csv')
data_harmony_pca.to_csv(save_dir+'panc2tech_harmony_pca.csv')
data_limma_pca.to_csv(save_dir+'panc2tech_limma_pca.csv')
data_mnn_pca.to_csv(save_dir+'panc2tech_mnn_pca.csv')
data_quant_pca.to_csv(save_dir+'panc2tech_quantnorm_pca.csv')
data_scAEQN_pca.to_csv(save_dir+'panc2tech_scAEQN_pca.csv')
data_scanorama_pca.to_csv(save_dir+'panc2tech_scanorama_pca.csv')
#data_scbatch_pca.to_csv(save_dir+'panc2tech_scbatch_pca.csv')
data_scvi_pca.to_csv(save_dir+'panc2tech_scvi_pca.csv')


