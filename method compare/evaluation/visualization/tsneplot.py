import pandas as pd
import anndata as ad
import scanpy as sc

##############读取数据#####################
data_path='E:/mymodel/methodcompare/datasets/brain/'
mnndata = pd.read_csv(data_path+'brain_mnn.csv', index_col=0)
limmadata = pd.read_csv(data_path+'brain_limma.csv', index_col=0)
combatdata = pd.read_csv(data_path+'brain_combat.csv', index_col=0)
quantdata = pd.read_csv(data_path+'brain_quant.csv', index_col=0)
scbatchdata = pd.read_csv(data_path+'brain_scbatch.csv', index_col=0)
bbknndata = pd.read_csv(data_path+'brain_bbknn.csv', index_col=0)
scaeqndata = pd.read_csv(data_path+'brain_scaeqn_latent.csv', index_col=0)
harmonydata = pd.read_csv(data_path+'brain_harmony.csv', index_col=0)
scanoramadata = pd.read_csv(data_path+'brain_scanorama.csv', index_col=0)
scvidata = pd.read_csv(data_path+'brain_scvi.csv', index_col=0)
adata = ad.read_h5ad('E:/mymodel/data/brain/brain.h5ad')

######################将数据转换为adata格式#############################
mnndata.columns = adata.obs.index
adata_mnn = ad.AnnData(mnndata.T, obs=adata.obs)
limmadata.columns = adata.obs.index
adata_limma = ad.AnnData(limmadata.T, obs=adata.obs)
combatdata.columns = adata.obs.index
adata_combat = ad.AnnData(combatdata.T, obs=adata.obs)
quantdata.columns = adata.obs.index
adata_quant = ad.AnnData(quantdata.T, obs=adata.obs)
scbatchdata.columns = adata.obs.index
adata_scbatch = ad.AnnData(scbatchdata.T, obs=adata.obs)
bbknndata.columns = adata.obs.index
adata_bbknn = ad.AnnData(bbknndata.T, obs=adata.obs)
adata_scaeqn = ad.AnnData(scaeqndata.T, obs=adata.obs)
harmonydata.index = adata.obs.index
adata_harmony = ad.AnnData(harmonydata, obs=adata.obs)
scanoramadata.index = adata.obs.index
adata_scanorama = ad.AnnData(scanoramadata, obs=adata.obs)
scvidata.index = adata.obs.index
adata_scvi = ad.AnnData(scvidata, obs=adata.obs)

###########################降维画聚类tsne图##############################
import matplotlib.pyplot as plt
with plt.rc_context({'figure.figsize':(6, 6)}):
    sc.tl.tsne(adata)
    fig1 = sc.pl.tsne(adata, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=50)
    sc.tl.tsne(adata_scaeqn)
    fig2 = sc.pl.tsne(adata_scaeqn, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=20)
    sc.tl.tsne(adata_mnn)
    fig3 = sc.pl.tsne(adata_mnn, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=20)
    sc.tl.tsne(adata_limma)
    fig4 = sc.pl.tsne(adata_limma, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=20)
    sc.tl.tsne(adata_combat)
    fig5 = sc.pl.tsne(adata_combat, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=20)
    sc.tl.tsne(adata_quant)
    fig6 = sc.pl.tsne(adata_quant, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=20)
    sc.tl.tsne(adata_scbatch)
    fig7 = sc.pl.tsne(adata_scbatch, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=20)
    sc.tl.tsne(adata_bbknn)
    fig8 = sc.pl.tsne(adata_bbknn, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True,size=20)
    sc.tl.tsne(adata_harmony)
    fig9 = sc.pl.tsne(adata_harmony, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=20)
    sc.tl.tsne(adata_scanorama)
    fig10 = sc.pl.tsne(adata_scanorama, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True, size=20)
    sc.tl.tsne(adata_scvi)
    fig11 = sc.pl.tsne(adata_scvi, color=['cell type', 'batch'], legend_fontsize='xx-small', return_fig=True,
                       size=20)

savedir = 'E:/mymodel/figure/tsne mouse panc/'
fig1.savefig(savedir+'uncorrected.jpg', dpi=300)
fig2.savefig(savedir+'scAEQN.jpg', dpi=600)
fig3.savefig(savedir+'mnn.jpg', dpi=600)
fig4.savefig(savedir+'limma.jpg', dpi=600)
fig5.savefig(savedir+'combat.jpg', dpi=600)
fig6.savefig(savedir+'quant.jpg', dpi=600)
fig7.savefig(savedir+'scbatch.jpg', dpi=600)
fig8.savefig(savedir+'bbknn.jpg', dpi=600)
fig9.savefig(savedir+'harmony.jpg', dpi=600)
fig10.savefig(savedir+'scanorama.jpg', dpi=600)
fig11.savefig(savedir+'scvi.jpg',dpi=600)



##################差异表达基因鉴定########################
scaeqndata = pd.read_csv(data_path+'mousepanc_scaeqn_output.csv',index_col=0)
adata_scaeqn = ad.AnnData(scaeqndata.T, obs=adata.obs)
sc.pp.highly_variable_genes(adata, n_top_genes=20)

sc.tl.rank_genes_groups(adata, 'cell type', n_genes=20, method='wilcoxon')


