import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import time
from datetime import timedelta
import scanpy as sc
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
print(sc.logging.print_versions())
import os
dirname = os.getcwd()
print(dirname)

from sklearn.metrics import silhouette_score


def silhouette_coeff_ASW(adata, method_use='raw', save_dir='', save_fn='', percent_extract=0.8):
    random.seed(0)
    asw_fscore = []
    asw_bn = []
    asw_bn_sub = []
    asw_ctn = []
    iters = []
    for i in range(20):
        iters.append('iteration_' + str(i + 1))
        rand_cidx = np.random.choice(adata.obs_names, size=int(len(adata.obs_names) * percent_extract), replace=False)
        print('nb extracted cells: ', len(rand_cidx))
        adata_ext = adata[rand_cidx, :]
        asw_batch = silhouette_score(adata_ext.X, adata_ext.obs['batch'])
        asw_celltype = silhouette_score(adata_ext.X, adata_ext.obs['cell_type'])
        min_val = -1
        max_val = 1
        asw_batch_norm = (asw_batch - min_val) / (max_val - min_val)
        asw_celltype_norm = (asw_celltype - min_val) / (max_val - min_val)

        fscoreASW = (2 * (1 - asw_batch_norm) * (asw_celltype_norm)) / (1 - asw_batch_norm + asw_celltype_norm)
        asw_fscore.append(fscoreASW)
        asw_bn.append(asw_batch_norm)
        asw_bn_sub.append(1 - asw_batch_norm)
        asw_ctn.append(asw_celltype_norm)

    #     iters.append('median_value')
    #     asw_fscore.append(np.round(np.median(fscoreASW),3))
    #     asw_bn.append(np.round(np.median(asw_batch_norm),3))
    #     asw_bn_sub.append(np.round(1 - np.median(asw_batch_norm),3))
    #     asw_ctn.append(np.round(np.median(asw_celltype_norm),3))
    asw_bn_mean = np.mean(asw_bn)
    asw_bn_sub_mean = np.mean(asw_bn_sub)
    asw_ctn_mean = np.mean(asw_ctn)
    asw_fscore_mean = np.mean(asw_fscore)
    asw_mean = [method_use, asw_bn_mean, asw_bn_sub_mean, asw_ctn_mean,asw_fscore_mean]
    df = pd.DataFrame({'asw_batch_norm': asw_bn, 'asw_batch_norm_sub': asw_bn_sub,
                       'asw_celltype_norm': asw_ctn, 'fscore': asw_fscore,
                       'method_use': np.repeat(method_use, len(asw_fscore))})
    df.to_csv(os.path.join(save_dir + save_fn + '.csv'))
    print('Save output of pca in: ', save_dir)
    print(df.values.shape)
    print(df.keys())
    return df, asw_mean


def createAnnData(data_dir, myDatafn):
    myData = pd.read_csv(os.path.join(data_dir, myDatafn), header=0, index_col=0)
    bex = ['batch', 'Batch', 'Batchlb', 'batchlb', 'BATCH']
    ib = np.isin(myData.keys(), bex)
    cex = ['celltype', 'CellType', 'cell_type', 'Cell_Type', 'ct']
    ict = np.isin(myData.keys(), cex)
    adata = sc.AnnData(myData.values[:, 0:9])
    adata.obs_names = myData.index
    adata.obs['batch'] = myData.values[:, np.where(ib)[0][0]]  # factor function in R
    adata.obs['cell_type'] = myData.values[:, np.where(ict)[0][0]]
    print(adata)
    return adata

data_dir = 'E:/mymodel/methodcompare/datasets/mesc_pca/'
save_dir = 'E:/mymodel/methodcompare/results/mesc/'
print(save_dir)
if not os.path.exists(save_dir+'asw_output/'): os.makedirs(os.path.join(save_dir,'asw_output/'))
fls = [ f for f in os.listdir(data_dir) if f.endswith("_pca.csv") & os.path.isfile(os.path.join(data_dir,f)) ]
fls

import random
final_ls = []
for f in fls:
    method_use = f[5:-8]
    print('Extract asw for ', method_use)
    save_fn = 'asw_' + method_use
    print(method_use)
    adata = createAnnData(data_dir, f)
    asw_val, asw_mean = silhouette_coeff_ASW(
        adata, method_use, os.path.join(save_dir, 'asw_output/'), save_fn, percent_extract=0.8)
    final_ls.append(asw_mean)


result = pd.DataFrame(final_ls, columns=['method','asw_bn_mean','asw_bn_sum_mean','asw_cln_mean','asw_fscore_mean'])
result.index = result['method']
result =  result.drop('method',axis=1)
result.to_csv(os.path.join(save_dir, 'asw_output/', 'final_asw_summary.csv'))