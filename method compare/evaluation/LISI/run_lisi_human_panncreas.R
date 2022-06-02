
library(ggplot2)
library(lisi)
library(ggpubr)

rm(list=ls())
this_dir <- 'E:/mymodel/methodcompare/datasets/human_pancreas_pca/'

out_dir<-"E:/mymodel/methodcompare/results/human pancreas/"

dataset_use = 'humanpanc'

setwd(this_dir)
eval_metric <- 'LISI_v2/' 

utils_dir <- 'E:/mymodel/methodcompare/evaluation/'
source(paste0(utils_dir,'LISI/lisi_utils.R'))

plx = 40



# Get output of LISI
methods_use <- c('Raw','mnn','limma',
                 'Harmony','combat','quantnorm',
                 'scbatch','bbknn','scVI',
                 'Scanorama','scAEQN')
method_dir <- c('raw','mnn','limma',
                'harmony','combat','quantnorm',
                'scbatch','bbknn','scvi',
                'scanorama','scAEQN')

fn_ls <- c()
for(i in rep(1:length(method_dir),1)){
  fn_ls <- c(fn_ls, paste0(dataset_use,'_', method_dir[i],'_','pca.csv'))
}
print(fn_ls)
dir.create(paste0(out_dir, eval_metric), showWarnings = F)

length(methods_use) == length(fn_ls)
for(i in rep(1:length(methods_use), 1)){
  print(methods_use[i])
  run_LISI_final(fn_ls[i], this_dir, out_dir, eval_metric, methods_use[i], plx)
}


#####################################
# iLISI batch
# Combine all results together
# 13 methods (except BBKNN) + raw data
#####################################


fn <- paste0(dataset_use,'_raw_pca.csv')
meta_ls <- get_celltype_common(this_dir, fn)
length(meta_ls$cells_common)
meta_ls$ct_common

get_cells_integration_iLISI_v2(dataset_use, meta_ls, out_dir, plx, eval_metric)


get_celltype_mixing_cLISI(dataset_use, out_dir, plx, eval_metric)


summary_LISI(out_dir, plottitle=paste0('LISI - ',dataset_use), plx, eval_metric)
