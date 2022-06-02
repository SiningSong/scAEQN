# reference : A benchmark of batch-effect correction methods for single-cell RNA sequencing data

# clear workspace
rm(list=ls())

# source relevant functions
utils_dir <- 'E:/mymodel/methodcompare/evaluation/'
source(paste0(utils_dir,"ARI/ARI_utils/run_ARISampled.R"))
source(paste0(utils_dir,"ARI/ARI_utils/ari_calcul_sampled.R"))
source(paste0(utils_dir,"ARI/ARI_utils/conclude_ARISampled.R"))

##############################
############ Set arguments

eval_metric <- 'ARISampled'

# read in files from this_dir
this_dir<-"E:/mymodel/methodcompare/datasets/mesc_pca/"  

# send output to out_dir
out_dir<-"E:/mymodel/methodcompare/results/mesc/"

# create relevant folder within the out_dir
dir.create(paste0(out_dir, eval_metric, "_CT"),showWarnings = FALSE)
#dir.create(paste0(out_dir, eval_metric, "_OP"),showWarnings = FALSE)

##############################
############ Run 'run_ARISampled' function on relevant files

methods_use <- 'ComBat'
fn <- 'mesc_combat_pca.csv'
Rcombat<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'MNN_Correct'
fn <- 'mesc_mnn_pca.csv'
Rmnncorrect<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'limma'
fn <- 'mesc_limma_pca.csv'
Rlimma<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'Harmony'
fn <- 'mesc_harmony_pca.csv'
Rharmony<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'bbknn'
fn <- 'mesc_bbknn_pca.csv'
Rbbknn<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'quantnorm'
fn <- 'mesc_quantnorm_pca.csv'
Rquantnorm<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'Raw'
fn <- 'mesc_raw_pca.csv'
Rraw<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'scbatch'
fn <- 'mesc_scbatch_pca.csv'
Rscbatch<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'Scanorama'
fn <- 'mesc_scanorama_pca.csv'
Rscanorama<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'scVI'
fn <- 'mesc_scvi_pca.csv'
Rscvi<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)

methods_use <- 'scAEQN'
fn <- 'mesc_scAEQN_pca.csv'
RscAEQN<-run_ARISampled(fn,this_dir,out_dir,eval_metric,methods_use)
####################################################################

##############################
############ Extracting all data from all methods in dataset 

# Reads files from dir_this 
dir_this<-paste0(out_dir, eval_metric, "_CT")
#dir_this<-paste0(out_dir, eval_metric, "_OP")

# Perform the following function to produce final CSV file
wholedf<-conclude_ARISampled(dir_this, "mesc") 

save.image(paste0("mesc", "_complete.RData"))

####################################################################

##############################
############ Consolidate all raw data into one file

setwd(dir_this)
rm(list=ls())

source("E:/mymodel/methodcompare/evaluation/ARI/ARI_utils/ARI_files_consolidate.R")
 
ari_consolidate()

### END
