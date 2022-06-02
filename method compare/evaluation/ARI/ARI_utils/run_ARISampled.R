# Author : Nicole Lee
# Date : 29/08/2019
# Purpose: First function to be called in ARI pipeline
#          Reads in relevant dataset and grabs essential columns
#          Calls the following function 'ari_calcul_sampled'

run_ARISampled <- function(fn, this_dir, out_dir, eval_metric, methods_use){
  setwd(this_dir)
  dataset_no<-rev(strsplit(this_dir, split = "/")[[1]])[1]
  
  thisData <- read.csv(paste0(this_dir,'/',fn), head=T, row.names = 1, check.names = FALSE)
  
  # Get relevant columns
  colPCA <- grep('(0)|(1)|(2)|(3)|(4)|(5)|(6)|(7)|(8)|(9)',colnames(thisData))
  colPCA <- colPCA[1:10]
  
  
  #setwd(paste0(out_dir, '/', eval_metric, "_OP"))
  #temp<-ari_calcul_sampled(myData=thisData, cpcs=colPCA, isOptimal=TRUE,
  #                         method_use = methods_use,  
  #                         base_name=paste0(dataset_no, eval_metric, '_OP_'))
  setwd(paste0(out_dir, '/', eval_metric, "_CT"))
  temp<-ari_calcul_sampled(myData=thisData, cpcs=colPCA, isOptimal=FALSE, 
                           method_use = methods_use,  
                           base_name=paste0(dataset_no, eval_metric, '_CT_'))
  return(temp)
}
