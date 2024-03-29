normalize_values <- function(myData, colns, min_val=1, max_val=2){
  med_val <- c()
  med_val_norm <- c()
  methods_use <- c()
  for(c in colnames(myData)){
    methods_use <- c(methods_use, c)
    mi <- median(myData[,c])
    med_val <- c(med_val, mi)
    med_val_norm <- c(med_val_norm, (mi - min_val) / (max_val - min_val))
  }
  myDataNorm <- data.frame('X1'=methods_use, 'X2'=med_val, 'X3'=med_val_norm)
  colnames(myDataNorm) <- colns
  return(myDataNorm)
}


get_celltype_common <- function(data_dir, fn){
  myData <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
  colnames(myData)[grep('[cC]ell_?[tT]ype',colnames(myData))] <- 'cell_type'
  colnames(myData)[grep('([bB]atch)|(BATCH)',colnames(myData))] <- 'batch'
  print(dim(myData))
  batches <- unique(myData$batch)
  celltypels <- unique(myData$cell_type)
  print(celltypels)
  if(length(batches)>1){
    ctls <- list()
    count <- 0
    for (b in batches){
      count <- count + 1
      ct <- unique(myData[which(myData$batch==b),'cell_type'])
      ctls[[count]] <- ct
      print("Batch")
      print(b)
      print(ct)
    }
    for(i in rep(1:length(ctls))){
      if(i==2){
        ct_common <- intersect(ctls[[i-1]], ctls[[i]])    
      }
      if(i>2){   #more than 2 batches
        ct_common <- intersect(ct_common, ctls[[i]])    
      }
    }
    ct_common <- unique(ct_common)
    print(ct_common)
    cells_common <- rownames(myData)[which(myData$cell_type %in% ct_common) ]
    return(list('ct_common'=ct_common, 'cells_common'=cells_common, 'batches'=batches, 'celltypels'=celltypels))
  } else{
    return(NULL)
  } 
  
}


run_LISI_final <- function(fn, data_dir, save_dir, eval_metric, methods_use, plx=40){
  
  myPCA <- read.csv(paste0(data_dir, fn), head=T, row.names = 1, check.names = FALSE)
  #myPCA <- ezTools::fast_read_table(paste0(data_dir, fn), row.names = T, col.names = T, sep = ",", check.names = FALSE)
  cpcs <- grep('(0)|(1)|(2)|(3)|(4)|(5)|(6)|(7)|(8)|(9)',colnames(myPCA))
  cpcs <- cpcs[1:10]
  lisi_embeddings <- myPCA[,cpcs]
  
  #colnames(myPCA)[grep('[cC]ell_?[tT]ype',colnames(myPCA))] <- 'cell_type'
  #colnames(myPCA)[grep('([bB]atch)|(BATCH)|(batchlb)',colnames(myPCA))] <- 'batch'
  
  lisi_meta_data <- subset(myPCA, select=c('batch','celltype'))
  
  lisi_label = c('batch', 'celltype')
  
  lisi_res <- lisi::compute_lisi(lisi_embeddings, lisi_meta_data, lisi_label,perplexity = plx)
  lisi_res$cell <- rownames(lisi_embeddings)
  
  lisi_batch <- subset(lisi_res,select=c('batch','cell'))
  lisi_celltype <- subset(lisi_res,select=c('celltype','cell'))

  write.table(lisi_batch, paste0(save_dir,eval_metric,'lisi_batch','_',plx,'_',methods_use,'.txt'), quote=F, sep=',', row.names=T, col.names=NA)
  write.table(lisi_celltype,paste0(save_dir,eval_metric,'lisi_celltype','_',plx,'_',methods_use,'.txt'), quote=F, sep=',', row.names=T, col.names=NA)
  #dir.create(paste0(this_dir,eval_metric,methods_use,'/'), showWarnings = F)
  #ezTools::fast_save_table(lisi_batch, paste0(this_dir,eval_metric,methods_use,'/'), paste0('lisi_batch_',plx,'.txt'),
  #                         row.names = T, col.names = T, quote = F, sep = "\t")
  
  #ezTools::fast_save_table(lisi_celltype, paste0(this_dir,eval_metric,methods_use,'/'), paste0('lisi_celltype_',plx,'.txt'),
  #                         row.names = T, col.names = T, quote = F, sep = "\t")
}


summary_LISI <- function(this_dir, plottitle='LISI - dataset', plx=40, eval_metric='LISI/',ht=400, wd=400){
  iLISI_df <- read.csv(paste0(this_dir, eval_metric,"result/",plx,"/iLISI_summary.csv"), head=T, check.names = F)
  colns <- c('methods_use', 'iLISI_median', 'iLISI_median_norm')
  # median_iLISI <- normalize_values(iLISI_df, colns, min_val=1, max_val=length(meta_ls$batches))
  median_iLISI <- normalize_values(iLISI_df, colns, min_val=min(iLISI_df), max_val=max(iLISI_df))
  mini <- min(median_iLISI$iLISI_median_norm)
  maxi <- max(median_iLISI$iLISI_median_norm)
  median_iLISI$iLISI_median_norm2 <- (median_iLISI$iLISI_median_norm - mini) / (maxi - mini)
  
  
  cLISI_df <- read.csv(paste0(this_dir,eval_metric,"result/",plx,"/cLISI_summary.csv"), head=T, check.names = F)
  colns <- c('methods_use', 'cLISI_median', 'cLISI_median_norm')
  # median_cLISI <- normalize_values(cLISI_df, colns, min_val=1, max_val=length(meta_ls$celltypels))
  median_cLISI <- normalize_values(cLISI_df, colns, min_val=min(cLISI_df), max_val=max(cLISI_df))
  median_cLISI$cLISI_median_norm_sub <- 1 - median_cLISI$cLISI_median_norm
  minc <- min(median_cLISI$cLISI_median_norm_sub)
  maxc <- max(median_cLISI$cLISI_median_norm_sub)
  median_cLISI$cLISI_median_norm_sub2 <- (median_cLISI$cLISI_median_norm_sub - minc) / (maxc - minc)

  
  
  final_df = merge(median_iLISI, median_cLISI, by="methods_use")
  final_df$sum_normXY <- final_df$iLISI_median_norm2 + final_df$cLISI_median_norm_sub2
  final_df$fscore <- (2 * final_df$iLISI_median_norm2 * final_df$cLISI_median_norm_sub2)/
                     (final_df$iLISI_median_norm2 + final_df$cLISI_median_norm_sub2)
  final_df <- final_df[order(final_df$fscore, decreasing = T),]
  
  # write.table(final_df,paste0(this_dir, eval_metric, "result/", plx, '/', 'summary_median_', plx, '.txt'), 
  #             quote=F, sep='\t', row.names=T, col.names=T)
  write.csv(final_df,paste0(this_dir, eval_metric, "result/", plx, '/', 'LISI_summary_median_', plx, '.csv'), 
              quote=F, row.names=F)
  
  
  # plot final LISI
  plot_final_LISI(final_df, plottitle, this_dir, plx, ht, wd, eval_metric, 
                  xstring = 'cLISI_median_norm_sub', ystring = 'iLISI_median_norm', plottype = 'methods_use')
  
}

#' Boxplot of LISI index
#' include library: library(ggpubr)
#' @param data : list of lisi index dataset
#' @param vect_names_method : name of each methods in the same order as data (label of y-axis) 
#' @param base_name : path directory
#' @param title : label of x-axis
#' @param toptitle : title
#' @return a horizontal boxplot to plot LISI values from different methods
#' @examples
LISI_boxplot_fun <- function(data,vect_names_method, base_name, title='iLISI', toptitle='Batch Mixing', save_results=F){
  
  library(ggpubr)
  
  # combine the LISI index from all the method in one dataset
  if(length(data)>1){
    for (i in seq_along(data)){
      colnames(data[[i]]) <- c(vect_names_method[i],'cell')
    }
    data <- Reduce(function(x, y) merge(x, y, by="cell"), data)
    data <- data[,-1]
  } else {
    data <- as.data.frame(data[[1]][,1])
    # colnames(data) <- vect_names_method
  }
  
  cn <- c()
  for (i in seq_along(data)){
    cn <- c(cn, vect_names_method[i])
  }
  totaldf <- as.data.frame(data)
  colnames(totaldf) <- cn
  if(save_results){
    write.csv(totaldf, file=paste0(base_name, title,"_summary.csv"), row.names=FALSE)
    med_val <- c()
    for (c in cn){
      med_val <- c(med_val, round(median(totaldf[,c]), digits = 3))
    }
    median_df <- data.frame('methods'=cn, 'median_value'=med_val)
    write.csv(median_df, file=paste0(base_name, title,"_median_val.csv"), row.names=FALSE)
    print('Results saved')
  }
  
  # melt to plot 
  plotdata <- tidyr::gather(data)
  
  # put the y-axis in order 
  plotdata$key <- factor(plotdata$key,levels = rev(vect_names_method), order=T)
  # write.csv(plotdata, file = paste0(base_name,'summary_',title,'.csv'), row.names = T, col.names = T)
  # print('Save the output summary in the folder: ')
  # print(base_name)
  # plot
  
  p <- ggplot(plotdata, aes(key, value))
  p <- p + geom_boxplot(color='black',outlier.size = 1) + coord_flip() # outlier.size=-1 remove outliers dots
  p <- p + labs(title=toptitle,x=" ", y = title)
  p <- p + theme(axis.title.y = element_blank(), 
                 axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
                 axis.title.x = element_text(size=15),
                 axis.text.x = element_text(size=11,colour = 'black'),
                 axis.text.y = element_text(size=11,colour = 'black'),
                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                 panel.border = element_blank(), 
                 plot.title = element_text(hjust = 0.5)) 
  # p <- p + theme_gray(base_size = 14)
  # axis.line.y=element_blank(), axis.ticks.y=element_blank()
  
  # save
  png(paste0(base_name, "boxplot_",title,".png"),height = 800, width=800, res = 2*72)
  print('Boxplot is saved')
  return(p)
}

scatter_plot_func <- function(plot_df, xstring, ystring, plottype, plottitle='Dataset', colorcode, xlabel='ARI(Cell type)', ylabel='1-ARI(batch)') {
  library(ggrepel)
  
  plot_df <- plot_df[order(as.character(plot_df[,plottype])),]
  
  plotobj <- ggplot(plot_df, aes_string(x = xstring, y = ystring, colour = plottype)) + 
    geom_point(shape=19, alpha = 1, size = 4) + theme_bw() #alpha = 0.6
  
  plotobj <- plotobj + geom_text_repel(aes_string(label = plottype),
                                       nudge_y      = 0.02,
                                       angle        = 0,
                                       direction    = "y",
                                       vjust        = 0,
                                       segment.size = 0.2)
  #ggplot_build(plotobj)$data
  
  # xmax = max(plot_df[[xstring]])
  # xmin = min(plot_df[[xstring]])
  # ymax = max(plot_df[[ystring]])
  # ymin = min(plot_df[[ystring]])
  # xmax = ceiling(xmax)
  # xmin = floor(xmin)
  # ymax = ceiling(ymax)
  # ymin = floor(ymin)
  # 
  if(plottitle==''){
    plotobj <- plotobj + labs(x=xlabel,y=ylabel) 
  } else {
    plotobj <- plotobj + labs(x=xlabel,y=ylabel,title=plottitle) 
  }
  
  plotobj <- plotobj + theme(legend.title = element_blank(), 
                             plot.title = element_text(color="black", size=18, hjust = 0.5),
                             legend.position = "none", 
                             axis.line = element_line(color = "black", size = 0.5, linetype = "solid"), 
                             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                             panel.border = element_blank()
                             # axis.text = element_blank()
                             #axis.title = element_blank(),
                             # axis.ticks = element_blank()
  ) + scale_fill_manual(values = colorcode)
  
  # + scale_color_brewer(palette=marker_palette)
  
  return(plotobj)
}


plot_final_LISI <- function(myData, plottitle='LISI', save_dir='', plx = 40, ht = 400, wd = 400, eval_metric='LISI/', 
                            xstring = 'cLISI_median_norm_sub', ystring = 'iLISI_median_norm', plottype = 'methods_use'){
  
  library(RColorBrewer)
  colourM = length(unique(myData[,plottype]))
  getPalette = colorRampPalette(brewer.pal(8, "Accent"))
  colorcode_M <- getPalette(colourM)
  lisi_plt <- scatter_plot_func(myData, xstring, ystring, 
                                plottype, paste0(plottitle,' plx ',plx), 
                                colorcode_M, 
                                xlabel='1 - cLISI cell type', 
                                ylabel='iLISI batch')
  
  
  print(lisi_plt)
  save(lisi_plt, file = paste0(save_dir, eval_metric, "result/", plx, "/mylisi_scatter.rda"))
  
  png(paste(save_dir, eval_metric,"result/",plx,"/LISI_batch_celltype_",plx,".png",sep=""),height = 2*ht, width=2*wd, res = 2*72)
  print(lisi_plt)
  dev.off()
  print('Save the output in: ')
  print(paste0(save_dir, eval_metric, "result/", plx))
}

get_cells_integration_iLISI_v2 <- function(dataset_use, meta_ls, this_dir='', plx = 40, eval_metric = 'LISI/'){
  
  mnn_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','mnn','.txt'), head=T, row.names = 1, sep=',',check.names = FALSE)
  quant_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','quantnorm','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  harmony_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Harmony','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  bbknn_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','bbknn','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scbatch_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scbatch','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scanorama_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Scanorama','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scvi_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scvi','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')  
  raw_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Raw','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')  
  combat_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','combat','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  limma_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','limma','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  scAEQN_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scAEQN','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  
  
  raw_df <- raw_df[meta_ls$cells_common,]
  mnn_df <- mnn_df[meta_ls$cells_common,]
  quant_df <- quant_df[meta_ls$cells_common,]
  harmony_df <- harmony_df[meta_ls$cells_common,]
  scbatch_df <- scbatch_df[meta_ls$cells_common,]
  scvi_df <- scvi_df[meta_ls$cells_common,]
  limma_df <- limma_df[meta_ls$cells_common,]
  combat_df <- combat_df[meta_ls$cells_common,]
  bbknn_df <- bbknn_df[meta_ls$cells_common,]
  scanorama_df <- scanorama_df[meta_ls$cells_common,]
  scAEQN_df <- scAEQN_df[meta_ls$cells_common,]
  
  dir.create(paste0(this_dir, eval_metric,"result/"), showWarnings = F)
  dir.create(paste0(this_dir, eval_metric,"result/",plx,"/"), showWarnings = F)
  methods_use <- c('raw','mnn','limma','harmony','combat','quantnorm','scbatch',
                   'scvi','scAEQN','scanorama','bbknn')
  piLISI <- LISI_boxplot_fun(data = list(raw_df, mnn_df, limma_df,combat_df,
                                         harmony_df, quant_df, scbatch_df,
                                         bbknn_df,scvi_df,  scanorama_df,  
                                         scAEQN_df),
                             vect_names_method = methods_use, 
                             title='iLISI', toptitle=paste0('iLISI Batch Mixing plx ',plx),
                             base_name=paste0(this_dir, eval_metric, "result/",plx,"/"), save_results=T)
  
  piLISI
  save(piLISI, file = paste0(this_dir, eval_metric, "result/",plx, "/piLISI_",dataset_use,"_plots_",plx,".rda"))
}



get_celltype_mixing_cLISI <- function(dataset_use, this_dir='', plx = 40, eval_metric = 'LISI/'){
  
  mnn_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','mnn','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  quant_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','quantnorm','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  harmony_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Harmony','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  bbknn_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','bbknn','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scbatch_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scbatch','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scanorama_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Scanorama','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scvi_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scvi','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')  
  raw_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Raw','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')  
  combat_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','combat','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  limma_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','limma','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  scAEQN_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scAEQN','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  
  dir.create(paste0(this_dir, eval_metric,"result/"), showWarnings = F)
  dir.create(paste0(this_dir, eval_metric,"result/",plx,"/"), showWarnings = F)
  methods_use <- c('raw','mnn','limma','harmony','combat','quantnorm','scbatch',
                   'scvi','scAEQN','scanorama','bbknn')
  pcLISI <- LISI_boxplot_fun(data =  list(raw_df, mnn_df, limma_df,combat_df,
                                          harmony_df, quant_df, scbatch_df,
                                          bbknn_df,scvi_df,  scanorama_df,  
                                          scAEQN_df),
                             vect_names_method = methods_use, 
                             title='cLISI', toptitle=paste0('cLISI Cell Type Mixing plx ',plx),
                             base_name=paste0(this_dir, eval_metric, "result/",plx,"/"), save_results=T)
  
  pcLISI
  
  save(pcLISI, file = paste0(this_dir, eval_metric, "result/",plx, "/pcLISI_",dataset_use,"_plots_",plx,".rda"))
}

get_cells_integration_iLISI_v2_panc2tech <- function(dataset_use, meta_ls, this_dir='', plx = 40, eval_metric = 'LISI/'){
  
  mnn_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','mnn','.txt'), head=T, row.names = 1, sep=',',check.names = FALSE)
  quant_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','quantnorm','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  harmony_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Harmony','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  bbknn_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','bbknn','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  #scbatch_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scbatch','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scanorama_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Scanorama','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scvi_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scvi','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')  
  raw_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Raw','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')  
  combat_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','combat','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  limma_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','limma','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  scAEQN_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scAEQN','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  
  
  raw_df <- raw_df[meta_ls$cells_common,]
  mnn_df <- mnn_df[meta_ls$cells_common,]
  quant_df <- quant_df[meta_ls$cells_common,]
  harmony_df <- harmony_df[meta_ls$cells_common,]
  #scbatch_df <- scbatch_df[meta_ls$cells_common,]
  scvi_df <- scvi_df[meta_ls$cells_common,]
  limma_df <- limma_df[meta_ls$cells_common,]
  combat_df <- combat_df[meta_ls$cells_common,]
  bbknn_df <- bbknn_df[meta_ls$cells_common,]
  scanorama_df <- scanorama_df[meta_ls$cells_common,]
  scAEQN_df <- scAEQN_df[meta_ls$cells_common,]
  
  dir.create(paste0(this_dir, eval_metric,"result/"), showWarnings = F)
  dir.create(paste0(this_dir, eval_metric,"result/",plx,"/"), showWarnings = F)
  methods_use <- c('raw','mnn','limma','harmony','combat','quantnorm',
                   'scvi','scAEQN','scanorama','bbknn')
  piLISI <- LISI_boxplot_fun(data = list(raw_df, mnn_df, limma_df,combat_df,
                                         harmony_df, quant_df, 
                                         bbknn_df,scvi_df,  scanorama_df,  
                                         scAEQN_df),
                             vect_names_method = methods_use, 
                             title='iLISI', toptitle=paste0('iLISI Batch Mixing plx ',plx),
                             base_name=paste0(this_dir, eval_metric, "result/",plx,"/"), save_results=T)
  
  piLISI
  save(piLISI, file = paste0(this_dir, eval_metric, "result/",plx, "/piLISI_",dataset_use,"_plots_",plx,".rda"))
}



get_celltype_mixing_cLISI_panc2tech <- function(dataset_use, this_dir='', plx = 40, eval_metric = 'LISI/'){
  
  mnn_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','mnn','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  quant_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','quantnorm','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  harmony_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Harmony','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  bbknn_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','bbknn','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  #scbatch_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scbatch','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scanorama_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Scanorama','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')
  scvi_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scvi','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')  
  raw_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','Raw','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',')  
  combat_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','combat','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  limma_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','limma','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  scAEQN_df <- read.table(paste0(out_dir,eval_metric,'lisi_batch','_',plx,'_','scAEQN','.txt'), head=T, row.names = 1, check.names = FALSE, sep = ',') 
  
  dir.create(paste0(this_dir, eval_metric,"result/"), showWarnings = F)
  dir.create(paste0(this_dir, eval_metric,"result/",plx,"/"), showWarnings = F)
  methods_use <- c('raw','mnn','limma','harmony','combat','quantnorm',
                   'scvi','scAEQN','scanorama','bbknn')
  pcLISI <- LISI_boxplot_fun(data =  list(raw_df, mnn_df, limma_df,combat_df,
                                          harmony_df, quant_df, 
                                          bbknn_df,scvi_df,  scanorama_df,  
                                          scAEQN_df),
                             vect_names_method = methods_use, 
                             title='cLISI', toptitle=paste0('cLISI Cell Type Mixing plx ',plx),
                             base_name=paste0(this_dir, eval_metric, "result/",plx,"/"), save_results=T)
  
  pcLISI
  
  save(pcLISI, file = paste0(this_dir, eval_metric, "result/",plx, "/pcLISI_",dataset_use,"_plots_",plx,".rda"))
}


run_LISI_final_celltype <- function(fn,data_dir, save_dir, eval_metric,methods_use, plx=30){
  #myPCA <- read.table(paste0(this_dir,eval_metric,methods_use,'/',fn), head=T, row.names = 1, check.names = FALSE)
  myPCA <- read.csv(paste0(data_dir,'/',fn), head=T, row.names = 1, check.names = FALSE)
  head(myPCA)
  
  cpcs <- grep('([Pp][Cc]_?)|(V)|(Harmony)|(W)',colnames(myPCA))
  cpcs <- cpcs[1:20]
  lisi_embeddings <- myPCA[,cpcs]
  
  colnames(myPCA)[grep('[cC]ell_?[tT]ype',colnames(myPCA))] <- 'cell_type'
  colnames(myPCA)[grep('([bB]atch)|(BATCH)',colnames(myPCA))] <- 'batch'
  
  lisi_meta_data <- myPCA[,c('batch','cell_type')]
  lisi_label = c('batch', 'cell_type')
  
  lisi_res <- lisi::compute_lisi(lisi_embeddings, lisi_meta_data, lisi_label,perplexity = plx)
  lisi_res$cell <- rownames(lisi_embeddings)
  
  # lisi_batch <- subset(lisi_res,select=c('batch','cell'))
  lisi_celltype <- subset(lisi_res,select=c('cell_type','cell'))
  dir.create(paste0(this_dir, eval_metric, methods_use,'/'), showWarnings = F)
  # write.table(lisi_batch, paste0(this_dir,eval_metric,methods_use,'/','lisi_batch_',plx,'.txt'), quote=F, sep='\t', row.names=T, col.names=NA)
  write.table(lisi_celltype,paste0(this_dir,eval_metric,methods_use,'/','lisi_celltype_',plx,'.txt'), quote=F, sep='\t', row.names=T, col.names=NA)
}
