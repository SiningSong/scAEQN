rm(list = ls())
cat("\014")

#######################导入数据######################
#导入数据
sce=readRDS('E:/mymodel/data/human/sce_human.rds')
uncorr_data=sce@assays@data@listData[["counts"]]
batch <- sce@colData$batch
celltype <- sce@colData$celltype
label=array(1:651)
for(i in 1:651) {
  if(celltype[i]%in%'beta'|celltype[i]%in%'beta.contaminated'=='TRUE'){
    label[i]='beta'
  }else if(celltype[i]%in%'alpha'|celltype[i]%in%'alpha.contaminated'=='TRUE'){
    label[i]='alpha'
  }else if(celltype[i]%in%'delta'|celltype[i]%in%'delta.contaminated'=='TRUE'){
    label[i]='delta'
  }else if(celltype[i]%in%'gamma'|celltype[i]%in%'gamma.contaminated'=='TRUE'){
    label[i]='gamma'
  }
}

scaeqn_data=read.csv('E:/mymodel/methodcompare/datasets/human pancreas/humanpanc_scaeqn_output.csv',row.names = 1)
########################seurat方法#####################
library(Seurat)
library(dplyr)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(readr)
####uncorrected
data <- CreateSeuratObject(counts = uncorr_data)   #构建seurat数据集
Idents(data) <- label
pbmc.markers <- FindAllMarkers(object = data,logfc.threshold=0.05)
#取每种细胞类型按log2FC值从大到小排序的前五个差异基因，用来绘制热图
top5 <- pbmc.markers%>%group_by(cluster)%>%top_n(n = 5, wt = avg_log2FC)
mat<-GetAssayData(data, slot = "counts")
gene_features <- top5
cluster_info <- sort(data@active.ident)
mat <- as.matrix(mat[top5$gene, names(cluster_info)])
#计算得到的差异基因在调整后的P值<1e-5且log2FC>2的基因个数
FC <- pbmc.markers$avg_log2FC[pbmc.markers$p_val_adj < 1e-2]
names(FC) <- rownames(pbmc.markers[pbmc.markers$p_val_adj < 1e-2,])
DEs_uncorr <- list()
#DEuncorr[] <- rownames(pbmc.markers[pbmc.markers$p_val_adj < 1e-5 & pbmc.markers$avg_log2FC > 2,])
DEs_uncorr <- pbmc.markers$gene[pbmc.markers$p_val_adj < 1e-2 & pbmc.markers$avg_log2FC > 0]
#1074

####scAEQN
data_scaeqn <- CreateSeuratObject(counts = scaeqn_data)
Idents(data_scaeqn) <- label
pbmc.markers_scaeqn <- FindAllMarkers(object = data_scaeqn,logfc.threshold=0.05)
##取每种细胞类型按log2FC值从大到小排序的前五个差异基因，用来绘制热图
top5_scaeqn <- pbmc.markers_scaeqn%>%group_by(cluster)%>%top_n(n = 5, wt = avg_log2FC)
mat_scaeqn<-GetAssayData(data_scaeqn, slot = "counts")
gene_features_scaeqn <- top5_scaeqn
cluster_info_scaeqn <- sort(data_scaeqn@active.ident)
mat_scaeqn <- as.matrix(mat_scaeqn[top5_scaeqn$gene, names(cluster_info_scaeqn)])
#计算得到的差异基因在调整后的P值<1e-5且log2FC>2的基因个数
FC_scaeqn <- pbmc.markers_scaeqn$avg_log2FC[pbmc.markers_scaeqn$p_val_adj < 1e-2]
names(FC_scaeqn) <- rownames(pbmc.markers_scaeqn[pbmc.markers_scaeqn$p_val_adj < 1e-2,])
DEs_scaeqn <- list()
#DEuncorr[] <- rownames(pbmc.markers[pbmc.markers$p_val_adj < 1e-5 & pbmc.markers$avg_log2FC > 2,])
DEs_scaeqn <- pbmc.markers_scaeqn$gene[pbmc.markers_scaeqn$p_val_adj < 1e-2 & pbmc.markers_scaeqn$avg_log2FC > 0]
#4892
###########################热图###########################
png('E:/mymodel/methodcompare/results/human pancreas/DEG_uncor.png',width = 2200, height = 1500,res=300)
p1 <- Heatmap(mat,                          #数据矩阵
              #width = unit(5, "cm"),        #图片宽度
              #height = unit(5, "cm"),     #图片高度
              name = "log2FC",
              border = TRUE,
              color_space = 'RGB',
              #rect_gp = gpar(col = "blue", lwd = 0.01),  
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              row_names_gp = gpar(family = 'Broman',fontsize = 5.5),
              column_split = cluster_info,
              top_annotation = HeatmapAnnotation(
                Cluster = cluster_info,
                annotation_legend_param =(cluster=list(anno_block(gp = gpar(fill = col),
                                                                  labels = levels(cluster_info))
                ))),
              column_title = 'Cell Type',
              column_title_gp = gpar(family = 'Broman', fontsize = 10),
              row_title = "Gene Makers",
              row_title_gp = gpar(family = 'Broman', fontsize = 10),)
print(p1)
dev.off()
png('E:/mymodel/methodcompare/results/human pancreas/DEG_scAEQN.png',width = 2200, height = 1500,res=300)
p2 <- Heatmap(mat_scaeqn,                          #数据矩阵
              name = "log2FC",
              border = TRUE,
              #rect_gp = gpar(col = "blue", lwd = 0.01),  
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              row_names_gp = gpar(family = 'Broman',fontsize = 5.5, fontface = "bold"),
              column_split = cluster_info_scaeqn,
              top_annotation = HeatmapAnnotation(
                Cluster = cluster_info_scaeqn,
                annotation_legend_param =(cluster=list(anno_block(gp = gpar(fill=col),
                                                                  labels = levels(cluster_info_scaeqn))
                ))),
              column_title = 'Cell Type',
              column_title_gp = gpar(family = 'Broman', fontsize = 10),
              row_title = "Gene Makers",
              row_title_gp = gpar(family = 'Broman', fontsize = 10),)
print(p2)
dev.off()
