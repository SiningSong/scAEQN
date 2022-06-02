rm(list = ls())
cat("\014")
#在人鼠脑细胞上的比较
#对比方法：MNN,limma,combat,scbatch,quantnorm

library(scater)
brain <- readRDS('E:/mymodel/data/brain/sce_brain.rds')
rawdat <- brain@assays@data@listData[["counts"]]
batch <- brain@colData$batch
celltype <- brain@colData$celltype

#MNN
t1=proc.time()
library(batchelor)
mnnmod = mnnCorrect(rawdat,batch=as.factor(batch),cos.norm.out = F)
mnnmod = mnnmod@assays@data@listData[['corrected']]
colnames(mnnmod) = colnames(rawdat)
rownames(mnnmod) = rownames(rawdat)
t2=proc.time()
tmnn=t2-t1
save(mnnmod,file = 'E:/mymodel/methodcompare/brain/brain_mnn.rdata')

#limma
t3=proc.time()
library(limma)
limmamod = removeBatchEffect(rawdat,batch)
colnames(limmamod) = colnames(rawdat)
t4=proc.time()
tlimma=t4-t3
save(limmamod,file = 'E:/mymodel/methodcompare/brain/brain_limma.rdata')

#ComBat
t5=proc.time()
library(sva)
combatmod <- ComBat(rawdat[rowSums(rawdat)>0,],as.numeric(as.factor(batch)))
colnames(combatmod) = colnames(rawdat)
t6=proc.time()
tcombat=t6-t5
save(combatmod,file = 'E:/mymodel/methodcompare/brain/brain_combat.rdata')

#QuantNorm
t7=proc.time()
quantmod <- QuantNorm(rawdat,as.numeric(as.factor(batch)),logdat=F,method='row/column',cor_method='pearson',max=50)
colnames(quantmod) = colnames(rawdat)
t8=proc.time()
tqn=t8-t7
save(quantmod,file = 'E:/mymodel/methodcompare/brain/brain_quant.rdata')


#scBatch
t9=proc.time()
library(scBatch)
distmod <- QuantNorm(rawdat,as.numeric(as.factor(batch)),logdat=F,method='row/column',cor_method='pearson',max=50)
scbatchmod <- scBatchCpp(c=rawdat,w=diag(ncol(rawdat)),d=distmod,m=1,max=50,tol=1e-5,step=0.00001,derif=scBatch::derif,verbose=T)
rownames(scbatchmod) = rownames(rawdat)
t10=proc.time()
tscbatch=t10-t9
#save(scbatchmod,file='E:/mymodel/brain_scbatch.rdata')
print(c('mnn:',tmnn))
print(c('limma:',tlimma))
print(c('combat:',tcombat))
print(c('quantnorm:',tqn))
print(c('scbatch:',tscbatch))

#Seurat V4
library(SeuratObject)
library(Seurat)
t11 = proc.time()
scdata <- as.Seurat(brain)
# split the object by dataset
scdata.list <- SplitObject(scdata, split.by = "batch")
# normalize and identify variable features for each dataset independently
scdata.list <- lapply(X = scdata.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = scdata.list)

### FindAnchors 
scdata.anchors <- FindIntegrationAnchors(object.list = scdata.list,
                                         normalization.method = "LogNormalize",  #如果前面是log标准化，这里改成LogNormalize
                                         dims = 1:20,
                                         k.anchor = 5,
                                         k.filter = 5,
                                         n.trees = 10,
                                         anchor.features = features)
### Integrate 运行速度慢
scdata.int <- IntegrateData(scdata.anchors, normalization.method="LogNormalize") #速度慢
seuratmod <- scdata.int@assays[['integrated']]@data
write.csv(seuratmod,file = 'E:/mymodel/methodcompare/brain/brain_seurat.csv')
t12 = proc.time()
tseurat = t12-t11
print(tseurat)

km <- function(cell.type,ccc,k){
  #计算10次的ARI平均值
  require(mclust)
  ARI = 0
  for (i in 1:10){
    kms <- kmeans(ccc,k)
    ARI = ARI + adjustedRandIndex(kms$cluster,cell.type)
  }
  ARI = ARI/10
  return(ARI)
}


#uncorrected data
RAWARI <- km(celltype,t(rawdat),6)


#ComBat
CombatARI <- km(celltype,t(combatmod),6)


#MNN
MNNARI <- km(celltype,t(mnnmod),6)


#limma
limmaARI <- km(celltype,t(limmamod),6)

#quantnorm
quantARI <- km(celltype,quanttmod,6)


#scBatch
scbatchARI <- km(celltype,t(scbatchmod),6)  #将数据转制后再聚类


