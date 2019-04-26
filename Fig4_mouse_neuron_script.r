#This figure generates Fig 4 A-F

library(scater)

scenet<-readRDS(url('https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/usoskin.rds'))

rawdat <- exprs(scenet)
batch <- scenet@colData$Library
cell.type <- scenet@colData$cell_type1

rawdat <- rawdat[,!(batch %in% c('L281','L282'))]
cell.type <- cell.type[!(batch %in% c('L281','L282'))]
batch <- batch[!(batch %in% c('L281','L282'))]

#ComBat
library(sva)
combatmod <- ComBat(rawdat[rowSums(rawdat)>0,],as.numeric(batch))

#MNN
library(scran)
batch128 <- rawdat[,batch=='L128']
batch129 <- rawdat[,batch=='L129']
batch130 <- rawdat[,batch=='L130']
batch141 <- rawdat[,batch=='L141']
batch142 <- rawdat[,batch=='L142']
batch226 <- rawdat[,batch=='L226']
batch227 <- rawdat[,batch=='L227']
batch228 <- rawdat[,batch=='L228']
mnnmod = mnnCorrect(batch128,batch129,batch130,batch141,batch142,batch226,batch227,batch228,cos.norm.out = F)
mnnmod = cbind(mnnmod$corrected[[1]],mnnmod$corrected[[2]],mnnmod$corrected[[3]],
                 mnnmod$corrected[[4]],mnnmod$corrected[[5]],mnnmod$corrected[[6]],
                 mnnmod$corrected[[7]],mnnmod$corrected[[8]])
colnames(mnnmod) = colnames(rawdat)

#scBatch
library(scBatch)
distmod <- QuantNorm(rawdat,as.numeric(batch),logdat=F,cor_method='pearson',max=5)
scbatchmod <- scBatchCpp(c=rawdat,w=diag(ncol(rawdat)),d=distmod,m=10,max=200,tol=1e-10,step=0.00001,derif=scBatch::derif,verbose=T)
rownames(scbatchmod) = rownames(rawdat)

#Clustering performance

km <- function(cell.type,ccc,k){
  #reports average ARI for 50 times of k-means clustering
  require(mclust)
  ARI = 0
  for (i in 1:50){
    kms <- kmeans(ccc,k)
    ARI = ARI + adjustedRandIndex(kms$cluster,cell.type)
  }
  ARI = ARI/50
  return(ARI)
}

#uncorrected data
km(cell.type,cor(rawdat),4)
#0.094

#ComBat
km(cell.type,cor(combatmod),4)
#0.111

#MNN
km(cell.type,cor(mnnmod),4)
#0.083

#QuantNorm
#km(cell.type,distmod,4)
#0.664

#scBatch
km(cell.type,cor(scbatchmod),4)
#0.719

library(ggplot2)
library(gridExtra)
library(Rtsne)


###########################################

# Figure 4A

tsne <- Rtsne(cor(rawdat))
plotdat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('Uncorrected data') + theme_classic()

tsne <- Rtsne(cor(combatmod))
plotcombat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat') + theme_classic()

tsne <- Rtsne(cor(mnnmod))
plotmnn <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('MNN') + theme_classic()

tsne <- Rtsne(cor(scbatchmod))
plotscbatch <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scBatch') + theme_classic()

grid.arrange(p1,p2,p3,p4,nrow=2)

###########################################

# Figure 4E

plotdat$Tac1 = rawdat['Tac1',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$Tac1 = combatmod['Tac1',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$Tac1 = mnnmod['Tac1',]
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$Tac1 = scbatchmod['Tac1',]
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

grid.arrange(p1,p2,p3,p4,nrow=2)

###########################################
# Figure 4F

plotdat$Th = rawdat['Th',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$Th = combatmod['Th',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$Th = mnnmod['Th',]
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$Th = scbatchmod['Th',]
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

grid.arrange(p1,p2,p3,p4,nrow=2)

###########################################
# Figure 4C

plotdat$Nefh = rawdat['Nefh',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$Nefh = combatmod['Nefh',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$Nefh = mnnmod['Nefh',]
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$Nefh = scbatchmod['Nefh',]
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

grid.arrange(p1,p2,p3,p4,nrow=2)

###########################################
# Figure 4D

plotdat$Mrgprd = rawdat['Mrgprd',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$Mrgprd = combatmod['Mrgprd',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$Mrgprd = mnnmod['Mrgprd',]
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$Mrgprd = scbatchmod['Mrgprd',]
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

grid.arrange(p1,p2,p3,p4,nrow=2)



#DE test with Seurat

DEmnn <- list()
DEraw <- list()
DEscbatch <- list()
DEcombat <- list()

colnames(rawdat) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(mnnmod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(combatmod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(scbatchmod) <- paste0(cell.type, "__", 1:ncol(rawdat))


library(Seurat)
for (i in 1:6){
  groups = combn(c('NP','NF','TH','PEP'),2)[,i]

  d = CreateSeuratObject(raw.data = (mnnmod-min(mnnmod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEmnn[[i]] <- rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(raw.data = (rawdat-min(rawdat))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEraw[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(raw.data = (scbatchmod-min(scbatchmod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEscbatch[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(raw.data = (combatmod-min(combatmod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEcombat[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
}

#save(DEmnn,file='DEmnn610_usoskin.rdata')
#save(DEscbatch,file='DEscbatch610_usoskin.rdata')
#save(DEcombat,file='DEcombat610_usoskin.rdata')
#save(DEraw,file='DEraw610_usoskin.rdata')


#######################################################################
# Figure 4B

library(eulerr)
V = list()
library(gridExtra)
for (i in 1:6){
  group = combn(c('NP','NF','TH','PEP'),2)[,i]
  Vstem <- euler(list(M = as.vector(DEmnn[[i]]),S = as.vector(DEscbatch[[i]]),
                      C = as.vector(DEcombat[[i]]),R = as.vector(DEraw[[i]])))
  
  V[[i]] <- plot(Vstem,quantities=T,lty=1:4,label=F,main=paste(group[1],'vs',group[2],sep=' '))
}
grid.arrange(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],nrow=2)

for (i in 1:6){
  group = combn(c('NP','NF','TH','PEP'),2)[,i]
  DElist<-list(M = as.vector(DEmnn[[i]]),S = as.vector(DEscbatch[[i]]),
               C = as.vector(DEcombat[[i]]),R = as.vector(DEraw[[i]]))
  save(DElist,file=paste('Usoskin610_',group[1],'_',group[2],'.rdata',sep=''))
}
