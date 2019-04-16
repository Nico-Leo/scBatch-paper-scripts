#This figure generates Fig 5 A-F

library(scater)

dat <- read.table('GSE81608_human_islets_rpkm.txt', header = T)
genes <- read.csv("human_gene_annotation.csv", header = T)
rownames(dat) <- genes[,2]
dat <- dat[,2:ncol(dat)]
filter <- rownames(dat[rowSums(dat >= 100) >= 10,])

scenet<-readRDS(url('https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/xin.rds'))

rawdat <- exprs(scenet)
batch <- scenet@colData$donor.id
cell.type <- scenet@colData$cell_type1
condition <- scenet@colData$condition

cell.type <- cell.type[condition=='Healthy']
batch <- batch[condition=='Healthy']
rawdat <- rawdat[filter,condition=='Healthy']

#ComBat
library(sva)
combatmod <- ComBat(rawdat[rowSums(rawdat)>0,],as.numeric(as.factor(batch)))

#MNN
library(scran)
batch1 <- rawdat[filter,batch=='Non T2D 1']
batch2 <- rawdat[filter,batch=='Non T2D 2']
batch3 <- rawdat[filter,batch=='Non T2D 3']
batch4 <- rawdat[filter,batch=='Non T2D 4']
batch5 <- rawdat[filter,batch=='Non T2D 5']
batch6 <- rawdat[filter,batch=='Non T2D 6']
batch7 <- rawdat[filter,batch=='Non T2D 7']
batch8 <- rawdat[filter,batch=='Non T2D 8']
batch9 <- rawdat[filter,batch=='Non T2D 9']
batch10 <- rawdat[filter,batch=='Non T2D 10']
batch11 <- rawdat[filter,batch=='Non T2D 11']
batch12 <- rawdat[filter,batch=='Non T2D 12']
mnnmod = mnnCorrect(batch1,batch2,batch3,batch4,batch5,batch6,batch7,batch8,
                    batch9,batch10,batch11,batch12,cos.norm.out = F)
mnnmod = cbind(mnnmod$corrected[[1]],mnnmod$corrected[[2]],mnnmod$corrected[[3]],
               mnnmod$corrected[[4]],mnnmod$corrected[[5]],mnnmod$corrected[[6]],
               mnnmod$corrected[[7]],mnnmod$corrected[[8]],mnnmod$corrected[[9]],
               mnnmod$corrected[[10]],mnnmod$corrected[[11]],mnnmod$corrected[[12]])
colnames(mnnmod) = colnames(rawdat)

#scBatch
library(scBatch)
distmod <- QuantNorm(rawdat,as.numeric(as.factor(batch)),logdat=F,cor_method='pearson',max=5)
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
#0.423

#ComBat
km(cell.type,cor(combatmod),4)
#0.442

#MNN
km(cell.type,cor(mnnmod),4)
#-0.005

#QuantNorm
#km(cell.type,distmod,4)
#0.555

#scBatch
km(cell.type,cor(scbatchmod),4)
#0.604


###########################################
# Figure 5A

tsne <- Rtsne(cor(rawdat))
plotdat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2],color=batch)
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=batch)) + ggtitle('Uncorrected data') + theme_classic()

tsne <- Rtsne(cor(combatmod))
plotcombat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2],color=batch)
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=batch)) + ggtitle('ComBat') + theme_classic()

tsne <- Rtsne(cor(mnnmod))
plotmnn <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2],color=batch)
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=batch)) + ggtitle('MNN') + theme_classic()

tsne <- Rtsne(cor(scbatchmod))
plotscbatch <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2],color=batch)
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=batch)) + ggtitle('scBatch') + theme_classic()

grid.arrange(p1,p2,p3,p4,nrow=2)

###########################################
# Figure 5C

plotdat$GCG = rawdat['GCG',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=GCG)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$GCG = combatmod['GCG',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=GCG)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$GCG = mnnmod['GCG',]
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=GCG)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$GCG = scbatchmod['GCG',]
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=GCG)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

grid.arrange(p1,p2,p3,p4,nrow=2)

###########################################
# Figure 5D

plotdat$INS = rawdat['INS',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=INS)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$INS = combatmod['INS',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=INS)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$INS = mnnmod['INS',]
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=INS)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$INS = scbatchmod['INS',]
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=INS)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

grid.arrange(p1,p2,p3,p4,nrow=2)

###########################################
# Figure 5E

plotdat$PPY = rawdat['PPY',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=PPY)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$PPY = combatmod['PPY',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=PPY)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$PPY = mnnmod['PPY',]
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=PPY)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$PPY = scbatchmod['PPY',]
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=PPY)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

grid.arrange(p1,p2,p3,p4,nrow=2)

###########################################
# Figure 5F

plotdat$SST = rawdat['SST',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=SST)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$SST = combatmod['SST',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=SST)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$SST = mnnmod['SST',]
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=SST)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$SST = scbatchmod['SST',]
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=SST)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

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
  groups = combn(c('alpha','beta','delta','gamma'),2)[,i]
  
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

#save(DEmnn,file='DEmnnfiltered_xin.rdata')
#save(DEscbatch,file='DEscbatchfiltered_xin.rdata')
#save(DEcombat,file='DEcombatfiltered_xin.rdata')
#save(DEraw,file='DErawfiltered_xin.rdata')

############################################################################
# Figure 5B

library(eulerr)
V = list()
library(gridExtra)
for (i in 1:6){
  groups = combn(c('alpha','beta','delta','gamma'),2)[,i]
  Vstem <- euler(list(M = as.vector(DEmnn[[i]]),S = as.vector(DEscbatch[[i]]),
                      C = as.vector(DEcombat[[i]]),R = as.vector(DEraw[[i]])))
  
  V[[i]] <- plot(Vstem,quantities=T,lty=1:4,label=F,main=paste(groups[1],'vs',groups[2],sep=' '))
}
grid.arrange(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],nrow=2)

for (i in 1:6){
  groups = combn(c('alpha','beta','delta','gamma'),2)[,i]
  DElist<-list(M = as.vector(DEmnn[[i]]),S = as.vector(DEscbatch[[i]]),
               C = as.vector(DEcombat[[i]]),R = as.vector(DEraw[[i]]))
  save(DElist,file=paste('xinfiltered_',groups[1],'_',groups[2],'.rdata',sep=''))
}
