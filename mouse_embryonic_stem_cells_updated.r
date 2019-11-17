#This figure generates Fig 5 A-F

library(scater)
scenet<-readRDS('https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/kolodziejczyk.rds')
rawdat <- exprs(scenet)
batch <- scenet@colData$batch
cell.type <- scenet@colData$cell_type1

i=1:704
batch_1 <- unlist(strsplit(as.character(batch),'_'))[2*i]
batch_1[batch_1 %in% c("4","5")] = "4"
rawdat <- rawdat[,batch_1 %in% c('2','3')]
batch <- batch_1[batch_1 %in% c('2','3')]
cell.type <- cell.type[batch_1 %in% c('2','3')]
dat.for.correct <- list(dat = rawdat,batch = batch, cell.type = cell.type)

#save(dat.for.correct,file='kolodziejczyk.rdata',version=2)
#load('kolodziejczyk.rdata')
#rawdat <- dat.for.correct$dat
#batch <- dat.for.correct$batch
#cell.type <- dat.for.correct$cell.type

#ComBat
library(sva)
combatmod <- ComBat(rawdat[rowSums(rawdat)>0,],as.numeric(as.factor(batch)))
#save(combatmod,file='kol_combat.rdata',version=2)

#limma
library(limma)
limmamod = removeBatchEffect(rawdat,batch)
#save(limmamod,file='kol_limma.rdata',version=2)

#MNN
library(batchelor)
mnnmod = mnnCorrect(rawdat,batch=as.factor(batch),cos.norm.out = F)
mnnmod = mnnmod@assays@.xData$data$corrected
colnames(mnnmod) = colnames(rawdat)
rownames(mnnmod) = rownames(rawdat)
#save(mnnmod,file='kol_mnn.rdata',version=2)

#rescaleBatches
batchelormod <- rescaleBatches(rawdat,batch=as.factor(batch))  
batchelormod <- batchelormod@assays@.xData$data$corrected
colnames(batchelormod) = colnames(rawdat)
rownames(batchelormod) = rownames(rawdat)
#save(batchelormod,file='kol_rescalebatch_noERCC.rdata',version=2)

#scBatch
library(scBatch)
distmod <- QuantNorm(rawdat,as.numeric(as.factor(batch)),logdat=F,cor_method='pearson',max=5)
scbatchmod <- scBatchCpp(c=rawdat,w=diag(ncol(rawdat)),d=distmod,m=1,max=200,tol=1e-5,step=0.00001,derif=scBatch::derif,verbose=T)
rownames(scbatchmod) = rownames(rawdat)
uniquebat <- unique(as.character(batch))
batch <- as.character(batch)
scbatchmod1 = scbatchmod
for (i in 1:nrow(scbatchmod)){
  for (j in 1:2){
    scbatchmod1[i,batch == uniquebat[j]] <- (scbatchmod1[i,batch == uniquebat[j]] - mean(scbatchmod1[i,batch == uniquebat[j]]) + mean(scbatchmod1[i,batch == uniquebat[1]]))
  }
}
rownames(scbatchmod1) = rownames(rawdat)
#save(scbatchmod1,file='kol_scbatch_max5_m1.rdata',version=2)

#scPLS
library(Citrus)
res <- scPLS(t(rawdat[1:38524,]), t(rawdat[38525:38616,]), k1 = 1, k2 = 2, iter = 300)
#save(res,file='scplsres_newk.rdata',version=2)
scplsmod <- res$Adjusted
scplsmod <- t(scplsmod)


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
km(cell.type,cor(rawdat),3)
#0.151

#ComBat
km(cell.type,cor(combatmod),3)
#0.231

#MNN
km(cell.type,cor(mnnmod),3)
#0.010

#QuantNorm
km(cell.type,distmod,3)
# 0.845

#scBatch
km(cell.type,cor(scbatchmod1),3)
#0.547

#rescaleBatches
km(cell.type,cor(batchelormod),3)
#0.252

#limma
km(cell.type,cor(limmamod),3)
#0.200

# scPLS
km(cell.type,cor(scplsmod),3)
#0.240


library(ggplot2)
library(gridExtra)
library(Rtsne)

###########################################
# Figure 5A

tsne <- Rtsne(cor(rawdat),is.distance=T)
plotdat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('Uncorrected data') + theme_classic()+theme(title = element_text(size = rel(1.5)),
                                                                                                                                           axis.title.x = element_text(size = rel(1)),
                                                                                                                                           axis.title.y = element_text(size = rel(1)),
                                                                                                                                           legend.text = element_text(size = rel(1.5)),
                                                                                                                                           legend.title = element_text(size = rel(1)))

tsne <- Rtsne(cor(combatmod),is.distance=T)
plotcombat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat') + theme_classic()+theme(title = element_text(size = rel(1.5)),
                                                                                                                                    axis.title.x = element_text(size = rel(1)),
                                                                                                                                    axis.title.y = element_text(size = rel(1)),
                                                                                                                                    legend.text = element_text(size = rel(1.5)),
                                                                                                                                    legend.title = element_text(size = rel(1)))

tsne <- Rtsne(cor(mnnmod),is.distance=T)
plotmnn <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p3 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('MNN') + theme_classic()+theme(title = element_text(size = rel(1.5)),
                                                                                                                              axis.title.x = element_text(size = rel(1)),
                                                                                                                              axis.title.y = element_text(size = rel(1)),
                                                                                                                              legend.text = element_text(size = rel(1.5)),
                                                                                                                              legend.title = element_text(size = rel(1)))

tsne <- Rtsne(cor(scbatchmod1),is.distance=T)
plotscbatch <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p4 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scBatch') + theme_classic()+theme(title = element_text(size = rel(1.5)),
                                                                                                                                      axis.title.x = element_text(size = rel(1)),
                                                                                                                                      axis.title.y = element_text(size = rel(1)),
                                                                                                                                      legend.text = element_text(size = rel(1.5)),
                                                                                                                                      legend.title = element_text(size = rel(1)))

tsne <- Rtsne(cor(limmamod),is.distance=T)
plotlimma <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p5 = ggplot(plotlimma,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('limma') + theme_classic()+theme(title = element_text(size = rel(1.5)),
                                                                                                                                  axis.title.x = element_text(size = rel(1)),
                                                                                                                                  axis.title.y = element_text(size = rel(1)),
                                                                                                                                  legend.text = element_text(size = rel(1.5)),
                                                                                                                                  legend.title = element_text(size = rel(1)))

tsne <- Rtsne(cor(batchelormod),is.distance=T)
plotbatchelore <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p6 = ggplot(plotbatchelore,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('rescaleBatches') + theme_classic()+theme(title = element_text(size = rel(1.5)),
                                                                                                                                                axis.title.x = element_text(size = rel(1)),
                                                                                                                                                axis.title.y = element_text(size = rel(1)),
                                                                                                                                                legend.text = element_text(size = rel(1.5)),
                                                                                                                                                legend.title = element_text(size = rel(1)))


tsne <- Rtsne(cor(scplsmod),is.distance=T)
plotscpls <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p7 = ggplot(plotscpls,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scPLS') + theme_classic()+theme(title = element_text(size = rel(1.5)),
                                                                                                                                  axis.title.x = element_text(size = rel(1)),
                                                                                                                                  axis.title.y = element_text(size = rel(1)),
                                                                                                                                  legend.text = element_text(size = rel(1.5)),
                                                                                                                                  legend.title = element_text(size = rel(1)))

ggsave("Review_Figures/kol_cell_tsne.png",plot=grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="png",height=8,width=16)
ggsave("Review_Figures/kol_cell_tsne.eps",plot=grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="eps",height=8,width=16,dpi=300)


################## PCA ####################

pca <- princomp(cor(rawdat))$scores[,1:2]
plotdat <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p1 = ggplot(plotdat,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('Uncorrected data') + theme_classic()

pca <- princomp(cor(combatmod))$scores[,1:2]
plotcombat <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p2 = ggplot(plotcombat,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat') + theme_classic()

pca <- princomp(cor(mnnmod))$scores[,1:2]
plotmnn <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p3 = ggplot(plotmnn,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('MNN') + theme_classic()

pca <- princomp(cor(scbatchmod1))$scores[,1:2]
plotscbatch <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p4 = ggplot(plotscbatch,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scBatch') + theme_classic()+
  scale_y_continuous(limits=c(-0.7, 0.6)) 

pca <- princomp(cor(limmamod))$scores[,1:2]
plotlimma <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p5 = ggplot(plotlimma,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('limma') + theme_classic()

pca <- princomp(cor(batchelormod))$scores[,1:2]
plotbatchelor <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p6 = ggplot(plotbatchelor,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('rescaleBatches - no ERCC') + theme_classic()

pca <- princomp(cor(scplsmod))$scores[,1:2]
plotscpls <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p7 = ggplot(plotscpls,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scPLS') + theme_classic()

ggsave("Review_Figures/kol_cell_PCA.png",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="png",height=8,width=16)
ggsave("Review_Figures/kol_cell_PCA.eps",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="eps",height=8,width=16,dpi=300)


#####################################################


#UMAP

library(umap)

UMAP <- umap(cor(rawdat))
plotdat <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p1 = ggplot(plotdat,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('Uncorrected data') + theme_classic()

UMAP <- umap(cor(combatmod))
plotcombat <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p2 = ggplot(plotcombat,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat') + theme_classic()

UMAP <- umap(cor(mnnmod))
plotmnn <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p3 = ggplot(plotmnn,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('MNN') + theme_classic()

UMAP <- umap(cor(scbatchmod1))
plotscbatch <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p4 = ggplot(plotscbatch,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scBatch') + theme_classic()

UMAP <- umap(cor(limmamod))
plotlimma <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p5 = ggplot(plotlimma,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('limma') + theme_classic()

UMAP <- umap(cor(batchelormod))
plotbatchelor <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p6 = ggplot(plotbatchelor,aes(x=UMAP1 , y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('rescaleBatches - no ERCC') + theme_classic()

UMAP <- umap(cor(scplsmod))
plotscpls <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p7 = ggplot(plotscpls,aes(x=UMAP1 , y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scPLS') + theme_classic()

ggsave("Review_Figures/kol_cell_umap.png",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="png",height=8,width=16)
ggsave("Review_Figures/kol_cell_umap.eps",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="eps",height=8,width=16,dpi=300)


################# kBET ##################

library(kBET)

kBETraw <- kBET(rawdat,batch,k0=20)
kBETcombat <- kBET(combatmod,batch,k0=20)
kBETmnn <- kBET(mnnmod,batch,k0=20)
kBETscbatch <- kBET(scbatchmod1,batch,k0=20)
kBETlimma <- kBET(limmamod,batch,k0=20)
kBETbatchelor<- kBET(batchelormod,batch,k0=20)
kBETscpls <- kBET(scplsmod,batch,k0=20)

kBETdisplay <- data.frame(kBET = c(kBETraw$stats$kBET.observed,
                                   kBETcombat$stats$kBET.observed,
                                   kBETmnn$stats$kBET.observed,
                                   kBETscbatch$stats$kBET.observed,
                                   kBETlimma$stats$kBET.observed,
                                   kBETbatchelor$stats$kBET.observed,
                                   kBETscpls$stats$kBET.observed),
                          Methods = c(rep('Raw data',100),
                                      rep('ComBat',100),
                                      rep('MNN',100),
                                      rep('scBatch',100),
                                      rep('limma',100),
                                      rep('rescaleBatches',100),
                                      rep('scPLS',100)))

kBETdisplay$Methods <- factor(kBETdisplay$Methods,levels=c('Raw data','ComBat',
                                                           'MNN','scBatch','limma',
                                                           'rescaleBatches',
                                                           'scPLS'))

pkbet = ggplot(data = kBETdisplay, 
               aes(x=Methods, y=kBET)) + geom_boxplot(aes(fill=Methods)) + 
  xlab('Methods') + ylab('Rejection Rates') + labs(title=paste("kBET rejection rates when k=20")) + ylim(0,1)
pkbet  + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Review_Figures/kol_kBET.png",plot=pkbet + theme(axis.text.x = element_text(angle = 90, hjust = 1)),device="png",height=6,width=6)
ggsave("Review_Figures/kol_kBET.eps",plot=pkbet + theme(axis.text.x = element_text(angle = 90, hjust = 1)),device="eps",height=6,width=6,dpi=300)


##########################################################
##########################################################
#DE test with Seurat

DEmnn <- list()
DEraw <- list()
DEscbatch <- list()
DEcombat <- list()
DElimma <- list()
DErescale <- list()
DEscpls <- list()

colnames(rawdat) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(mnnmod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(combatmod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(scbatchmod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(limmamod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(batchelormod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(scplsmod) <- paste0(cell.type, "__", 1:ncol(rawdat))



library(Seurat)
for (i in 1:3){
  groups = combn(c('2i','a2i','lif'),2)[,i]
  
  d = CreateSeuratObject(counts = (mnnmod-min(mnnmod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEmnn[[i]] <- rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(counts = (rawdat-min(rawdat))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEraw[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(counts = (scbatchmod-min(scbatchmod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEscbatch[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(counts = (combatmod-min(combatmod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEcombat[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(counts = (limmamod-min(limmamod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DElimma[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(counts = (batchelormod-min(batchelormod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DErescale[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
  
  d = CreateSeuratObject(counts = (scplsmod-min(scplsmod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEscpls[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])
}

save(DEmnn,file='Review_Figures/DEmnnfiltered_kol.rdata')
save(DEscbatch,file='Review_Figures/DEscbatchfiltered_kol.rdata')
save(DEcombat,file='Review_Figures/DEcombatfiltered_kol.rdata')
save(DEraw,file='Review_Figures/DErawfiltered_kol.rdata')
save(DElimma,file='Review_Figures/DElimmafiltered_kol.rdata')
save(DErescale,file='Review_Figures/DErescalenoerccfiltered_kol.rdata')
save(DEscpls,file='Review_Figures/DEscplsfiltered_kol.rdata')


library(eulerr)
V = list()
library(gridExtra)
for (i in 1:3){
  groups = combn(c('2i','a2i','lif'),2)[,i]
  Vstem <- euler(list(Uncorrected = as.vector(DEraw[[i]]), ComBat = as.vector(DEcombat[[i]]), 
                      MNN = as.vector(DEmnn[[i]]),scBatch = as.vector(DEscbatch[[i]]),
                      limma = as.vector(DElimma[[i]]), 
                      rescaleBatches = as.vector(DErescale[[i]]), scPLS = as.vector(DEscpls[[i]])))
  
  if (i %in% c(2)){
    V[[i]] <- plot(Vstem,quantities=T,lty=1:6,label=F,main=paste(groups[1],'vs',groups[2],sep=' '),legend=list(nrow=4,ncol=2,fontsize=15,font=1))
  }else{
    V[[i]] <- plot(Vstem,quantities=T,lty=1:6,label=F,main=paste(groups[1],'vs',groups[2],sep=' '),legend=F)
  }}

grid.arrange(V[[1]],V[[2]],V[[3]],nrow=1)
ggsave("Review_Figures/kol_venn.png",plot=grid.arrange(V[[1]],V[[2]],V[[3]],nrow=1)
       ,device="png",height=4,width=16,dpi=300)
ggsave("Review_Figures/kol_venn.eps",plot=grid.arrange(V[[1]],V[[2]],V[[3]],nrow=1)
       ,device="eps",height=4,width=16,dpi=300)


########################################################################
# Gene Ontology analysis using GOstats

# First save DE gene list for the cell type pairs
# Each file contains detected DE genes from Uncorrected, MNN, ComBat and scBatch
# Please consider save data in a new directory in convenience of next steps
for (i in 1:3){
  groups = combn(c('2i','a2i','lif'),2)[,i]
  DElist<-list(Uncorrected = as.vector(DEraw[[i]]), ComBat = as.vector(DEcombat[[i]]), 
               MNN = as.vector(DEmnn[[i]]),scBatch = as.vector(DEscbatch[[i]]),
               limma = as.vector(DElimma[[i]]), rescaleBatches = as.vector(DErescale[[i]]),
               scPLS = as.vector(DEscpls[[i]]))
  save(DElist,file=paste('Revision_Figures/kol/DElist/kol_',groups[1],'_',groups[2],'.rdata',sep=''))
}

library(org.Mm.eg.db)
library(GOstats)

# Match gene names from raw data with database
all.genes<-rownames(rawdat)

library(mygene)
all.entrez <- getGenes(all.genes)
all.entrez <- all.entrez@listData$entrezgene
all.entrez.2<-all.genes
for(i in 1:length(all.entrez.2)) all.entrez.2[i]<-all.entrez[[i]][1]
all.entrez<-all.entrez.2

# Main function to conduct GO analysis
get.table<-function(sel.genes, all.genes, all.entrez)
{
  sel.entrez<-all.entrez[which(all.genes %in% sel.genes)]
  params <- new("GOHyperGParams", geneIds=sel.entrez[!is.na(sel.entrez)], universeGeneIds=all.entrez[!is.na(all.entrez)], ontology="BP", pvalueCutoff=0.005,conditional=F, testDirection="over", annotation="org.Mm.eg.db")
  over = hyperGTest(params)
  ov<-summary(over)
  ov<-ov[ov[,6]<=500 & ov[,6]>=10,]
  for(i in 2:4) ov[,i]<-signif(ov[,i], 3)
  ov
}


setwd('Revision_Figures/kol/DElist')
# Obtain the file names of saved DE gene lists (should have 3 .rdata files)
files<-dir(pattern=".rdata")

# Conduct and save GO analysis results for all 3 cell type pairs
for(i in 1:3)
{
  # Read saved data
  load(files[i])
  GOlist<-new("list")
  for(k in 1:8)
  {
    GOlist[[k]]<-get.table(DElist[[k]], all.genes, all.entrez)
  }
  names(GOlist)<-names(DElist)
  
  b<-GOlist[[1]]
  b<-rbind(c(names(GOlist)[1],rep("",6)), b)
  for(k in 2:8)
  {
    b<-rbind(b, c(names(GOlist)[k],rep("",6)))
    b<-rbind(b, GOlist[[k]])
  }
  
  write.table(b, paste(strsplit(files[i], ".rdata")[[1]][1],".GO 005.txt"), sep="\t", quote=F, col.names=F, row.names=F)
}
