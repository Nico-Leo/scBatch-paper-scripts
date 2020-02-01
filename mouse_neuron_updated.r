d <- read.csv("Usoskin et al. External resources Table 1.csv", stringsAsFactors = F)
exprs_data <- data.matrix(d[11:nrow(d), 10:ncol(d)])
ann <- t(d[c(1:3, 5:8), 10:ncol(d)])
colnames(ann) <- d[c(1:3, 5:8), 9]
colnames(ann)[5:7] <- c("cell_type1", "cell_type2", "cell_type3")

# filter bad cells
filt <- !grepl("Empty well", ann[,5]) &
  !grepl("NF outlier", ann[,5]) &
  !grepl("TH outlier", ann[,5]) &
  !grepl("NoN outlier", ann[,5]) &
  !grepl("NoN", ann[,5]) &
  !grepl("Central, unsolved", ann[,5]) &
  !grepl(">1 cell", ann[,5]) &
  !grepl("Medium", ann[,5])
colnames(exprs_data) <- rownames(ann)
rownames(exprs_data) <- d[11:nrow(d), 1]

### ANNOTATIONS
ann <- as.data.frame(ann[filt,])
exprs_data <- exprs_data[,filt]

library(scater)

scenet<-readRDS(url('https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/usoskin.rds'))

rawdat <- exprs(scenet)
batch <- scenet@colData$Library
cell.type <- scenet@colData$cell_type1

rawdat <- rawdat[,!(batch %in% c('L281','L282'))]
exprs_data <- exprs_data[,!(batch %in% c('L281','L282'))]
cell.type <- cell.type[!(batch %in% c('L281','L282'))]
batch <- batch[!(batch %in% c('L281','L282'))]

#scBatch
library(scBatch)
distmod <- QuantNorm(rawdat,as.numeric(batch),logdat=F,cor_method='pearson',max=5)
scbatchmod <- scBatchCpp(c=rawdat,w=diag(ncol(rawdat)),d=distmod,m=10,max=1200,tol=1e-10,step=0.00001,derif=scBatch::derif,verbose=T)
uniquebat <- unique(as.character(batch))
batch <- as.character(batch)
scbatchmod1 = scbatchmod
for (i in 1:nrow(scbatchmod)){
  for (j in 2:8){
    scbatchmod1[i,batch == uniquebat[j]] <- (scbatchmod[i,batch == uniquebat[j]] - mean(scbatchmod[i,batch == uniquebat[j]]) + mean(scbatchmod[i,batch == uniquebat[1]]))
  }
}
rownames(scbatchmod1) = rownames(rawdat)

#ComBat & ComBat-seq
library(sva)
combatmod <- ComBat(rawdat[rowSums(rawdat)>0,],as.numeric(as.factor(batch)))
#save(combatmod,file='usoskin_combat.rdata')

combatseqmod <- ComBat_seq(exprs_data,as.numeric(as.factor(batch)),NULL)
sce <- SingleCellExperiment(assays=list(counts=combatseqmod))
sce <- normalize(sce)
combatseqmod <- exprs(sce)
save(combatseqmod,file='usoskin_combatseq.rdata')

#MNN
library(batchelor)
mnnmod = mnnCorrect(rawdat,batch=as.factor(batch),cos.norm.out = F)
mnnmod = mnnmod@assays@.xData$data$corrected
colnames(mnnmod) = colnames(rawdat)
rownames(mnnmod) = rownames(rawdat)
#save(mnnmod,file='usoskin_mnn.rdata')

#rescalebatch
batchelormod <- rescaleBatches(rawdat,batch=as.factor(batch))  
batchelormod <- batchelormod@assays@.xData$data$corrected
colnames(batchelormod) = colnames(rawdat)
rownames(batchelormod) = rownames(rawdat)
#save(batchelormod,file='usoskin_rescalebatch.rdata')

#limma
library(limma)
limmamod = removeBatchEffect(rawdat,batch)
#save(limmamod,file='usoskin_limma.rdata')

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

#ComBat-seq
km(cell.type,cor(combatseqmod),4)
#0.076

#MNN
km(cell.type,cor(mnnmod),4)
#0.149

#QuantNorm
#km(cell.type,distmod,4)
#0.664

#scBatch
km(cell.type,cor(scbatchmod1),4)
#0.636

#rescaleBatches
km(cell.type,cor(batchelormod),4)
#0.085

#limma
km(cell.type,cor(limmamod),4)
#0.125

library(ggplot2)
library(gridExtra)
library(Rtsne)


###########################################

# Figure 4A

tsne <- Rtsne(cor(rawdat))
km(cell.type,tsne$Y,4)
#0.24
plotdat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('Uncorrected data') + theme_classic()+
  theme(title = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1)))

tsne <- Rtsne(cor(combatmod))
km(cell.type,tsne$Y,4)
#0.44
plotcombat <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat') + theme_classic()+
  theme(title = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1)))

tsne <- Rtsne(cor(combatseqmod))
km(cell.type,tsne$Y,4)
#0.472
plotcombatseq <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p3 = ggplot(plotcombatseq,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat-seq') + theme_classic()+
  theme(title = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1)))


tsne <- Rtsne(cor(mnnmod))
km(cell.type,tsne$Y,4)
#0.15
plotmnn <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p4 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('MNN') + theme_classic()+
  theme(title = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1)))


tsne <- Rtsne(cor(scbatchmod1))
km(cell.type,tsne$Y,4)
#0.72
plotscbatch <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p5 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scBatch') + theme_classic()+
  theme(title = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1)))


tsne <- Rtsne(cor(limmamod))
km(cell.type,tsne$Y,4)
#0.30
plotlimma <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p6 = ggplot(plotlimma,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('limma') + theme_classic()+
  theme(title = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1)))


tsne <- Rtsne(cor(batchelormod))
km(cell.type,tsne$Y,4)
#0.28
plotbatchelor <- data.frame(TSNE1=tsne$Y[,1],TSNE2=tsne$Y[,2])
p7 = ggplot(plotbatchelor,aes(x=TSNE1, y=TSNE2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('rescaleBatches') + theme_classic()+
  theme(title = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1)))


ggsave("Revision_Figures/seq_Usoskin/usoskin_cell_tsne.png",plot=grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="png",height=8,width=16,dpi=300)
ggsave("Revision_Figures/seq_Usoskin/usoskin_cell_tsne.eps",plot=grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="eps",height=8,width=16,dpi=300)

###########################################

# Figure 4E

plotdat$Tac1 = rawdat['Tac1',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$Tac1 = combatmod['Tac1',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotcombatseq$Tac1 = combatseqmod['Tac1',]
p3 = ggplot(plotcombatseq,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('ComBat-seq') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$Tac1 = mnnmod['Tac1',]
p4 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$Tac1 = scbatchmod1['Tac1',]
p5 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotlimma$Tac1 = limmamod['Tac1',]
p6 = ggplot(plotlimma,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('limma') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotbatchelor$Tac1 = batchelormod['Tac1',]
p7 = ggplot(plotbatchelor,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('rescaleBatches') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")


ggsave("Revision_Figures/seq_Usoskin/usoskin_Tac1_tsne.png",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="png",height=8,width=16)
#ggsave("Revision_Figures/seq_Usoskin/usoskin_Tac1_tsne.eps",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="eps",height=8,width=16,dpi=300)

###########################################
# Figure 4F

plotdat$Th = rawdat['Th',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$Th = combatmod['Th',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotcombatseq$Th = combatseqmod['Th',]
p3 = ggplot(plotcombatseq,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('ComBat-seq') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$Th = mnnmod['Th',]
p4 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$Th = scbatchmod1['Th',]
p5 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=Th)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotlimma$Th = limmamod['Th',]
p6 = ggplot(plotlimma,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('limma') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotbatchelor$Th = batchelormod['Th',]
p7 = ggplot(plotbatchelor,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('rescaleBatches') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

ggsave("Revision_Figures/seq_Usoskin/usoskin_Th_tsne.png",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="png",height=8,width=16)
#ggsave("Revision_Figures/seq_Usoskin/usoskin_Th_tsne.eps",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="eps",height=8,width=16,dpi=300)

###########################################
# Figure 4C

plotdat$Nefh = rawdat['Nefh',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$Nefh = combatmod['Nefh',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotcombatseq$Nefh = combatseqmod['Nefh',]
p3 = ggplot(plotcombatseq,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('ComBat-seq') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$Nefh = mnnmod['Nefh',]
p4 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$Nefh = scbatchmod1['Nefh',]
p5 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=Nefh)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotlimma$Nefh = limmamod['Nefh',]
p6 = ggplot(plotlimma,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('limma') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotbatchelor$Nefh = batchelormod['Nefh',]
p7 = ggplot(plotbatchelor,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('rescaleBatches') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

ggsave("Revision_Figures/seq_Usoskin/usoskin_Nefh_tsne.png",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="png",height=8,width=16)
#ggsave("Revision_Figures/seq_Usoskin/usoskin_Nefh_tsne.eps",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="eps",height=8,width=16,dpi=300)

###########################################
# Figure 4D

plotdat$Mrgprd = rawdat['Mrgprd',]
p1 = ggplot(plotdat,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('Uncorrected data') + theme_classic() + scale_color_distiller(palette = "YlOrRd")

plotcombat$Mrgprd = combatmod['Mrgprd',]
p2 = ggplot(plotcombat,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('ComBat') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotcombatseq$Mrgprd = combatseqmod['Mrgprd',]
p3 = ggplot(plotcombatseq,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('ComBat-seq') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotmnn$Mrgprd = mnnmod['Mrgprd',]
p4 = ggplot(plotmnn,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('MNN') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotscbatch$Mrgprd = scbatchmod1['Mrgprd',]
p5 = ggplot(plotscbatch,aes(x=TSNE1, y=TSNE2, color=Mrgprd)) + geom_point(size=1) + ggtitle('scBatch') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotlimma$Mrgprd = limmamod['Mrgprd',]
p6 = ggplot(plotlimma,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('limma') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

plotbatchelor$Mrgprd = batchelormod['Mrgprd',]
p7 = ggplot(plotbatchelor,aes(x=TSNE1, y=TSNE2, color=Tac1)) + geom_point(size=1) + ggtitle('rescaleBatches') + theme_classic()+ scale_color_distiller(palette = "YlOrRd")

ggsave("Revision_Figures/seq_Usoskin/usoskin_Mrgprd_tsne.png",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="png",height=8,width=16)
#ggsave("Revision_Figures/seq_Usoskin/usoskin_Mrgprd_tsne.eps",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="eps",height=8,width=16,dpi=300)


################## PCA ####################

pca <- princomp(cor(rawdat))$scores[,1:2]
plotdat <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p1 = ggplot(plotdat,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('Uncorrected data') + theme_classic()

pca <- princomp(cor(combatmod))$scores[,1:2]
plotcombat <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p2 = ggplot(plotcombat,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat') + theme_classic()

pca <- princomp(cor(combatseqmod))$scores[,1:2]
plotcombatseq <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p3 = ggplot(plotcombatseq,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat-seq') + theme_classic()

pca <- princomp(cor(mnnmod))$scores[,1:2]
plotmnn <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p4 = ggplot(plotmnn,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('MNN') + theme_classic()

pca <- princomp(cor(scbatchmod1))$scores[,1:2]
plotscbatch <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p5 = ggplot(plotscbatch,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scBatch') + theme_classic()

pca <- princomp(cor(limmamod))$scores[,1:2]
plotlimma <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p6 = ggplot(plotlimma,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('limma') + theme_classic()

pca <- princomp(cor(batchelormod))$scores[,1:2]
plotbatchelor <- data.frame(PC1=pca[,1],PC2=pca[,2],color=batch)
p7 = ggplot(plotbatchelor,aes(x=PC1 , y=PC2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('rescaleBatches') + theme_classic()


ggsave("Revision_Figures/seq_Usoskin/usoskin_cell_PCA.png",plot=grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="png",height=8,width=16)
#ggsave("Revision_Figures/seq_Usoskin/usoskin_cell_PCA.eps",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="eps",height=8,width=16,dpi=300)


#####################################################


#UMAP

library(umap)

UMAP <- umap(cor(rawdat))
plotdat <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p1 = ggplot(plotdat,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('Uncorrected data') + theme_classic()

UMAP <- umap(cor(combatmod))
plotcombat <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p2 = ggplot(plotcombat,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat') + theme_classic()

UMAP <- umap(cor(combatseqmod))
plotcombatseq <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p3 = ggplot(plotcombatseq,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('ComBat-seq') + theme_classic()

UMAP <- umap(cor(mnnmod))
plotmnn <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p4 = ggplot(plotmnn,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('MNN') + theme_classic()

UMAP <- umap(cor(scbatchmod1))
plotscbatch <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p5 = ggplot(plotscbatch,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('scBatch') + theme_classic()

UMAP <- umap(cor(limmamod))
plotlimma <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p6 = ggplot(plotlimma,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('limma') + theme_classic()

UMAP <- umap(cor(batchelormod))
plotbatchelor <- data.frame(UMAP1=UMAP$layout[,1],UMAP2=UMAP$layout[,2])
p7 = ggplot(plotbatchelor,aes(x=UMAP1, y=UMAP2)) + geom_point(size=1,aes(color=cell.type)) + ggtitle('rescaleBatches') + theme_classic()

ggsave("Revision_Figures/seq_Usoskin/usoskin_cell_umap.png",plot=grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6,p7,nrow=2,ncol=4),device="png",height=8,width=16,dpi=300)
#ggsave("Revision_Figures/seq_Usoskin/usoskin_cell_umap.eps",plot=grid.arrange(p1,p2,p3,p4,p5,p6,p7,nrow=2),device="eps",height=8,width=16,dpi=300)


################# kBET ##################

library(kBET)

kBETraw <- kBET(rawdat,batch,k0=20)
kBETcombat <- kBET(combatmod,batch,k0=20)
kBETcombatseq <- kBET(combatseqmod,batch,k0=20)
kBETmnn <- kBET(mnnmod,batch,k0=20)
kBETscbatch <- kBET(scbatchmod1,batch,k0=20)
kBETlimma <- kBET(limmamod,batch,k0=20)
kBETbatchelor <- kBET(batchelormod,batch,k0=20)
boxplot(kBETraw$stats$kBET.observed,
        kBETcombat$stats$kBET.observed,
        kBETcombatseq$stats$kBET.observed,
        kBETmnn$stats$kBET.observed,
        kBETscbatch$stats$kBET.observed,
        kBETlimma$stats$kBET.observed,
        kBETbatchelor$stats$kBET.observed)
kBETdisplay <- data.frame(kBET = c(kBETraw$stats$kBET.observed,
                                   kBETcombat$stats$kBET.observed,
                                   kBETcombatseq$stats$kBET.observed,
                                   kBETmnn$stats$kBET.observed,
                                   kBETscbatch$stats$kBET.observed,
                                   kBETlimma$stats$kBET.observed,
                                   kBETbatchelor$stats$kBET.observed),
                          Methods = c(rep('Raw data',100),
                                      rep('ComBat',100),
                                      rep('ComBat-seq',100),
                                      rep('MNN',100),
                                      rep('scBatch',100),
                                      rep('limma',100),
                                      rep('rescaleBatches',100)))

kBETdisplay$Methods <- factor(kBETdisplay$Methods,levels=c('Raw data','ComBat','ComBat-seq',
                                                           'MNN','scBatch','limma',
                                                           'rescaleBatches'))

library(RColorBrewer)
pkbet = ggplot(data = kBETdisplay, 
               aes(x=Methods, y=kBET)) + geom_boxplot(aes(fill=Methods)) + 
  xlab('Methods') + ylab('Rejection Rates') + labs(title=paste("kBET rejection rates when k=20")) + ylim(0,1)+ 
  scale_fill_manual(values=brewer.pal(n = 8, name = 'Set3')[1:7])

ggsave("Revision_Figures/seq_Usoskin/usoskin_kBET.png",plot=pkbet,device="png",height=8,width=8,dpi=300)
ggsave("Revision_Figures/seq_Usoskin/usoskin_kBET.eps",plot=pkbet,device="eps",height=4,width=6,dpi=300)

###################################################
###################################################

#DE test with Seurat


DEmnn <- list()
DEraw <- list()
DEscbatch <- list()
DEcombat <- list()
DEcombatseq <- list()
DElimma <- list()
DErescale <- list()

colnames(rawdat) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(mnnmod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(combatmod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(combatseqmod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(scbatchmod1) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(limmamod) <- paste0(cell.type, "__", 1:ncol(rawdat))
colnames(batchelormod) <- paste0(cell.type, "__", 1:ncol(rawdat))


library(Seurat)
for (i in 1:6){
  groups = combn(c('NP','NF','TH','PEP'),2)[,i]
  
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

  d = CreateSeuratObject(counts = (scbatchmod1-min(scbatchmod1))[,cell.type %in% groups],names.field = 1, names.delim = "__")
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

  d = CreateSeuratObject(counts = (combatseqmod-min(combatseqmod))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")

  FC <- res$avg_logFC[res$p_val_adj < 1e-6]
  names(FC) <- rownames(res[res$p_val_adj < 1e-6,])
  DEcombatseq[[i]] = rownames(res[res$p_val_adj < 1e-6 & res$avg_logFC > 2,])

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
}

save(DEmnn,file='Revision_Figures/seq_Usoskin/DEmnn610_usoskin.rdata')
save(DEscbatch,file='Revision_Figures/seq_Usoskin/DEscbatch610_usoskin.rdata')
save(DEcombat,file='Revision_Figures/seq_Usoskin/DEcombat610_usoskin.rdata')
save(DEcombatseq,file='Revision_Figures/seq_Usoskin/DEcombatseq610_usoskin.rdata')
save(DEraw,file='Revision_Figures/seq_Usoskin/DEraw610_usoskin.rdata')
save(DElimma,file='Revision_Figures/seq_Usoskin/DElimma610_usoskin.rdata')
save(DErescale,file='Revision_Figures/seq_Usoskin/DErescale610_usoskin.rdata')



library(eulerr)
V = list()
eulerr_options(legend=list(font=2,byrow=F,side='bottom',nrow=2,ncol=4))
#library(gridExtra)
for (i in 1:6){
  groups = combn(c('NP','NF','TH','PEP'),2)[,i]
  Vstem <- euler(list(Uncorrected = as.vector(DEraw[[i]]), ComBat = as.vector(DEcombat[[i]]),"ComBat-seq" = as.vector(DEcombatseq[[i]]),
                      MNN = as.vector(DEmnn[[i]]),scBatch = as.vector(DEscbatch[[i]]),
                      limma = as.vector(DElimma[[i]]), RescaleBatches = as.vector(DErescale[[i]])))
  
  if (i %in% c(5)){
    V[[i]] <- plot(Vstem,fills=brewer.pal(n = 8, name = 'Set3')[1:7],quantities=T,lty=1:6,label=F,main=paste(groups[1],'vs',groups[2],sep=' '),legend=list(nrow=4,ncol=2,fontsize=15,font=1))
  }else{
    V[[i]] <- plot(Vstem,fills=brewer.pal(n = 8, name = 'Set3')[1:7],quantities=T,lty=1:6,label=F,main=paste(groups[1],'vs',groups[2],sep=' '),legend=F)
  }
}

grid.arrange(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],nrow=2)
ggsave("Revision_Figures/seq_Usoskin/usoskin_venn.png",plot=grid.arrange(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],nrow=2)
       ,device="png",height=8,width=16,dpi=300)
ggsave("Revision_Figures/seq_Usoskin/usoskin_venn.eps",plot=grid.arrange(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]],nrow=2)
       ,device="eps",height=8,width=16,dpi=300)


for(i in 1:6){
  ggsave(paste("Revision_Figures/seq_Usoskin/usoskin_venn_",i,".eps",sep=""),plot=V[[i]],device='eps',height=4,width=4,dpi=300)
}

ggsave(paste("Revision_Figures/seq_Usoskin/usoskin_venn_",5,".eps",sep=""),plot=V[[5]],device='eps',height=6,width=4,dpi=300)


########################################################################
# Gene Ontology analysis using GOstats

# First save DE gene list for the cell type pairs
# Each file contains detected DE genes from Uncorrected, MNN, ComBat and scBatch
# Please consider save data in a new directory in convenience of next steps
for (i in 1:6){
  group = combn(c('NP','NF','TH','PEP'),2)[,i]
  DElist<-list(Uncorrected = as.vector(DEraw[[i]]), ComBat = as.vector(DEcombat[[i]]),"ComBat-seq" = as.vector(DEcombatseq[[i]]), 
               MNN = as.vector(DEmnn[[i]]),scBatch = as.vector(DEscbatch[[i]]),
               limma = as.vector(DElimma[[i]]), RescaleBatches = as.vector(DErescale[[i]]))
  save(DElist,file=paste('Revision_Figures/seq_Usoskin/DElist/seq_Usoskin610_',group[1],'_',group[2],'.rdata',sep=''))
}

library(org.Hs.eg.db)
library(GOstats)

# Match gene names from raw data with database
all.genes<-rownames(rawdat)
all.entrez<-mget(toupper(all.genes), org.Hs.egSYMBOL2EG, ifnotfound=NA)
all.entrez.2<-all.genes
for(i in 1:length(all.entrez.2)) all.entrez.2[i]<-all.entrez[[i]][1]
all.entrez<-all.entrez.2

# Main function to conduct GO analysis
get.table<-function(sel.genes, all.genes, all.entrez)
{
  sel.entrez<-all.entrez[which(all.genes %in% sel.genes)]
  params <- new("GOHyperGParams", geneIds=sel.entrez[!is.na(sel.entrez)], universeGeneIds=all.entrez[!is.na(all.entrez)], ontology="BP", pvalueCutoff=0.005,conditional=F, testDirection="over", annotation="org.Hs.eg.db")
  over = hyperGTest(params)
  ov<-summary(over)
  ov<-ov[ov[,6]<=500 & ov[,6]>=10,]
  for(i in 2:4) ov[,i]<-signif(ov[,i], 3)
  ov
}

setwd('Revision_Figures/seq_Usoskin/DElist')
# Obtain the file names of saved DE gene lists (should have 6 .rdata files)
files<-dir(pattern=".rdata")

# Conduct and save GO analysis results for all 6 cell type pairs
for(i in 1:6)
{
  # Read saved data
  load(files[i])
  GOlist<-new("list")
  for(k in 1:7)
  {
    GOlist[[k]]<-get.table(DElist[[k]], all.genes, all.entrez)
  }
  names(GOlist)<-names(DElist)
  
  b<-GOlist[[1]]
  b<-rbind(c(names(GOlist)[1],rep("",6)), b)
  for(k in 2:7)
  {
    b<-rbind(b, c(names(GOlist)[k],rep("",6)))
    b<-rbind(b, GOlist[[k]])
  }
  
  write.table(b, paste(strsplit(files[i], ".rdata")[[1]][1],".GO 005.txt"), sep="\t", quote=F, col.names=F, row.names=F)
}
