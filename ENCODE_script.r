library(scBatch)
library(pheatmap)
library(scran)
library(sva)
library(mclust)

data("ENCODE")
batches <- c(rep(1,13),rep(2,13))
cell.type <- c(1:13,1:13)

#raw data heatmap
pheatmap(cor(ENCODE))

#standardization

standardization <- function(dat,batch)
{
  if(is.matrix(dat)==FALSE){
    dat <- as.matrix(dat)
  }
  
  for (i in 1:dim(dat)[1]){
    for (j in unique(batch)){
      mean <- mean(dat[i,batch==j])
      var <- var(dat[i,batch==j])
      dat[i,batch==j] <- (dat[i,batch==j] - mean)/sqrt(var)
    }
  }
  dat <- stats::na.omit(dat)
  return(dat)
}

dat <- standardization(ENCODE,batches)

#ComBat

combatmod = ComBat(dat,batches)
pheatmap(cor(combatmod))

#MNN

batch1 = dat[,1:13]
batch2 = dat[,14:26]
mnnmod = mnnCorrect(batch1,batch2,k=13)
mnnmod = cbind(mnnmod$corrected[[1]],mnnmod$corrected[[2]])
colnames(mnnmod) = colnames(ENCODE)
pheatmap(cor(mnnmod))

#scBatch

set.seed(1234)
distmod <- QuantNorm(dat,batches,method='row/column', cor_method='pearson', logdat=F,tol=1e-4)
scbatchmod <- scBatchCpp(c=dat,w=diag(26),d=distmod,m=1,max=1000,step=1e-6,tol=1e-20,derif=scBatch::derif,verbose = T)
colnames(scbatchmod) = colnames(ENCODE)
pheatmap(cor(scbatchmod))

#Calculate ARIs

ARI <- function(cell.type,ccc,k){
  hc <- hclust(dist(ccc),method = 'complete')
  pred <- cutree(hc,k)
  ARI.complete <- adjustedRandIndex(cell.type,pred)
  
  hc <- hclust(dist(ccc),method = 'average')
  pred <- cutree(hc,k)
  ARI.average <- adjustedRandIndex(cell.type,pred)
  
  hc <- hclust(dist(ccc),method = 'single')
  pred <- cutree(hc,k)
  ARI.single <- adjustedRandIndex(cell.type,pred)
  
  return(max(ARI.complete,ARI.average,ARI.single))
}

ARI(cell.type,cor(scbatchmod),13)
#0.884
ARI(cell.type,cor(combatmod),13)
#0.701
ARI(cell.type,cor(mnnmod),13)
#0.675
ARI(cell.type,cor(ENCODE),13)
#-0.056