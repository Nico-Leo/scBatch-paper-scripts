# The following script generates AUCs and PRAUCs used in Fig. 2(IV), where the batch location and scale parameters are varied among batches.
# The scripts for the other 4 configurations can be easily obtained based on the instructions in this script.

# Using parallel computing for simulation

library(doSNOW)
n.node<-25
cl<-makeSOCKcluster(rep("localhost",n.node))
registerDoSNOW(cl)

trueDE <- function(object,group1,group2){
  #Output DE genes, defined by > 1.5 logFC, for the two target groups
 
  output <- list()
  rowdata <- rowData(object)
  fc <- rowdata[,group1]/rowdata[,group2]
  fc[fc < 1] <- 1/fc[fc < 1]
  output$DEidx <- as.numeric(fc > 1.5)
  output$DE <- as.character(rowdata$Gene)[output$DEidx==1]
  return(output)
}


DEanalysis <- function(dat,cell.type,DE,groups){
  #Conduct Seurat DE analysis to find DE genes and return AUC and PRAUC based on truth
  
  require(PRROC)

  colnames(dat) <- paste0(cell.type, "__", 1:ncol(dat))
  d = CreateSeuratObject(raw.data = (dat-min(dat))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     ident.2 = levels(factor(groups))[2], test.use = "bimod")
  pval.edgeR = res$p_val
  names(pval.edgeR) = rownames(res)
  q = p.adjust(res$p_val, method = "BH")

  qAUC <- roc.curve(1-q[names(pval.edgeR) %in% DE$DE],1-q[!(names(pval.edgeR) %in% DE$DE)])$auc
  qPRAUC <- pr.curve(1-q[names(pval.edgeR) %in% DE$DE],1-q[!(names(pval.edgeR) %in% DE$DE)])$auc.integral
 
  return(c(qAUC,qPRAUC))
}

library(scBatch)
library(Seurat)
library(scran)
library(splatter)
library(QuantNorm)
library(sva)
library(mclust)
library(PRROC)

#Sample size
n = 1000

#Number of genes
p = 10000

seeds <- sample(1:10000,50)

# Using parallel computing to collect results
results <- foreach(i=1:50,.combine="rbind",.verbose=TRUE,.packages=c('splatter','QuantNorm','sva','mclust','scran','scBatch','Seurat')) %dopar%{
      
      sim.groups <- splatSimulate(nGenes=p, batchCells = c(n/4,n/4,n/4,n/4), seed = seeds[i],group.prob = c(0.4,0.3,0.2,0.1),
                                  de.prob = c(0.1,0.1,0.1,0.1), batch.facLoc=c(0.1,0.2,0.05,0.15), 
                                  batch.facScale=c(0.05,0.1,0.25,0.3), method = "groups", verbose = F)
                                  
      #For configurations (I)-(III), we can let batch.facScale or batch.facLoc undefined in splatSimulate function.
      #For configuration (V), we need to modify batchCells = n since there is no batch effects. In addition, no batch.facScale and batch.facLoc should be specified.
     
      # Normalize the count matrix
      sim.groups <- normalise(sim.groups)
      
      # Find true DE genes for the six pairs of cell types
      DE12 <- trueDE(sim.groups,'DEFacGroup1','DEFacGroup2')
      DE13 <- trueDE(sim.groups,'DEFacGroup1','DEFacGroup3')
      DE14 <- trueDE(sim.groups,'DEFacGroup1','DEFacGroup4')
      DE23 <- trueDE(sim.groups,'DEFacGroup2','DEFacGroup3')
      DE24 <- trueDE(sim.groups,'DEFacGroup2','DEFacGroup4')
      DE34 <- trueDE(sim.groups,'DEFacGroup3','DEFacGroup4')
      
      # Extract batch information
      batch = sim.groups@colData$Batch
      # In configuration (V), we randomly assign batches using batch = sample(1:4,n,replace = T,prob = c(0.25,0.25,0.25,0.25))

      # Extract cell types
      cell.type = sim.groups@colData$Group
      
      # Extract count matrix
      exp <- exprs(sim.groups)

      # Obtain the modified distance matrix by QuantNorm
      ccc <- QuantNorm(exp,as.numeric(as.factor(batch)),logdat=F,method='row/column',cor_method='pearson',max=5)

      # Conduct scBatch to obtain modified count matrix
      GD5 <-scBatchCpp(c=exp,d=ccc,w=diag(n),m=5,max=1200,tol=1e-10,step=0.0001,derif=scBatch::derif,verbose=F)
      colnames(GD5) = colnames(exp)
      rownames(GD5) = rownames(exp)

      # Obtain MNN correction
      batches = as.numeric(as.factor(batch))
      mnn.out <- mnnCorrect(exp[,batches==1],exp[,batches==2],exp[,batches==3],exp[,batches==4],k=20,cos.norm.in=T,cos.norm.out=F)
      X.mnn<-do.call(cbind, mnn.out$corrected)
      colnames(X.mnn) = colnames(exp)

      # Obtain ComBat correction
      cleandat <- ComBat(exp[rowSums(exp)>0,],batch)

      # Conduct DE analysis for different cell types and report AUC and PRAUC
      c(DEanalysis(exp,cell.type,DE12,c('Group1','Group2')),DEanalysis(GD5,cell.type,DE12,c('Group1','Group2')),
        DEanalysis(X.mnn,cell.type,DE12,c('Group1','Group2')),DEanalysis(cleandat,cell.type,DE12,c('Group1','Group2')),
        DEanalysis(exp,cell.type,DE13,c('Group1','Group3')),DEanalysis(GD5,cell.type,DE13,c('Group1','Group3')),
        DEanalysis(X.mnn,cell.type,DE13,c('Group1','Group3')),DEanalysis(cleandat,cell.type,DE13,c('Group1','Group3')),
        DEanalysis(exp,cell.type,DE14,c('Group1','Group4')),DEanalysis(GD5,cell.type,DE14,c('Group1','Group4')),
        DEanalysis(X.mnn,cell.type,DE14,c('Group1','Group4')),DEanalysis(cleandat,cell.type,DE14,c('Group1','Group4')),
        DEanalysis(exp,cell.type,DE23,c('Group2','Group3')),DEanalysis(GD5,cell.type,DE23,c('Group2','Group3')),
        DEanalysis(X.mnn,cell.type,DE23,c('Group2','Group3')),DEanalysis(cleandat,cell.type,DE23,c('Group2','Group3')),
        DEanalysis(exp,cell.type,DE24,c('Group2','Group4')),DEanalysis(GD5,cell.type,DE24,c('Group2','Group4')),
        DEanalysis(X.mnn,cell.type,DE24,c('Group2','Group4')),DEanalysis(cleandat,cell.type,DE24,c('Group2','Group4')),
        DEanalysis(exp,cell.type,DE34,c('Group3','Group4')),DEanalysis(GD5,cell.type,DE34,c('Group3','Group4')),
        DEanalysis(X.mnn,cell.type,DE34,c('Group3','Group4')),DEanalysis(cleandat,cell.type,DE34,c('Group3','Group4')))
}

# The script to plot boxplots are omitted
