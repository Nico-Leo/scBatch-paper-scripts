#####################################################

## functions for simulating RNA-seq count data.

## Take an option object and sample sizes.

#####################################################



simRNAseq1 <- function(simOptions, n1, n2, DEid) {
  
  ##  set.seed(simOptions$sim.seed)
  
  if(simOptions$design == "2grp")
    
    data <- simRNAseq.2grp1(simOptions, n1, n2, DEid)
  
  
  
  data
  
}



#####################################################

## simulate 2-group RNA-seq data

#####################################################

simRNAseq.2grp1 <- function(simOptions, n1, n2, DEid) {
  
  
  
  ## make design vector
  
  design <- c(rep(0, n1), rep(1, n2))
  
  
  
  ## generate ID for DE genes
  
  ngenes = simOptions$ngenes
  
  p.DE = simOptions$p.DE
  
  
  
  ## generate lfc for all genes
  
  lfc = rep(0, ngenes)
  
  lfc[DEid] = simOptions$lfc
  
  
  
  ## generate mean expressions for all replicates
  
  lmeanExpr <- makeMeanExpr.2grp(simOptions$lBaselineExpr, lfc, n1, n2)
  
  
  
  ## generate counts
  
  allmu <- lmeanExpr
  
  phi <- exp(simOptions$lOD)
  
  x <- rnegbinom(length(allmu), allmu, phi)
  
  x <- matrix(x,ncol=n1+n2)
  
  
  
  ## return
  
  list(counts=x, designs=design, DEid=DEid, simOptions=simOptions)
  
}



###############################################################

## generate mean expression values in two group comparison.

###############################################################

makeMeanExpr.2grp <- function(lBaselineExpr, lfc, n1, n2) {
  
  ## generate mean expressions in two groups
  
  result <- matrix(exp(lBaselineExpr+c(rep(lfc/2,n1),rep(-lfc/2,n2))),
                   
                   ncol=n1+n2)
  
  ##result <- sweep(result, 2, sizefactor, FUN="*")
  
  result
  
}





###############################################################

## run simulation and DE detection

###############################################################

runSims <- function(Nreps=c(3,5,7,10), Nreps2, nsims=100, sim.opts, DEmethod=c("edgeR", "DSS", "DESeq"), verbose=TRUE) {
  
  
  
  DEmethod = match.arg(DEmethod)
  
  ## generate size factor if not given
  
  ##   if(missing(sizefactor)) {
  
  ##     sizefactor=rep(1, n1+n2)
  
  ##   } else if(!is.vector(sizefactor) | length(sizefactor) != (m.n1+m.n2) ) {
  
  ##     stop("sizefactor must be a vector of length n1+n2!\n")
  
  ##   }
  
  
  
  if(missing(Nreps2))
    
    Nreps2 = Nreps
  
  else {
    
    if(length(Nreps2) != length(Nreps))
      
      stop("Nreps and Nreps2 must be vectors of the same length")
    
  }
  
  n1 = max(Nreps)
  
  n2 = max(Nreps2)
  
  
  
  ## start simulation
  
  set.seed(sim.opts$sim.seed)
  
  pvalue = fdrs = xbar = array(NA,dim=c(sim.opts$ngenes,length(Nreps), nsims))
  
  DEids = lfcs = NULL
  
  for(i in 1:nsims) {
    
    if(verbose)
      
      cat("Simulation number", i, "\n")
    
    
    
    ## update the simulatino option. Reregenerate the DEid and lfc
    
    sim.opts$sim.seed = sim.opts$sim.seed + 1
    
    sim.opts=update.RNAseq.SimOptions.2grp(sim.opts)
    
    
    
    ## generate data
    
    ##sim.opts$sim.seed = sim.opts$sim.seed + 1
    
    dat.sim.big = simRNAseq(sim.opts, n1, n2)
    
    DEids[[i]] = dat.sim.big$DEid
    
    lfcs[[i]] = dat.sim.big$simOptions$lfc
    
    
    
    ##  for different sample sizes
    
    for(j in seq(along=Nreps)) {
      
      nn1 = Nreps[j]
      
      nn2 = Nreps2[j]
      
      
      
      ## take a subsample of the simulated counts
      
      idx = c(1:nn1, n1+(1:nn2))
      
      this.design = dat.sim.big$designs[idx]
      
      this.X = dat.sim.big$counts[,idx]
      
      this.simOpts = sim.opts
      
      ## filter out genes with all 0 counts
      
      ss = rowSums(this.X)
      
      ix.valid = ss>0
      
      this.X.valid = this.X[ix.valid,, drop=FALSE]
      
      
      
      ## create an object and pass into DE detection
      
      data0=list(counts=this.X.valid, designs=this.design)
      
      if (DEmethod=="edgeR")
        
        res1=run.edgeR(data0)
      
      if (DEmethod=="DESeq")
        
        res1=run.DESeq(data0)
      
      if (DEmethod=="DSS")
        
        res1=run.DSS(data0)
      
      
      
      ## store results. Be careful here about filtering
      
      pval = fdr = rep(1, nrow(this.X))
      
      X.bar1 = rep(0, nrow(this.X))
      
      pval[ix.valid] = res1[, "pval"]
      
      fdr[ix.valid] = res1[, "fdr"]
      
      sizeF = colSums(data0$count); sizeF = sizeF/median(sizeF)#size factor
      
      X.bar1[ix.valid] = rowMeans(sweep(data0$count,2,sizeF,FUN="/"))
      
      pvalue[,j,i] = pval
      
      fdrs[,j,i] = fdr
      
      xbar[,j,i]=X.bar1
      
    }
    
  }
  
  
  
  ## return
  
  list(pvalue=pvalue, fdrs=fdrs, xbar=xbar, DEid=DEids, lfcs=lfcs, Nreps1=Nreps, Nreps2=Nreps2, sim.opts=sim.opts)
  
}





###########################################################################

## update the simulation option object.

## This is to regenreate a new set of DE genes.

## They have different (regenreated) parameters.

###########################################################################

update.RNAseq.SimOptions.2grp <- function(sim.opts) {
  
  ## update DEid and lfc, but keep the others
  
  RNAseq.SimOptions.2grp(ngenes = sim.opts$ngenes,
                         
                         lBaselineExpr=sim.opts$lBaselineExpr,
                         
                         lOD=sim.opts$lOD,
                         
                         p.DE=sim.opts$p.DE,lfc=sim.opts$lfc,
                         
                         sim.seed=sim.opts$sim.seed)
  
}

######### some utility functions for simulation



######################################################################

## function to set baseline expression

######################################################################

setBaselineExpr <- function(input, ngenes) {
  
  
  
  param = NULL
  
  if(is.numeric(input)) { ## numeric
    
    if(length(input)==1 & is.numeric(input)) { ## constant
      
      lmeanExpr=rep(input, ngenes)
      
    } else if (length(input)!=ngenes) { ## vector
      
      #            stop("The length of lmeanExpr doesn't equal to ngenes!\n")
      
      lmeanExpr=sample(input, ngenes, replace=TRUE)
      
    } else
      
      lmeanExpr = input
    
  } else if (is.function(input)) { # a function
    
    lmeanExpr=input(ngenes)
    
  } else if (is.character(input)) { ## a character. Sample from real data.
    
    ## Note I save the mean expr in log-scale so have to exp them.
    
    datatsets <- c("cheung", "gilad","maqc","bottomly")
    
    if(input %in% datatsets ) {
      
      eval(parse(text=paste0("data(",input,", envir=environment())")))
      
    } else {
      
      stop("Unrecognized string of lmeanExpr. It must be one of 'cheung', 'gilad', 'bottomly', or 'maqc'!\n")
      
    }
    
    lmeanExpr=sample(param$lmean, ngenes, replace=TRUE)
    
  }
  
  else {
    
    stop("Unrecognized form of lmeanExpr!\n")
    
  }
  
  lmeanExpr
  
}



######################################################################

## function to set baseline expression, given sequencing depth

######################################################################

setBaselineExpr.seqDepth <- function(seqDepth, ngenes) {
  
  GE.human = NULL
  
  data(GE.human, envir=environment())
  
  GE.sample = sample(GE.human, ngenes, replace=TRUE)
  
  ntotal = sum(GE.sample)
  
  p0 = seqDepth / ntotal
  
  lmeanExpr=log(GE.sample*p0)
  
  lmeanExpr
  
}



######################################################################

## function to set over dispersion parameter

######################################################################

setOD <- function(input, ngenes) {
  
  param <- NULL
  
  if(is.numeric(input)) {
    
    if(length(input)==1) { ## constant
      
      lOD=rep(input, ngenes)
      
    } else if (length(input)!=ngenes) { ## vector
      
      #stop("The length of OD doesn't equal to ngenes!\n")
      
      lOD=sample(input, ngenes, replace=TRUE)
      
    } else
      
      lOD = input
    
  } else if (is.function(input)) { # a function
    
    lOD=input(ngenes)
    
  } else if (is.character(input)) { ## a character
    
    datatsets <- c("cheung", "gilad","maqc","bottomly")
    
    if(input %in% datatsets ) {
      
      eval(parse(text=paste0("data(",input,", envir=environment())")))
      
    } else {
      
      stop("Unrecognized string of lmeanExpr. It must be one of 'cheung', 'gilad', 'bottomly', or 'maqc'!\n")
      
    }
    
    lOD=sample(param$lOD, ngenes, replace=TRUE)
    
  }
  
  else {
    
    stop("Unrecognized form of lOD!\n")
    
  }
  
  lOD
  
}





##############################################################

## default function to generate null and

## alternative log fold change

##############################################################

lfc.null <- function(n0) {
  
  rnorm(n0, mean=0, sd=0.0)
  
}



lfc.alt <- function(nDE) {
  
  nDE1 = round(nDE/2); nDE2 = nDE - nDE1
  
  ##lfc = c(rnorm(nDE1, mean=1, sd=0.2), rnorm(nDE2, mean=-1, sd=0.2))
  
  ##lfc = c(runif(nDE1, 0.5, 3), runif(nDE2, -3, -0.5))
  
  lfc = rnorm(nDE, 0, 1.5)
  
  lfc
  
}





######################################################################

## generate negative binomial rv, given mu and phi (over-dispersion)

######################################################################

rnegbinom <- function (n, mu=1, phi=0.01){
  
  rpois(n, rgamma(n, shape=1/phi,scale=mu*phi))
  
}

#Please run the blocks ONE and TWO first for some functions and packages.
######################################################################################################
#BLOCK ONE#
###########


library(doSNOW)

n.node<-50

cl<-makeSOCKcluster(rep("localhost",n.node))

registerDoSNOW(cl)

DEanalysis <- function(dat,cell.type,id,groups,scaledata=F){
  
  #Conduct Seurat DE analysis to find DE genes and return AUC and PRAUC based on truth
  require(PRROC)
  require(Seurat)
  
  colnames(dat) <- paste0(cell.type, "__", 1:ncol(dat))
  
  d = CreateSeuratObject(counts = (dat-min(dat))[,cell.type %in% groups],names.field = 1, names.delim = "__")
  
  if (scaledata == T){
    d = NormalizeData(d)
  }
  
  res <- FindMarkers(d, ident.1 = levels(factor(groups))[1],
                     
                     ident.2 = levels(factor(groups))[2], test.use = "wilcox")#,min.pct = 0,logfc.threshold = 0)
  
  pval.edgeR = res$p_val
  
  names(pval.edgeR) = rownames(res)
  
  pval.edgeR = pval.edgeR[!is.na(pval.edgeR)]
  
  q = p.adjust(pval.edgeR, method = "BH")
  
  qAUC <- roc.curve(1-q[names(q) %in% id],1-q[!(names(q) %in% id)])$auc
  qPRAUC <- pr.curve(1-q[names(q) %in% id],1-q[!(names(q) %in% id)])$auc.integral
  
  signif <- q[q<0.2]
  nDE <- sum(names(signif) %in% id)
  return(c(qAUC,qPRAUC,nDE))
  
}

CompARI <- function(dat,cell.type){
  ccc <- 1-cor(dat,method='pearson')
  require(mclust)
  ARI=0
  for (i in 1:50){
    kms <- kmeans(ccc,6)
    ARI = ARI + adjustedRandIndex(kms$cluster,cell.type)
  }
  ARI/50
}


#Number of genes

p=20000

seeds <- sample(1:10000,100)

#sample size
for (size in c(270,540,1080)){
  
  #de probability
  for (de in c(0.05)){
    
    #sequencing depth
    for (depth in c(0.1)){
      
      results <- foreach(i=1:50,.combine="rbind",.multicombine = T,.maxcombine = 100,.verbose=TRUE,.packages=c('PROPER','sva','batchelor','limma','scBatch','Seurat','kBET','scater')) %dopar%{
       
        #Data simulation 
        sim.opts.Cheung1 = RNAseq.SimOptions.2grp(ngenes = p, p.DE=de,
                                                  lOD="cheung", lBaselineExpr="cheung",sim.seed = seeds[i])
        
        N <- sum(exp(sim.opts.Cheung1$lBaselineExpr))
        deps <- floor(depth*N)
        sim.opts.Cheung1$lBaselineExpr <- as.vector(log(rmultinom(1,deps,exp(sim.opts.Cheung1$lBaselineExpr)/N)))
        sim.opts.Cheung1$lOD <- sim.opts.Cheung1$lOD + runif(p,0,0.5)
        
        sim1 <- simRNAseq(sim.opts.Cheung1, n1=30, n2=90)
        sim2 <- simRNAseq(sim.opts.Cheung1, n1=30, n2=30)
        sim3 <- simRNAseq(sim.opts.Cheung1, n1=30, n2=150)
        
        id1 <- sim1$DEid
        id2 <- sim2$DEid
        id3 <- sim3$DEid
        
        a1 = sim1$counts
        a2 = sim2$counts
        a3 = sim3$counts
        
        sim.opts.Cheung1b <- sim.opts.Cheung1
        sim.opts.Cheung1b$lOD <- sim.opts.Cheung1$lOD + runif(p,0,0.5)
        sim.opts.Cheung1b$lBaselineExpr <- sim.opts.Cheung1$lBaselineExpr + runif(p,-3,3)
        
        
        b1 = simRNAseq1(sim.opts.Cheung1b, n1=30, n2=100, id1)$counts
        b2 = simRNAseq1(sim.opts.Cheung1b, n1=30, n2=40, id2)$counts
        b3 = simRNAseq1(sim.opts.Cheung1b, n1=30, n2=130, id3)$counts
        
        sim.opts.Cheung1c <- sim.opts.Cheung1
        sim.opts.Cheung1c$lOD <- sim.opts.Cheung1$lOD + runif(p,4,5)
        sim.opts.Cheung1c$lBaselineExpr <- sim.opts.Cheung1$lBaselineExpr + runif(p,-3,3)
        
        
        c1 = simRNAseq1(sim.opts.Cheung1c, n1=30, n2=85, id1)$counts
        c2 = simRNAseq1(sim.opts.Cheung1c, n1=35, n2=40, id2)$counts
        c3 = simRNAseq1(sim.opts.Cheung1c, n1=30, n2=140, id3)$counts
        
        bat1 = cbind(a1,a2,a3)
        bat2 = cbind(b1,b2,b3)
        bat3 = cbind(c1,c2,c3)
        

        batch <- cbind(bat1,bat2,bat3)
        rownames(batch) <- 1:nrow(batch)
        
        batch <- batch[rowSums(batch) > 0,]
        colnames(batch) <- c(rep(1,30),rep(2,90),rep(3,30),rep(4,30),rep(5,30),rep(6,150),rep(1,30),rep(2,100),rep(3,30),rep(4,40),rep(5,60),rep(6,100),rep(1,30),rep(2,85),rep(3,35),rep(4,40),rep(5,60),rep(6,110))
        
        bat1 <- batch[,1:360] 
        bat2 <- batch[,361:720]
        bat3 <- batch[,721:1080]
        
        randnum <- sample(1:dim(bat1)[2],size/3,replace=FALSE)
        batch1 <- bat1[,randnum]
        
        randnum <- sample(1:dim(bat2)[2],size/3,replace=FALSE)
        batch2 <- bat2[,randnum]
        
        randnum <- sample(1:dim(bat3)[2],size/3,replace=FALSE)
        batch3 <- bat3[,randnum]
        
        newbatch <- cbind(batch1[,sample(ncol(batch1))],batch2[,sample(ncol(batch2))],batch3[,sample(ncol(batch3))])
        newbatch <- newbatch[complete.cases(newbatch),]
        batches <- c(rep(1,dim(batch1)[2]),rep(2,dim(batch2)[2]),rep(3,dim(batch3)[2]))
        
        raw = log2(newbatch+1)
        raw = raw[rowSums(raw) > 0,]
        
        # Normalize the count matrix
        
        colnames(newbatch) = 1:ncol(newbatch)
        sce <- SingleCellExperiment(assays=list(counts=newbatch))
        sce <- normalize(sce)
        exp <- exprs(sce)
        
        # Obtain the modified distance matrix by QuantNorm
        start_time <- proc.time()
        ccc <- scBatch::QuantNorm(exp,as.numeric(as.factor(batches)),logdat=F,method='row/column',cor_method='pearson',max=10)
        
        # Conduct scBatch to obtain modified count matrix
        
        GD1 <-scBatchCpp(c=exp,d=ccc,w=diag(size),m=1,max=20,tol=1e-8,step=0.0001,derif=scBatch::derif,verbose=F)
        GD1new <- GD1
        for (k in 1:nrow(GD1)){
          GD1new[k,batches==1] <- (GD1[k,batches==1] - mean(GD1[k,batches==1]))/sqrt(var(GD1[k,batches==1]))
          GD1new[k,batches==2] <- (GD1[k,batches==2] - mean(GD1[k,batches==2]))/sqrt(var(GD1[k,batches==2]))
          GD1new[k,batches==3] <- (GD1[k,batches==3] - mean(GD1[k,batches==3]))/sqrt(var(GD1[k,batches==3]))
        }
        colnames(GD1new) = colnames(raw)
        rownames(GD1new) = rownames(raw)
        end_time <- proc.time()
        timescbatch <- (end_time - start_time)[3]
        
        # Obtain MNN correction
        
        start_time <- proc.time()
        mnn.out <- mnnCorrect(sce,batch=as.factor(batches),k=20,cos.norm.in=T,cos.norm.out=F)
        X.mnn<-mnn.out@assays@.xData$data$corrected
        colnames(X.mnn) = colnames(raw)
        rownames(X.mnn) = rownames(raw)
        end_time <- proc.time()
        timemnn <- (end_time - start_time)[3]
        
        # Obtain ComBat correction
        
        start_time <- proc.time()
        cleandat <- ComBat(exp[rowSums(exp)>0,],batches)
        end_time <- proc.time()
        timecombat <- (end_time - start_time)[3]
        
        # Obtain limma correction
        start_time <- proc.time()
        limmadat <- removeBatchEffect(exp,batches)
        end_time <- proc.time()
        timelimma <- (end_time - start_time)[3]
        
        # Obtain batchelor correction
        start_time <- proc.time()
        batchelordat <- rescaleBatches(sce,batch=as.factor(batches))  
        batchelordat <- batchelordat@assays@.xData$data$corrected
        colnames(batchelordat) = colnames(exp)
        rownames(batchelordat) = rownames(exp)
        end_time <- proc.time()
        timerescale <- (end_time - start_time)[3]
        
        # Conduct DE analysis for different cell types and report AUC and PRAUC
        
        cell.type = colnames(raw)
        
        c(DEanalysis(raw,cell.type,id1,c(1,2)), DEanalysis(exp,cell.type,id1,c(1,2)),DEanalysis(GD1new,cell.type,id1,c(1,2)),
          
          DEanalysis(X.mnn,cell.type,id1,c(1,2)),DEanalysis(cleandat,cell.type,id1,c(1,2)),
          
          DEanalysis(limmadat,cell.type,id1,c(1,2)),DEanalysis(batchelordat,cell.type,id1,c(1,2)),
          
          DEanalysis(raw,cell.type,id2,c(3,4)), DEanalysis(exp,cell.type,id2,c(3,4)),DEanalysis(GD1new,cell.type,id2,c(3,4)),
          
          DEanalysis(X.mnn,cell.type,id2,c(3,4)),DEanalysis(cleandat,cell.type,id2,c(3,4)),
          
          DEanalysis(limmadat,cell.type,id2,c(3,4)),DEanalysis(batchelordat,cell.type,id2,c(3,4)),
          
          DEanalysis(raw,cell.type,id3,c(5,6)), DEanalysis(exp,cell.type,id3,c(5,6)),DEanalysis(GD1new,cell.type,id3,c(5,6)),
          
          DEanalysis(X.mnn,cell.type,id3,c(5,6)),DEanalysis(cleandat,cell.type,id3,c(5,6)),
          
          DEanalysis(limmadat,cell.type,id3,c(5,6)),DEanalysis(batchelordat,cell.type,id3,c(5,6)),
          
          CompARI(raw,cell.type),CompARI(exp,cell.type),CompARI(GD1new,cell.type),
          CompARI(X.mnn,cell.type),CompARI(cleandat,cell.type),CompARI(limmadat,cell.type),
          CompARI(batchelordat,cell.type),
          
          timescbatch,timemnn,timecombat,timelimma,timerescale)
      }
      
      save(results,file=paste(p,"_",size,"_",de,"_",depth,".rdata",sep=""),version=2)
      
    }
  }
}