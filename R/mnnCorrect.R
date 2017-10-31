mnnCorrect<- function(..., inquiry.genes=NULL, hvg.genes=NULL, k=20, sigma=1, cos.norm=TRUE, svd.dim=0, order=NULL,k.clara=0, withQC=FALSE,varCare=TRUE) 
# Performs correction based on the batches specified in the ellipsis.
#    
# written by Laleh Haghverdi
# with modifications by Aaron Lun
# created 7 April 2017
# last modified 22 June 2017
{
  #library(FNN)
  library(cluster)
  #library(scran)
  library(Matrix)
  library(matrixStats)
  
    batches <- batches0 <- list(...) 
    nbatches <- length(batches) 
    if (nbatches < 2L) { stop("at least two batches must be specified") }
   
    # Checking for identical number of rows (and rownames).
    first <- batches[[1]]
    ref.nrow <- nrow(first)
    ref.rownames <- rownames(first)
    for (b in 2:nbatches) {
        current <- batches[[b]]
        if (!identical(nrow(current), ref.nrow)) {
            stop("number of rows is not the same across batches")
        } else if (!identical(rownames(current), ref.rownames)) {
            stop("row names are not the same across batches")
        }
    }

    if (k.clara<k & k.clara!=0) {
      stop("no. pre clusters should either be zero or larger than k")}
      
    # Subsetting to the desired subset of genes.
    if (!is.null(inquiry.genes)) {
        batches0 <- lapply(batches0, "[", i=inquiry.genes, , drop=FALSE) # Need the extra comma!
    }
    if (!is.null(hvg.genes)) { 
        batches <- lapply(batches, "[", i=hvg.genes, , drop=FALSE)
    }
    inquiry.genes <- scran:::.subset_to_index(inquiry.genes, first, byrow=TRUE)
    hvg.genes <- scran:::.subset_to_index(hvg.genes, first, byrow=TRUE)
    inquiry.in.hvg <- inquiry.genes %in% hvg.genes 

    # Applying cosine normalization for MNN identification. 
    if (cos.norm) { batches <- lapply(batches, cosine.norm) }

    # Setting up the order.
    if (is.null(order)) {
        order <- seq_len(nbatches)
    } else {
        order <- as.integer(order)
        if (!identical(seq_len(nbatches), sort(order))) { 
            stop(sprintf("'order' should contain values in 1:%i", nbatches))
        }
    }
   
    # Setting up the variables.
    ref <- order[1]
    ref.batch <- batches[[ref]]    #for hidden calculations
    ref.batch0 <- batches0[[ref]]  #for output G.N
    
    output0 <- vector("list", nbatches)
    output0[[ref]] <- ref.batch0        #output first G.N batch before Gn modifications by aggregate
    
    if ( (ncol(ref.batch)> 3000) & (k.clara>0) ){
    #k.clara=1000#round(ncol(ref.batch)/1000)+1
    X<-clara(t(ref.batch),k=k.clara)
    ref.batch<-t(X$medoids)  #gn
    ref.batch0<-t(aggregate(t(ref.batch0),by=list(X$clustering),FUN="median"))[-1,]  #Gn  ##gives error
    }
    
    #mnns.list <- list(NA_integer_, nbatches, 2)
    mnns.list <- matrix(list(), nrow=nbatches, ncol=2)
    
    print(paste0("The first batch is taken as the first reference."))
    
    for (b in 2:nbatches) { 
        print(paste0("Mergeing batch ", b ,"..." ))
      
        other.batch <- batches[[order[b]]] #hidden layer  gN
        other.batch0 <- batches0[[order[b]]] #makes output GN
        
        if ( (ncol(other.batch)>3000) & (k.clara>0) ){
        #k.clara=1000#round(ncol(other.batch)/1000)+1
        X<-clara(t(other.batch),k=k.clara)
        other.batch<-t(X$medoids)           #gn
        other.clust<-X$clustering
        }
        else {other.clust<-NULL#seq_len(ncol(other.batch))
        }
        # Finding pairs of mutual nearest neighbours.
        sets <- find.mutual.nn(ref.batch, other.batch, ref.batch0, other.batch0, clust2=other.clust, k1=k, k2=k, sigma=sigma, withQC=withQC)
        #on the dimension of result: list(gn.correction.vector , GN.correction0.vector)<-find.mutual.nn(gn,gn,Gn,GN)
        s1 <- sets$set1
        s2 <- sets$set2

        if (svd.dim==0){
            correction <- t(sets$vect)   #gn
            correction0 <- t(sets$vect0) #GN
        } else {
            ## Computing the biological subspace in both batches.
            ndim <- min(c(svd.dim, dim(ref.batch), dim(other.batch)))
            span1 <- get.bio.span(ref.batch0[,s1,drop=FALSE], inquiry.in.hvg, min(ndim, length(s1)))
            span2 <- get.bio.span(other.batch0[,s2,drop=FALSE], inquiry.in.hvg, min(ndim, length(s2)))
            #nshared <- find.shared.subspace(span1, span2, assume.orthonormal=TRUE, get.angle=FALSE)$nshared
            #if (nshared==0L) { warning("batches not sufficiently related") }
    
            #reduce the component in each span from the batch correction vector, span1 span2 order does not matter
            
            bv0 <- sets$vect0      
            bio.comp0 <- bv0 %*% span1 %*% t(span1)
            correction0 <- t(bv0) - t(bio.comp0)
            bio.comp <- t(correction0) %*% span2 %*% t(span2)
            correction0 <- correction0 - t(bio.comp0)
            
            ndim <- min(c(svd.dim, dim(ref.batch), dim(other.batch)))
            hvg.in.hvg<-rep(TRUE,length(hvg.genes)) #hvg.genes %in% hvg.genes
            span1 <- get.bio.span(ref.batch[,s1,drop=FALSE], hvg.in.hvg, min(ndim, length(s1)))
            span2 <- get.bio.span(other.batch[,s2,drop=FALSE], hvg.in.hvg, min(ndim, length(s2)))
            
            bv <- sets$vect
            bio.comp <- bv %*% span1 %*% t(span1)
            correction <- t(bv) - t(bio.comp)
            bio.comp <- t(correction) %*% span2 %*% t(span2)
            correction <- correction - t(bio.comp)
            #correction <- t(sets$vect)
        } 
        
        # Applying the correction and storing the numbers of nearest neighbors.
        other.batch <- other.batch + correction #gn+gn 
        other.batch0 <- other.batch0 + correction0 #GN+GN
        
        #mnns.list[[b]] <- list(s1,s2)
        #mnns.list[b,] <- c(length(s1), length(s2))
        mnns.list[[b,1]] <- s1
        mnns.list[[b,2]] <- s2
        output0[[b]] <- other.batch0 #GN

        # Expanding the reference batch to include the new, corrected data.
        ref.batch <- cbind(ref.batch, other.batch) #gn,gn
        ref.batch0 <- cbind(ref.batch0, other.batch0) #G.n,G.n
    }

    # Formatting output to be consistent with input.
    names(output0) <- names(batches0)
    list(corrected=output0, mnns.list=mnns.list, batch.vects=sets$vect0)
}

find.mutual.nn <- function(exprs1, exprs2, exprs10, exprs20,clust2=NULL, k1, k2, sigma=1, withQC=FALSE,varCare=TRUE)
# Finds mutal neighbors between data1 and data2.
# Computes the batch correction vector for each cell in data2.
{
   
  if (!is.null(clust2)) {
  exprs20<-t(aggregate(t(exprs20),by=list(clust2),FUN="median"))[-1,] #internal of this function G.n
  }
    data1 <- t(exprs1)
    data2 <- t(exprs2)
    
    data10 <- t(exprs10)
    data20 <- t(exprs20)
    
    n1 <- nrow(data1)
    n2 <- nrow(data2)
    n.total <- n1 + n2
   
    W21 <- FNN::get.knnx(data2, query=data1, k=k1)
    W12 <- FNN::get.knnx(data1, query=data2, k=k2)
    
    inds12<-cbind( as.vector(rep(seq_len(n1), each=k1)),  as.vector(t(W21$nn.index)) )
    inds21<-cbind(as.vector(t(W12$nn.index)) , as.vector(rep(seq_len(n2), each=k2)) )
    
    A<-rbind(inds12,inds21)
    keeps=duplicated(A)  ##duplicate rows of A are mutaul nns
    A<-A[keeps,]
    
    if (withQC) {
    kn=min(nrow(data1),nrow(data2),100)
    #norm.thr=round(kn-4*sqrt(kn)) #round( (kn+sqrt(kn))/2) #between sum of kn random (orthogonal) vectors and max
    A<-qcMNNs(A,data1,data2,kk=kn)
    }
    # Computing the batch correction vector between MNN pairs.
    A1 <- A[,1]
    A2 <- A[,2] 
    vect <- data1[A1,] - data2[A2,]    
    vect0 <- data10[A1,] - data20[A2,] 
    
    ###matched shift variance in data1
    if (varCare==TRUE){
    tvecte<-cosine.norm(t(vect))
    data2p<-data2 %*% tvecte #n*m (m vects) 
    var2=colVars(data2p) #variance of data on the direction of each vect
    r.var2=runif(nrow(vect),min=0,max=sqrt(var2)/2) #random size elongated shift
    vect<-vect+r.var2 * t(tvecte)
    mnn.seeds<-data1[A1,]-vect   #change mnn seeds in data2 to adjust for variance
    
    samples<-sample(1:5,length(A2),replace=TRUE)
    W2n<- FNN::get.knnx(data2, query=mnn.seeds, k=5)$nn.index
    A2n<-sapply(seq_len(nrow(W2n)), function(i) W2n[i,samples[i]]) # update A2

    #print(sum(A2==A2n)/length(A2))
    A2<-A2n
    vect <- data1[A1,] - data2[A2,]    # update vect
    vect0 <- data10[A1,] - data20[A2,] # update vect0
    ###matched shift variance in data2
    tvecte<-cosine.norm(t(vect))
    data1p<-data1 %*% tvecte #n*m (m vects) 
    var1=colVars(data1p) #variance of data on the direction of each vect
    r.var1=runif(nrow(vect),min=0,max=sqrt(var1)/2) #random size elongated shift
    vect<-vect+r.var1 * t(tvecte)

    tvect0e<-cosine.norm(t(vect0))
    data10p<-data10 %*% tvect0e #n*m (m vects) 
    var10=colVars(data10p) #variance of data on the direction of each vect
    r.var10=runif(nrow(vect0),min=0,max=sqrt(var10)/2)
    vect0<-vect0+r.var10* t(tvect0e)
    }
      
    # Gaussian smoothing of individual correction vectors for MNN pairs.
    if (sigma==0) {
        G <- matrix(1, n2, n2)
    } else if (n2< 3000) {
      dd2 <- as.matrix(dist(data2))  #colsums of G give the density !until else
      G <- exp(-dd2^2/sigma)
      # G <- matrix(0,n2,n2)
      # dd2 <- unclass(proxy::dist(data2[A2,], data2)) 
      # G[A2,seq_len(nrow(data2))]=exp(-dd2^2/sigma) #colsums of G give the density
      # G<- t(G) #rowsums of G give the density
    } else {
        kk <- min(length(A2),100)
        W <- get.knnx(data2[A2,], query=data2, k=kk)
        
        G <- matrix(0,n2,n2)
        # for (i in seq_len(n2)) { 
        #     #G[i,A2[W$nn.index[i,]]]=W$}
        #     G[i,A2[W$nn.index[i,]]]=exp(-(W$nn.dist[i,])^2/sigma) 
        # }
        nonlist3<-t(exp(-(W$nn.dist)^2/sigma) )#lapply(seq_len(nrow(W$nn.dist)), function(i) exp(-(W$nn.dist[i,])^2/sigma)  )
        list2<-lapply(seq_len(nrow(W$nn.dist)), function(i) A2[W$nn.index[i,]] )
        nonlist1<-rep(seq_len(n2),each=kk)
        G[cbind(nonlist1,unlist(list2))]=nonlist3  #rowsums of G give the density 
    }
    #G <- (G+t(G)) /2
    
    #################
    D <- rowSums(G)
    nA2 <- tabulate(A2, nbins=n2)
    norm.dens <- t(G/(D*nA2))[,A2,drop=FALSE] # density normalized to avoid domination from dense parts
    batchvect <- norm.dens %*% vect 
    partitionf <- rowSums(norm.dens)
    partitionf[partitionf==0]<-1  # to avoid nans (instead get 0s)
    batchvect <- batchvect/partitionf

    batchvect0 <- norm.dens %*% vect0 
    batchvect0 <- batchvect0/partitionf
    if(!is.null(clust2)) {
    batchvect00<- batchvect0[clust2,]
    }else{batchvect00<- batchvect0}    
    # Report cells that are MNNs, and the correction vector per cell in data2.
    list(set1=unique(A1), set2=unique(A2), vect=batchvect, vect0=batchvect00)
}

get.bio.span <- function(exprs, inquiry.in.hvg, ndim) 
# Computes the basis matrix of the biological subspace of 'exprs'.
# The first 'ndim' dimensions are assumed to capture the biological subspace.
# Avoids extra dimensions dominated by technical noise, which will result in both 
# trivially large and small angles when using find.shared.subspace().
{
    keeph <- numeric(nrow(exprs))
    keeph[inquiry.in.hvg] <- 1
    exprs <- exprs * keeph
    exprs <- exprs - rowMeans(exprs) 
    S <- svd(exprs, nu=ndim, nv=0)
    S$u   
}

####
cosine.norm <- function(X)
# Computes the cosine norm, with some protection from zero-length norms.
#G*N  
{
    cellnorm <- pmax(1e-8, sqrt(colSums(X^2)))
    X/matrix(cellnorm, nrow(X), ncol(X), byrow=TRUE)
}

####
qcMNNs<-function (mnns,data1,data2,kk=100) {
  #data n*g
  #kk is no. nearest neighbours to calculate sum of vectors
  #library(FNN)
  
  ###set the norm.thr to 80% quantile of all cell's vects.norm in data1
#  if (is.null(norm.thr)) {
  data1 <- data1 - rowMeans(data1) 
  S <- svd(data1, nu=2, nv=0)
  data1<-S$u 
  
  data2 <- data2 - rowMeans(data2) 
  S <- svd(data2, nu=2, nv=0)
  data2<-S$u 
  
  W <- FNN::get.knnx(data1, query=data1, k=kk)
  id.list<-lapply(seq_len(nrow(W$nn.index)), function(i) W$nn.index[i,])
  center.list<-lapply(seq_len(nrow(W$nn.index)), function(i) matrix( rep(data1[i,],each=kk),nrow=kk,ncol=ncol(data1) ) )
  
  vects<-(lapply(seq_len(nrow(W$nn.index)), function(i) data1[W$nn.index[i,],] - center.list[[i]] )) #connecting vectors to each neighbour
  vects.n<-lapply((vects), function(v) (colSums(t(cosine.norm(t(v))))) ) #sum of normalized vectors 
  vects.norm<-sapply((vects.n), function(v) sqrt(sum(v^2)) ) #norm of the sum vector 
  
    #norm.thr=round( (kk-4*sqrt(kk)) )
    norm.thr=quantile( vects.norm ,0.80)
 #   }
  ##check data1
  MNNS=mnns[,1]
  dat1<-data1[MNNS,]
  
  W <- FNN::get.knnx(data1, query=dat1, k=kk)
  id.list<-lapply(seq_len(nrow(W$nn.index)), function(i) W$nn.index[i,])
  center.list<-lapply(seq_len(nrow(W$nn.index)), function(i) matrix( rep(data1[MNNS[i],],each=kk),nrow=kk,ncol=ncol(data1) ) )
  
  vects<-(lapply(seq_len(nrow(W$nn.index)), function(i) data1[W$nn.index[i,],] - center.list[[i]] )) #connecting vectors to each neighbour
  vects.n<-lapply((vects), function(v) (colSums(t(cosine.norm(t(v))))) ) #sum of normalized vectors 
  vects.norm<-sapply((vects.n), function(v) sqrt(sum(v^2)) ) #norm of the sum vector 
  
  goodMNN=matrix(0,nrow=nrow(mnns),ncol=2)
  goodMNN[,1]<-(vects.norm<norm.thr)
  
  hist(vects.norm,100, main=paste0("norm (upper)threshold for sum of 100 vectors= ",norm.thr) )
  
  ratio=sum(goodMNN[,1])/nrow(goodMNN)
  if (ratio==0) {
    warning("All mnns in the reference batch failed the QC edge test. results are shown without QC.Try running the function with a lower norm.thr for withQC.")
    #keeps=seq_len(nrow(mnns))
    goodMNN[,1]=seq_len(nrow(mnns))
    } else {print(paste0("Accepted ratio of mnns in the reference batch is ",round(ratio*100), "%.") ) 
      #keeps=rowSums(goodMNN)>0 # if one of the pairs is a good (not on the edge) the pair is accepted
    }
  
  ###set the norm.thr to 80% quantile of all cell's vects.norm in data2
  #if (is.null(norm.thr)) {
    
    W <- FNN::get.knnx(data2, query=data2, k=kk)
    id.list<-lapply(seq_len(nrow(W$nn.index)), function(i) W$nn.index[i,])
    center.list<-lapply(seq_len(nrow(W$nn.index)), function(i) matrix( rep(data2[i,],each=kk),nrow=kk,ncol=ncol(data2) ) )
    
    vects<-(lapply(seq_len(nrow(W$nn.index)), function(i) data2[W$nn.index[i,],] - center.list[[i]] )) #connecting vectors to each neighbour
    vects.n<-lapply((vects), function(v) (colSums(t(cosine.norm(t(v))))) ) #sum of normalized vectors 
    vects.norm<-sapply((vects.n), function(v) sqrt(sum(v^2)) ) #norm of the sum vector 
    
    #norm.thr=round( (kk-4*sqrt(kk)) )
    norm.thr=quantile( vects.norm ,0.80)
 # }
  ##check data2
  MNNS=mnns[,2]
  dat2<-data2[MNNS,]
  
  W <- FNN::get.knnx(data2, query=dat2, k=kk)
  id.list<-lapply(seq_len(nrow(W$nn.index)), function(i) W$nn.index[i,])
  center.list<-lapply(seq_len(nrow(W$nn.index)), function(i) matrix( rep(data2[MNNS[i],],each=kk),nrow=kk,ncol=ncol(data2) ) )
  
  vects<-(lapply(seq_len(nrow(W$nn.index)), function(i) data2[W$nn.index[i,],]- center.list[[i]] ))
  vects.n<-lapply((vects), function(v) (colSums(t(cosine.norm(t(v))))) ) #goes wrong
  vects.norm<-sapply((vects.n), function(v) sqrt(sum(v^2)) )
  
  goodMNN[,2]<-(vects.norm<norm.thr)
  
  hist(vects.norm,100, main=paste0("norm (upper)threshold for sum of 100 vectors= ",norm.thr) )
  
  ratio=sum(goodMNN[,2])/nrow(goodMNN)
  if (ratio==0) {
    warning("All mnns in the second batch failed the QC edge test. results are shown without QC.Try running the function with a lower norm.thr for withQC.")
    #keeps=seq_len(nrow(mnns))  
    goodMNN[,2]=seq_len(nrow(mnns))
    } else {print(paste0("Accepted ratio of mnns in the second batch is ",round(ratio*100), "%.") ) 
      #keeps=rowSums(goodMNN)>0 # if one of the pairs is a good (not on the edge) the pair is accepted
    }
  keeps=rowSums(goodMNN)>0 #if one of the pairs is a good (not on the edge) the pair is accepted
  #keeps=rowSums(goodMNN)>1 #only if both are not edge
  return(mnns[keeps,])
  
}
