library(Matrix)
library(SummarizedExperiment)
library(tidyverse)
library(uwot)
library(edgeR)
library(FNN)
library(matrixStats)
library(Rcpp)
set.seed(1)


#Optimized LSI for scRNA-seq analysis
optimizeLSI <- function(inSCE, mat, scaleTo = 10000, priorCount = 3, pcsUse = 1:25, 
                        resolution = c(0.2, 0.4, 0.8), varFeatures = c(2500, 2500, 2500), seed = 1){
  
  set.seed(seed)
  stopifnot(length(resolution) > 1)
  stopifnot(length(resolution) == length(varFeatures))
  
  #Initialize List
  lsiOut <- list()
  
  #Initial LSI uses variances that are across all single cells and will have larger batch relationships
  i <- 1
  message("Initial LSI...")
  matNorm <- t(t(mat)/Matrix::colSums(mat)) * scaleTo
  matNorm@x <- log2(matNorm@x + 1)
  idVarFeatures <- head(order(sparseRowVariances(matNorm),decreasing=TRUE), varFeatures[i])
  lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
  # clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)
  pca <- lsiObj$matSVD
  rownames(pca) <- gsub("_", "-", rownames(pca))
  new_pca <- Seurat::CreateDimReducObject(embeddings = pca, key = "PC_")
  inSCE <- runSeuratFindClusters(inSCE, externalReduction = new_pca, dims = pcsUse, resolution = resolution[i])
  clusters <- colData(inSCE)[[paste0("Seurat_louvain_Resolution", resolution[i])]]
  
  #Store
  lsiOut[[paste0("iter", i)]] <- list(
    lsiMat = lsiObj$matSVD, 
    varFeatures = idVarFeatures, 
    clusters = clusters
  )
  
  for(i in seq(2, length(varFeatures))){
    
    message(sprintf("Additional LSI %s...", i))
    
    #Run LSI
    clusterMat <- edgeR::cpm(groupSums(mat, clusters, sparse = TRUE), log=TRUE, prior.count = priorCount)
    idVarFeatures <- head(order(rowVars(clusterMat), decreasing=TRUE), varFeatures[i])
    lsiObj <- calcLSI(mat[idVarFeatures,], binarize = FALSE, nComponents = max(pcsUse))
    # clusters <- seuratSNN(lsiObj$matSVD, dims.use = pcsUse, resolution = resolution[i], n.start = 10, print.output = FALSE)
    pca <- lsiObj$matSVD
    rownames(pca) <- gsub("_", "-", rownames(pca))
    new_pca <- Seurat::CreateDimReducObject(embeddings = pca, key = "PC_")
    inSCE <- runSeuratFindClusters(inSCE, externalReduction = new_pca, dims = pcsUse, resolution = resolution[i])
    clusters <- colData(inSCE)[[paste0("Seurat_louvain_Resolution", resolution[i])]]
    
    
    if(i == length(varFeatures)){
      #Save All Information from LSI Attempt
      lsiOut[[paste0("iter", i)]] <- list(
        lsiObj = lsiObj, 
        varFeatures = idVarFeatures, 
        clusters = clusters,
        matNorm = matNorm
      )
    }else{
      lsiOut[[paste0("iter", i)]] <- list(
        lsiMat = lsiObj$matSVD, 
        varFeatures = idVarFeatures, 
        clusters = clusters
      )
    }
    
  }
  
  return(lsiOut)
  
}

#Compute Fast Sparse Row Variances
sparseRowVariances <- function (m){
  rM <- Matrix::rowMeans(m)
  rV <- computeSparseRowVariances(m@i + 1, m@x, rM, ncol(m))
  return(rV)
}

#Sparse Variances Rcpp
sourceCpp(code='
  #include <Rcpp.h>
  using namespace Rcpp;
  using namespace std;
  // [[Rcpp::export]]
  Rcpp::NumericVector computeSparseRowVariances(IntegerVector j, NumericVector val, NumericVector rm, int n) {
    const int nv = j.size();
    const int nm = rm.size();
    Rcpp::NumericVector rv(nm);
    Rcpp::NumericVector rit(nm);
    int current;
    // Calculate RowVars Initial
    for (int i = 0; i < nv; ++i) {
      current = j(i) - 1;
      rv(current) = rv(current) + (val(i) - rm(current)) * (val(i) - rm(current));
      rit(current) = rit(current) + 1;
    }
    // Calculate Remainder Variance
    for (int i = 0; i < nm; ++i) {
      rv(i) = rv(i) + (n - rit(i))*rm(i)*rm(i);
    }
    rv = rv / (n - 1);
    return(rv);
  }'
)

calcLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
  
  set.seed(1)
  
  #TF IDF LSI adapted from flyATAC
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    idx <- head(order(Matrix::rowSums(mat), decreasing = TRUE), nFeatures)
    mat <- mat[idx,] 
  }else{
    idx <- which(Matrix::rowSums(mat) > 0)
    mat <- mat[idx,]
  }
  
  #Calc RowSums and ColSums
  colSm <- Matrix::colSums(mat)
  rowSm <- Matrix::rowSums(mat)
  
  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/colSm)
  idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  
  #Return Object
  out <- list(
    matSVD = matSVD, 
    rowSm = rowSm, 
    colSm = colSm, 
    idx = idx, 
    svd = svd, 
    binarize = binarize, 
    nComponents = nComponents,
    date = Sys.Date(),
    seed = 1)
  
  out
  
}

#Helper function for summing sparse matrix groups
groupSums <- function (mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
  stopifnot(!is.null(groups))
  stopifnot(length(groups) == ncol(mat))
  gm <- lapply(unique(groups), function(x) {
    if (sparse) {
      Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
    else {
      rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
    }
  }) %>% Reduce("cbind", .)
  colnames(gm) <- unique(groups)
  return(gm)
}
