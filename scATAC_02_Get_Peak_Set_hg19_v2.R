#Creating a peak set, summarized experiment and LSI clustering
#07/31/19
#Adapted from Satpathy*, Granja*, et al. 
#Massively parallel single-cell chromatin landscapes of human immune 
#cell development and intratumoral T cell exhaustion (2019)
#Created by Jeffrey Granja
library(Matrix)
library(SummarizedExperiment)
library(matrixStats)
library(readr)
library(GenomicRanges)
library(magrittr)
library(edgeR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg19)
set.seed(1)

countInsertions <- function(query, fragments, by = "RG"){
  #Count By Fragments Insertions
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  by <- "RG"
  overlapDF <- DataFrame(findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))
  #Calculate Overlap Stats
  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total
  #Summarize
  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1], 
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)), 
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  return(out)
}

seuratLSI <- function(mat, nComponents = 50, binarize = TRUE, nFeatures = NULL){
  #TF IDF LSI adapted from flyATAC
  cs <- Matrix::colSums(mat)
  if(binarize){
    message(paste0("Binarizing matrix..."))
    mat@x[mat@x > 0] <- 1 
  }
  if(!is.null(nFeatures)){
    message(paste0("Getting top ", nFeatures, " features..."))
    mat <- mat[head(order(Matrix::rowSums(mat),decreasing = TRUE),nFeatures),] 
  }
  #Calc TF IDF
  message("Computing Term Frequency IDF...")
  freqs <- t(t(mat)/Matrix::colSums(mat))
  idf   <- as(log(1 + ncol(mat) / Matrix::rowSums(mat)), "sparseVector")
  tfidf <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% freqs
  #Calc SVD then LSI
  message("Computing SVD using irlba...")
  svd <- irlba::irlba(tfidf, nComponents, nComponents)
  svdDiag <- matrix(0, nrow=nComponents, ncol=nComponents)
  diag(svdDiag) <- svd$d
  matSVD <- t(svdDiag %*% t(svd$v))
  rownames(matSVD) <- colnames(mat)
  colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))
  #Make Seurat Object
  message("Making Seurat Object...")
  mat <- mat[1:100,] + 1
  obj <- CreateSeuratObject(mat, project='scATAC', min.cells=0, min.genes=0)
  
  # # new setDimReduction code
  # dimObj <- CreateDimReducObject(
  #   embeddings = matSVD,
  #   key = "PC"
  # )
  # obj[["pca"]] <- dimObj
  
  # outdated code for Seurat 
  obj <- SetDimReduction(object = obj, reduction.type = "pca", slot = "cell.embeddings", new.data = matSVD)
  obj <- SetDimReduction(object = obj, reduction.type = "pca", slot = "key", new.data = "PC")
  return(obj)
}

addClusters <- function(obj, minGroupSize = 50, dims.use = seq_len(50), initialResolution = 0.8){
  #First Iteration of Find Clusters
  currentResolution <- initialResolution
  
  # below is outdated
  obj <- FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, resolution = currentResolution, print.output = FALSE)
  # obj <- FindNeighbors(object = obj, reduction = "pca", dims = dims.use) # new code
  # obj <- FindClusters(object = obj, resolution = currentResolution) # new code
  
  
  minSize <- min(table(obj@meta.data[[paste0("res.",currentResolution)]]))
  # minSize <- min(table(obj@meta.data[[paste0("RNA_snn_res.",currentResolution)]]))
  
  nClust <- length(unique(paste0(obj@meta.data[[paste0("res.",currentResolution)]])))
  # nClust <- length(unique(paste0(obj@meta.data[[paste0("RNA_snn_res.",currentResolution)]])))
  
  message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s", currentResolution, nClust, minSize))
  #If clusters are smaller than minimum group size
  while(minSize <= minGroupSize){
    obj@meta.data <- obj@meta.data[,-which(colnames(obj@meta.data)==paste0("res.",currentResolution))]
    # obj@meta.data <- obj@meta.data[,-which(colnames(obj@meta.data)==paste0("RNA_snn_res.",currentResolution))]
    
    currentResolution <- currentResolution*initialResolution
    
    # obj <- FindNeighbors(object = obj, reduction = "pca", dims = dims.use) # new code
    # obj <- FindClusters(object = obj, resolution = currentResolution) # new code
    
    obj <- FindClusters(object = obj, reduction.type = "pca", dims.use = dims.use, resolution = currentResolution, print.output = FALSE, force.recalc = TRUE)
    
    minSize <- min(table(obj@meta.data[[paste0("res.",currentResolution)]]))
    nClust <- length(unique(paste0(obj@meta.data[[paste0("res.",currentResolution)]])))
    
    # minSize <- min(table(obj@meta.data[[paste0("RNA_snn_res.",currentResolution)]]))
    # nClust <- length(unique(paste0(obj@meta.data[[paste0("RNA_snn_res.",currentResolution)]])))
    
    message(sprintf("Current Resolution = %s, No of Clusters = %s, Minimum Cluster Size = %s", currentResolution, nClust, minSize))
  }
  return(obj)
}

extendedPeakSet <- function(df, BSgenome = NULL, extend = 250, blacklist = NULL, nSummits = 100000){
  #Helper Functions
  readSummits <- function(file){
    df <- suppressMessages(data.frame(readr::read_tsv(file, col_names = c("chr","start","end","name","score"))))
    df <- df[,c(1,2,3,5)] #do not keep name column it can make the size really large
    return(GenomicRanges::makeGRangesFromDataFrame(df=df,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE))
  }
  nonOverlappingGRanges <- function(gr, by = "score", decreasing = TRUE, verbose = FALSE){
    stopifnot(by %in% colnames(mcols(gr)))
    clusterGRanges <- function(gr, filter = TRUE, by = "score", decreasing = TRUE){
      gr <- sort(sortSeqlevels(gr))
      r <- GenomicRanges::reduce(gr, min.gapwidth=0L, ignore.strand=TRUE)
      o <- findOverlaps(gr,r)
      mcols(gr)$cluster <- subjectHits(o)
      gr <- gr[order(mcols(gr)[,by], decreasing = decreasing),]
      gr <- gr[!duplicated(mcols(gr)$cluster),]
      gr <- sort(sortSeqlevels(gr))
      mcols(gr)$cluster <- NULL
      return(gr)
    }
    if(verbose){
      message("Converging", appendLF = FALSE)
    }
    i <-  0
    gr_converge <- gr
    while(length(gr_converge) > 0){
      if(verbose){
        message(".", appendLF = FALSE)
      }
      i <-  i + 1
      gr_selected <- clusterGRanges(gr = gr_converge, filter = TRUE, by = by, decreasing = decreasing)
      gr_converge <- subsetByOverlaps(gr_converge ,gr_selected, invert=TRUE) #blacklist selected gr
      if(i == 1){ #if i=1 then set gr_all to clustered
        gr_all <- gr_selected
      }else{
        gr_all <- c(gr_all, gr_selected)
      } 
    }
    if(verbose){
      message("\nSelected ", length(gr_all), " from ", length(gr))
    }
    gr_all <- sort(sortSeqlevels(gr_all))
    return(gr_all)
  }
  #Check-------
  stopifnot(extend > 0)
  stopifnot("samples" %in% colnames(df))
  stopifnot("groups" %in% colnames(df))
  stopifnot("summits" %in% colnames(df))
  stopifnot(!is.null(BSgenome))
  stopifnot(all(apply(df,1,function(x){file.exists(paste0(x[3]))})))
  #------------
  #Deal with blacklist
  if(is.null(blacklist)){
    blacklist <- GRanges()
  }else if(is.character(blacklist)){
    blacklist <- rtracklayer::import.bed(blacklist)
  }
  stopifnot(inherits(blacklist,"GenomicRanges"))
  #------------
  #Time to do stuff
  chromSizes <- GRanges(names(seqlengths(BSgenome)), IRanges(1, seqlengths(BSgenome)))
  chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
  groups <- unique(df$groups)
  groupGRList <- GenomicRanges::GenomicRangesList(lapply(seq_along(groups), function(i){
      df_group = df[which(df$groups==groups[i]),]
      grList <- GenomicRanges::GenomicRangesList(lapply(paste0(df_group$summits), function(x){
        extended_summits <- readSummits(x) %>%
          resize(., width = 2 * extend + 1, fix = "center") %>%     
          subsetByOverlaps(.,chromSizes,type="within") %>%
          subsetByOverlaps(.,blacklist,invert=TRUE) %>%
          nonOverlappingGRanges(., by="score", decreasing=TRUE)
        extended_summits <- extended_summits[order(extended_summits$score,decreasing=TRUE)]
        if(!is.null(nSummits)){
          extended_summits <- head(extended_summits, nSummits)
        }
        mcols(extended_summits)$scoreQuantile <- trunc(rank(mcols(extended_summits)$score))/length(mcols(extended_summits)$score)
        extended_summits
      }))
      #Non Overlapping
      grNonOverlapping <- nonOverlappingGRanges(unlist(grList), by = "scoreQuantile", decreasing = TRUE)
      #Free Up Memory
      remove(grList)
      gc()
      grNonOverlapping
    }))
  grFinal <- nonOverlappingGRanges(unlist(groupGRList), by = "scoreQuantile", decreasing = TRUE)
  grFinal <- sort(sortSeqlevels(grFinal))
  return(grFinal)
}

groupSums <- function(mat, groups = NULL, na.rm = TRUE, sparse = FALSE){
    stopifnot(!is.null(groups))
    stopifnot(length(groups) == ncol(mat))
    gm <- lapply(unique(groups), function(x) {
        if (sparse) {
            Matrix::rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }else {
            rowSums(mat[, which(groups == x), drop = F], na.rm = na.rm)
        }
    }) %>% Reduce("cbind", .)
    colnames(gm) <- unique(groups)
    return(gm)
}

#-------------------------------------------------------------------------------------------------
# Start
#-------------------------------------------------------------------------------------------------
# this is using the fragments from the filtered cells in rds file previously generated from script 01
# will read all, healthy disease and all rds files, selecte carefully below
fragmentFiles <- list.files("/projectnb/paxlab/isarfraz/Data", pattern = ".rds", full.names = TRUE)
# if using next time, it will read all RDS files, so filter to include only fragments.rds files from script 01
# fragmentFiles <- fragmentFiles[1]

#-------------------------------------------------------------------------------------------------
# Get Counts In Windows
#-------------------------------------------------------------------------------------------------
# REF #
# ref hg19 genome
genome <- BSgenome.Hsapiens.UCSC.hg19
# making granges object from ref genome (each row is one chr)
chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
# possibly fitering some chr in ref
chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
# making ranges of each chr in ref 
windows <- unlist(tile(chromSizes, width = 2500)) 

# NOW OUR DATA #
# from fragments file, it creates a counts matrix
# columns are cells
# rows are chr windows (ranges)
countsList <- lapply(seq_along(fragmentFiles), function(i){
	message(sprintf("%s of %s", i, length(fragmentFiles)))
	counts <- countInsertions(windows, readRDS(fragmentFiles[i]), by = "RG")[[1]]
  counts
})
mat <- lapply(countsList, function(x) x) %>% Reduce("cbind",.)
remove(countsList)
gc()
# so using ref genome, it creates counts of each chr window/region from our data

#-------------------------------------------------------------------------------------------------
# Run LSI Clustering with Seurat
#-------------------------------------------------------------------------------------------------
set.seed(1)
message("Making Seurat LSI Object...")
# LSI is latent semantic indexing (dimred)
obj <- seuratLSI(mat, nComponents = 25, nFeatures = 20000) 
# possibly groups 20000 rows into one (windows), using LSI? maybe this is dimred
# output is 100 features x full cells

# # setting rownames to mat
# rownames(mat) <- paste0(rep("f", nrow(mat)), rep(1:nrow(mat), each = 1))
# # run LSI now
# obj <- seuratLSI(mat, nComponents = 25, nFeatures = 20000)

message("Adding Graph Clusters...")
obj <- addClusters(obj, dims.use = 2:25, minGroupSize = 200, initialResolution = 0.8)
# starts clustering from initial resolution of 0.8 and keeps going until every cluster >= 200 cells

saveRDS(obj, "/projectnb/paxlab/isarfraz/Data/Save-LSI-Windows-Seurat.rds")
clusterResults <- split(rownames(obj@meta.data), paste0("Cluster",obj@meta.data[,ncol(obj@meta.data)]))
remove(obj)
gc()

#-------------------------------------------------------------------------------------------------
# Get Cluster Beds
#-------------------------------------------------------------------------------------------------

# for each cluster it is now separately cells from that cluster into each clusters own bed file (like a granges df)

dirClusters <- "/projectnb/paxlab/isarfraz/Data/LSI-Cluster-Beds/"
dir.create(dirClusters)
for(i in seq_along(fragmentFiles)){
	fragments <-readRDS(fragmentFiles[i])
	for(j in seq_along(clusterResults)){
	  message(sprintf("%s of %s", j, length(clusterResults)))
	  fragmentsj <- fragments[fragments$RG %in% clusterResults[[j]]]
	  if(length(fragmentsj) > 0){
	    out <- data.frame(
	      chr = c(seqnames(fragmentsj), seqnames(fragmentsj)), 
	      start = c(as.integer(start(fragmentsj) - 1), as.integer(end(fragmentsj) - 1)), 
	      end = c(as.integer(start(fragmentsj)), as.integer(end(fragmentsj)))
	      ) %>% readr::write_tsv(
          x = ., 
          append = TRUE, 
          path = paste0(dirClusters, paste0(names(clusterResults)[j], ".bed")), 
          col_names = FALSE)
	    }
	}
}

#-------------------------------------------------------------------------------------------------
# Run MACS2
#-------------------------------------------------------------------------------------------------

# uses previous bed files for each cluster
# then for each cluster it, macs was used for peak calling (identify accessible regions)


dirPeaks <- "/projectnb/paxlab/isarfraz/Data/LSI-Cluster-Peaks/"
method <- "q"
cutoff <- 0.05
shift <- -75
extsize <- 150
genome_size <- 2.7e9
for(j in seq_along(clusterResults)){
	message(sprintf("%s of %s", j, length(clusterResults)))
	clusterBedj <- paste0(dirClusters,names(clusterResults)[j],".bed")
	cmdPeaks <- sprintf(
	    "macs2 callpeak -g %s --name %s --treatment %s --outdir %s --format BED --nomodel --call-summits --nolambda --keep-dup all", 
	    genome_size, 
	    names(clusterResults)[j], 
	    clusterBedj, 
	    dirPeaks
	  )
	if (!is.null(shift) & !is.null(extsize)) {
	  cmdPeaks <- sprintf("%s --shift %s --extsize %s", cmdPeaks, shift, extsize)
	}
	if (tolower(method) == "p") {
	  cmdPeaks <- sprintf("%s -p %s", cmdPeaks, cutoff)
	}else {
	  cmdPeaks <- sprintf("%s -q %s", cmdPeaks, cutoff)
	}
	message("Running Macs2...")
	message(cmdPeaks)
	system(cmdPeaks, intern = TRUE)
}

#-------------------------------------------------------------------------------------------------
# Make Non-Overlapping Peak Set
#-------------------------------------------------------------------------------------------------

# overlapping peaks are merged

dirPeaks <- "/projectnb/paxlab/isarfraz/Data/LSI-Cluster-Peaks" # added by me because of /
df <- data.frame(
  samples = gsub("\\_summits.bed","",list.files(dirPeaks, pattern = "\\_summits.bed", full.names = FALSE)),
  groups = "scATAC",
  summits = list.files(dirPeaks, pattern = "\\_summits.bed", full.names = TRUE)
  )

# downloaded hg19.blacklist.bed from https://github.com/Boyle-Lab/Blacklist/blob/master/lists/Blacklist_v1/hg19-blacklist.bed.gz
# problematic regions of genome (hat have anomalous, unstructured, or high signal) - removal required
# extendedPeakSet reads peaks previously identified and possibly removes blacklisted ones
unionPeaks <- extendedPeakSet(
    df = df,
    BSgenome = genome, 
    extend = 250,
    blacklist = "/projectnb/paxlab/isarfraz/Data/hg19-blacklist.bed",
    nSummits = 200000
  )
unionPeaks <- unionPeaks[seqnames(unionPeaks) %in% paste0("chr",c(1:22,"X"))]
unionPeaks <- keepSeqlevels(unionPeaks, paste0("chr",c(1:22,"X")))

#Create Counts list from peaks 
countsPeaksList <- lapply(seq_along(fragmentFiles), function(i){
  message(sprintf("%s of %s", i, length(fragmentFiles)))
  gc()
  countInsertions(unionPeaks, readRDS(fragmentFiles[i]), by = "RG")
})

#CountsMatrix from peaks
mat <- lapply(countsPeaksList, function(x) x[[1]]) %>% Reduce("cbind",.)
frip <- lapply(countsPeaksList, function(x) x[[2]]) %>% unlist
total <- lapply(countsPeaksList, function(x) x[[3]]) %>% unlist

dim(mat) # filtered peaks = (identified by macs and filtered from blacklisted)

# object for downstream analysis
se <- SummarizedExperiment(
  assays = SimpleList(counts = mat), 
  rowRanges = unionPeaks
  )
rownames(se) <- paste(seqnames(se),start(se),end(se),sep="_")
colData(se)$FRIP <- frip
colData(se)$uniqueFrags <- total / 2
saveRDS(se, "/projectnb/paxlab/isarfraz/Data/scATAC-Summarized-Experiment.rds")


