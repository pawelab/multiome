library(SingleCellExperiment)
library(singleCellTK)

#TODO:
# see if they do any filtering with their rds file
# see if their rds has same dims
# run lsi too (is it run on counts matrix) https://quanteda.io/articles/pkgdown/examples/lsa.html
# where are they getting clusters from?
# they did remove ribosomal and mito genes (see if already removed, otherwise remove)

# Input Data 
################
rna_exp <- readRDS("/projectnb/paxlab/isarfraz/Data/GSM4138876_scRNA_PBMC_D4T1.rds")
# one healthy rna exp sample 
sce_hl <- SingleCellExperiment(list(counts = rna_exp))

# Pre-processing (they don't do)
###############################
sce_hl <- runSeuratNormalizeData(sce_hl)
sce_hl <- runSeuratFindHVG(sce_hl)

# Run DimRed (they use LSI)
################
sce_hl <- runSeuratPCA(sce_hl, useFeatureSubset = "hvf")
plotDimRed(sce_hl, useReduction = "seuratPCA") #looks pretty representative of pbmc


# Need Clusters (more res increases clusters)
#################
sce_hl <- runSeuratFindClusters(sce_hl, useReduction = "pca", resolution = 0.8)
sce_hl <- runSeuratFindClusters(sce_hl, useReduction = "pca", resolution = 0.5)
unique(colData(sce_hl)[["Seurat_louvain_Resolution0.8"]])
unique(colData(sce_hl)[["Seurat_louvain_Resolution0.5"]])
# same 12/13 at both res 0.8 and 0.5, probably!


# Run UMAP & visualize
####################
sce_hl <- runSeuratUMAP(sce_hl, useReduction = "pca")
plotSeuratReduction(sce_hl, useReduction = "umap", groupBy = "Seurat_louvain_Resolution0.8", showLegend = TRUE)

# just to find celltypes (just for fun)
sce_hl <- runSingleR(inSCE = sce_hl, useAssay = "seuratNormData", useBltinRef = "hpca")
unique(colData(sce_hl)[["SingleR_hpca_main_labels"]])
plotSeuratReduction(sce_hl, useReduction = "umap", groupBy = "SingleR_hpca_main_labels", showLegend = TRUE)
# the above look pretty representative compared to seurat tutorial too, but these are not divided into cell subtypes
# https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

#####################
## second script
######################

# load healthy cells (use same as above) 
# but loading again because cannot cbind hl and mpal with stored computations
sce_hl <- SingleCellExperiment(list(counts = rna_exp))

# load disease cells
mpal_rna_exp <- readRDS("/projectnb/paxlab/isarfraz/Data/GSM4138878_scRNA_MPAL1_T1.rds")
sce_mpal <- SingleCellExperiment(list(counts = mpal_rna_exp))
# very less cells

# find common features (and remove MT features) (yes there are MT, also remove MRP ribosomal)
sce_hl <- sce_hl[!grepl("^MT-", rownames(sce_hl)), ]
sce_mpal <- sce_mpal[!grepl("^MT-", rownames(sce_mpal)), ]

sce_hl <- sce_hl[!grepl("^MRP", rownames(sce_hl)), ]
sce_mpal <- sce_mpal[!grepl("^MRP", rownames(sce_mpal)), ]

# using common features cbind healthy and disease cells into one object/matrix
common_features <- intersect(rownames(sce_hl), rownames(sce_mpal)) #only few that are not
sce_hl <- sce_hl[common_features, ]
sce_mpal <- sce_mpal[common_features, ]
sce_all <- cbind(sce_hl, sce_mpal) #combined hl and mpal object


# compute dimred (lsi or pca?) from combined matrix (not sure if in lsi func they process them separately of each other)
## lets process combined first
sce_all <- runSeuratNormalizeData(sce_all)
sce_all <- runSeuratFindHVG(sce_all) #would with sce_all it make sense to find hvg and use with pca
sce_all <- runSeuratPCA(sce_all, useFeatureSubset = "hvf", nPCs = 50)
plotDimRed(sce_all, "seuratPCA") #still looks kinda representative of pbmc with slight changes here there

# run umap of combined matrix
sce_all <- runSeuratUMAP(sce_all, useReduction = "pca", dims = 25, minDist = 0.5, nNeighbors = 30)
plotSeuratReduction(sce_all, useReduction = "umap")


# try LSI
#########
lsiObj <- optimizeLSI(inSCE = sce_all, mat = assay(sce_all, "counts"))
sceLSI <- sce_all

pca <- lsiObj[[length(lsiObj)]]$lsiObj$matSVD[,1:25]
rownames(pca) <- gsub("_", "-", rownames(pca))
new_pca <- Seurat::CreateDimReducObject(embeddings = pca, key = "PC_")

sceLSI <- runSeuratUMAP(sceLSI, externalReduction = new_pca, dims = 25, minDist = 0.5, nNeighbors = 30)
sceLSI <- runSeuratFindClusters(sceLSI, externalReduction = new_pca)

plotSCEDimReduceColData(inSCE = sceLSI, colorBy = "Seurat_louvain_Resolution0.8", reducedDimName = "seuratUMAP")
# not same but they are using all samples, me only 1 healthy and 1 mpal

# classify cells







