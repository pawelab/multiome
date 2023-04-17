# module load macs2
# module load R
# R

setwd("/projectnb/paxlab/isarfraz/RProjects/Non-Coding-Mutations-Multiome-Analysis")
.libPaths("/usr2/collab/isarfraz/R/x86_64-pc-linux-gnu-library/4.2")
library(ArchR)
addArchRThreads(threads = 16) 
addArchRGenome("hg19")
inputFiles <- getTutorialData("Hematopoiesis")
ArrowFiles <- createArrowFiles(inputFiles = inputFiles, sampleNames = names(inputFiles),filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE)
doubScores <- addDoubletScores(input = ArrowFiles, k = 10, knnMethod = "UMAP", LSIMethod = 1)
projHeme1 <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = "HemeTutorial", copyArrows = FALSE)
projHeme2 <- filterDoublets(projHeme1)
projHeme2 <- addIterativeLSI(ArchRProj = projHeme2, useMatrix = "TileMatrix", name = "IterativeLSI", iterations = 2, clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10), varFeatures = 25000, dimsToUse = 1:30)
projHeme2 <- addHarmony(ArchRProj = projHeme2,reducedDims = "IterativeLSI",name = "Harmony", groupBy = "Sample")
projHeme2 <- addClusters(input = projHeme2,reducedDims = "IterativeLSI",method = "Seurat",name = "Clusters",resolution = 0.8)
projHeme2 <- addUMAP(ArchRProj = projHeme2, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine")
projHeme2 <- addTSNE(ArchRProj = projHeme2, reducedDims = "IterativeLSI", name = "TSNE", perplexity = 30)
projHeme2 <- addUMAP(ArchRProj = projHeme2, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, metric = "cosine")
projHeme2 <- addTSNE(ArchRProj = projHeme2, reducedDims = "Harmony", name = "TSNEHarmony", perplexity = 30)
projHeme2 <- addImputeWeights(projHeme2)
if(!file.exists("scRNA-Hematopoiesis-Granja-2019.rds")){
  download.file(
    url = "https://jeffgranja.s3.amazonaws.com/ArchR/TestData/scRNA-Hematopoiesis-Granja-2019.rds",
    destfile = "scRNA-Hematopoiesis-Granja-2019.rds"
  )
}

seRNA <- readRDS("scRNA-Hematopoiesis-Granja-2019.rds")
projHeme2 <- addGeneIntegrationMatrix(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix",matrixName = "GeneIntegrationMatrix",reducedDims = "IterativeLSI",seRNA = seRNA,addToArrow = FALSE,groupRNA = "BioClassification",nameCell = "predictedCell_Un",nameGroup = "predictedGroup_Un",nameScore = "predictedScore_Un")
projHeme3 <- addGeneIntegrationMatrix(ArchRProj = projHeme2, useMatrix = "GeneScoreMatrix",matrixName = "GeneIntegrationMatrix",reducedDims = "IterativeLSI",seRNA = seRNA,addToArrow = TRUE,force= TRUE,groupRNA = "BioClassification",nameCell = "predictedCell",nameGroup = "predictedGroup",nameScore = "predictedScore")
projHeme3 <- addImputeWeights(projHeme3)
cM <- confusionMatrix(projHeme3$Clusters, projHeme3$predictedGroup)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
remapClust <- c(
  "01_HSC" = "Progenitor",
  "02_Early.Eryth" = "Erythroid",
  "03_Late.Eryth" = "Erythroid",
  "04_Early.Baso" = "Basophil",
  "05_CMP.LMPP" = "Progenitor",
  "06_CLP.1" = "CLP",
  "07_GMP" = "GMP",
  "08_GMP.Neut" = "GMP",
  "09_pDC" = "pDC",
  "10_cDC" = "cDC",
  "11_CD14.Mono.1" = "Mono",
  "12_CD14.Mono.2" = "Mono",
  "13_CD16.Mono" = "Mono",
  "15_CLP.2" = "CLP",
  "16_Pre.B" = "PreB",
  "17_B" = "B",
  "18_Plasma" = "Plasma",
  "19_CD8.N" = "CD8.N",
  "20_CD4.N1" = "CD4.N",
  "21_CD4.N2" = "CD4.N",
  "22_CD4.M" = "CD4.M",
  "23_CD8.EM" = "CD8.EM",
  "24_CD8.CM" = "CD8.CM",
  "25_NK" = "NK"
)
remapClust <- remapClust[names(remapClust) %in% labelNew]
labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters, newLabels = labelNew2, oldLabels = labelOld)
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Clusters2")
pathToMacs2 <- findMacs2()

# save workspace
# save.image(file = "archR_workflow_before_peaks.RData")
#

# load workspace
# load("archR_workflow_before_peaks.RData")

print("=========")
print("==peaks start==")
print(pathToMacs2)
projHeme4 <- addReproduciblePeakSet(
  ArchRProj = projHeme4, 
  groupBy = "Clusters2", 
  pathToMacs2 = findMacs2()
)
print("==peaks end==")
getPeakSet(projHeme4)
projHeme5 <- addPeakMatrix(projHeme4)


