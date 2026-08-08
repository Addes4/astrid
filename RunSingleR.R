options(warn = -1)
suppressMessages(library(tidyverse))
suppressMessages(library(Matrix))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(SingleR))

args <- commandArgs(TRUE)

input_adata <- args[1]
output_csv <- args[2]
species <- args[3]

pseudobulk_matrix <- read.csv(input_adata, row.names=1)

reference_file <- ifelse(length(args) >= 4, args[4], "")

if (reference_file != "" && file.exists(reference_file)) {
  newRef <- readRDS(reference_file)
} else if (species == "human") {
  excludeSamples <- read.table(file = "data/excludeSingleR.csv", header = T, sep = ",")
  newRef <- readRDS("data/ASTRID_SingleR_Reference_20240701.Rds")
  if(nrow(excludeSamples > 0)){
    newRef <- newRef[,!(colnames(newRef) %in% excludeSamples$samples)]
  }
} else if (species == "mouse") {
  newRef <- readRDS("data/ASTRID_SingleR_Reference_mouse_20260116.Rds")
} else {
  print("Invalid species entered!")
  return(0)
}

newRefCounts <- assays(newRef)[["counts"]]

pseudobulk_matrix <- t(pseudobulk_matrix)

common_genes <- intersect(rownames(pseudobulk_matrix), rownames(newRefCounts))

pseudobulk_matrix <- pseudobulk_matrix[common_genes,]
pseudobulk_matrix <- sweep(pseudobulk_matrix,2,colSums(pseudobulk_matrix), `/`)
pseudobulk_matrix <- log10(pseudobulk_matrix *1e5 + 1)
colnames(pseudobulk_matrix) <- gsub("clustering\\_level\\_2\\.", replacement = "", colnames(pseudobulk_matrix))

newRefCounts <- newRefCounts[common_genes,]
newRefCounts <- sweep(newRefCounts,2,colSums(newRefCounts), `/`)
newRefCounts <- log10(newRefCounts *1e5 + 1)

predictions <- SingleR(pseudobulk_matrix, ref = newRefCounts, labels = newRef$cellTypeGRINT, de.n = 50)

write.table(predictions, file = output_csv, sep = ",", quote = F)
