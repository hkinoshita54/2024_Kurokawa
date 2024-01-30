####################
# description ----
# reanalysis of scRNA-seq of Tff1-Rspo3OE-Lgr4ff with BD Rhapsody (not 10xgenomics)
# using rda file provided by the company (ImmunoGeneTeqs)


####################
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")


####################
# load data ----
load("~/WORKSPACE/2024_Kurokawa/RDSfiles/HayakawaWTA1_Seurat.rda")    # contains a seurat object


####################
# set conditions ----
sample_name = "reanalysis_rda"
file_no = "010"
npcs = 30

path = paste0("plots/", sample_name, "/", file_no, "_lognorm/", "npc", as.character(npcs))
dir.create(path, recursive = TRUE)
RDSfile = paste0("RDSfiles/seu_", sample_name, "_", file_no, "_lognorm_npc", as.character(npcs), ".RDS")


####################
# convert the seurat onject to v5 ----
seu <- CreateSeuratObject(counts = mBC@raw.data, min.cells = 3, min.features = 200)
seu[["RNA"]]$data <- mBC@data
seu[["RNA"]]$scale.data <- mBC@scale.data
seu@meta.data <- mBC@meta.data


####################
# cluster without integration ----
seu <- FindVariableFeatures(seu)    # use data and scale.data already calculated
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("cluster.png", path = path, width = 5, height = 5, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes()
ggsave("id.png", path = path, width = 5, height = 5, units = "in", dpi = 150)


####################
#### feature plots ----
files <- list.files(path = "gene_set/annotation/", full.names = TRUE)
features <- lapply(files, FUN = readLines) %>% unlist()
for(i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) +
      NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), plot = p, 
           path = path, 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

feature = "Fzd5"
FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + NoAxes() + NoLegend()
ggsave(paste0(feature, ".png"), path = path, width = 5, height = 5, units = "in", dpi = 150)


####################
#### save seurat object ----
saveRDS(seu, file = RDSfile)
