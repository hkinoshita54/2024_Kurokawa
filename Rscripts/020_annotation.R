####################
# description ----
# proceed from 010_reanalysis_rda.R, lognorm, npc=30, res=1
# manual annotation dependent on marker expression


####################
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SeuratWrappers)
options(Seurat.object.assay.version = "v5")


####################
# load data, set parameters ----
seu <- readRDS(file = "RDSfiles/seu_reanalysis_rda_010_lognorm_npc30.RDS")
description = "annotation"
file_no = "020"
path = file.path("plots", description)
dir.create(path = path, recursive = TRUE)
RDSfile = paste0("RDSfiles/", "seu_", file_no, "_", description, ".RDS")


####################
# annotation ----
Idents(seu) <- "seurat_clusters"
seu <- RenameIdents(
  seu, 
  `0` = "Parietal", `1` = "Parietal", `2` = "Parietal", `3` = "Pit", `4` = "Parietal", 
  `5` = "Parietal", `6` = "Neck", `7` = "Isthmus", `8` = "Tuft", `9` = "RBC", 
  `10` = "Endocrine", `11` = "Chief", `12` ="Immune"
)
seu$celltype <- Idents(seu)
seu$celltype <- factor(seu$celltype, levels = c("Isthmus", "Pit", "Neck", "Chief", "Parietal", "Tuft", "Endocrine", "RBC", "Immune"))
Idents(seu) <- seu$celltype
DimPlot(seu, label = TRUE, repel = TRUE) + NoAxes()
ggsave("annotation.png", path = path, width = 6, height = 5, units = "in", dpi = 150)
DimPlot(seu, split.by = "orig.ident") + NoAxes()
ggsave("split.png", path = path, width = 16, height = 4, units = "in", dpi = 150)


####################
# recheck feature plot ----
features = readLines("gene_set/04_epi_markers_2.txt")
for(i in 1:length(features)){
  tryCatch({
    p <- FeaturePlot(seu, features = features[i], cols = c("lightgrey","darkred"), label = TRUE, repel = TRUE) + 
      NoAxes() + NoLegend()
    ggsave(paste0(as.character(features[i]), ".png"), 
           plot = p, path = path, 
           width = 5, height = 5, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}


####################
# check dot plot ----
features = readLines("gene_set/04_epi_markers_2.txt")
DotPlot(seu, features = features) + RotatedAxis()
ggsave("dotplot.png", path = path, width = 5, height = 4, units = "in", dpi = 150) 


####################
# organize meta data and save ----
seu$orig.ident <- as.character(seu$orig.ident)
seu$orig.ident[seu$TagIDs == "A951"] <- "WT"
seu$orig.ident[seu$TagIDs == "A952"] <- "Lgr4KO"
seu$orig.ident[seu$TagIDs == "A953"] <- "Rspo3OE"
seu$orig.ident[seu$TagIDs == "A954"] <- "L+R"
seu$orig.ident <- factor(seu$orig.ident, levels = c("WT", "Lgr4KO", "Rspo3OE", "L+R"))
seu@meta.data <- seu@meta.data[,c(1:4, 11:13, 26, 27, 31)]
saveRDS(seu, file = RDSfile)
