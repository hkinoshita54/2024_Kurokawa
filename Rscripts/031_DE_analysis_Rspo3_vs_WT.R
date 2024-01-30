####################
# description ----
# proceed from 020_annotation.R
# DE analysis, WT vs Lgr4KO, WT vs Rspo3OE, WT vs L+R


####################
# load packages ----
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(presto)
library(gprofiler2)
library(msigdbr)
library(fgsea)


####################
# load data, set parameters ----
seu <- readRDS(file = "RDSfiles/seu_020_annotation.RDS")

file_no = "030"
description = "DE"
path = file.path("results", description)
dir.create(path = path, recursive = TRUE)


####################
# prepare gene sets for GSEA ----
collections <- list()
collections$HALLMARKS <- msigdbr(species = "Mus musculus", category = "H")
collections <- lapply(collections, function(x) {
  out <- split(x = x$gene_symbol, f = x$gs_name)
})


####################
#### preparation for DE analysis ----
seu$celltype_id <- paste(seu$celltype, seu$orig.ident, sep = "_")
Idents(seu) <- "celltype_id"

# how many cells in each celltype and celltype-id
table(seu$celltype)
table(seu$celltype_id)


####################
#### DE: Isthmus Rspo3 vs WT ----

# DE analysis by using presto::wilcoxauc()
res <- wilcoxauc(X = seu[["RNA"]]$data, y = seu$celltype_id, groups_use = c("Isthmus_WT", "Isthmus_Rspo3OE"))
res <- res %>% filter(group == "Isthmus_Rspo3OE" & (pct_in != 0 | pct_out != 0))    # removing genes with no expression in both groups
write.table(res, file = file.path(path, "Isthmus_Rspo3OE_vs_WT.txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
ggplot(res, aes(logFC, -log10(pval))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, feature,"")), colour = "red", size = 3)
ranks <- res %>% select(feature, auc) %>% deframe()

# run fgsea
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = ranks, eps=0.0, minSize = 15, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, file = file.path(path, "GSEA_H_Isthmus_Rspo3_vs_WT.txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


####################
#### DE: Pit Rspo3 vs WT ----

# DE analysis by using presto::wilcoxauc()
res <- wilcoxauc(X = seu[["RNA"]]$data, y = seu$celltype_id, groups_use = c("Pit_WT", "Pit_Rspo3OE"))
res <- res %>% filter(group == "Pit_Rspo3OE" & (pct_in != 0 | pct_out != 0))    # removing genes with no expression in both groups
write.table(res, file = file.path(path, "Pit_Rspo3OE_vs_WT.txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
ggplot(res, aes(logFC, -log10(pval))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, feature,"")), colour = "red", size = 3)
ranks <- res %>% select(feature, auc) %>% deframe()

# run fgsea
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = ranks, eps=0.0, minSize = 15, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, file = file.path(path, "GSEA_H_Pit_Rspo3_vs_WT.txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


####################
#### DE: Neck Rspo3 vs WT ----

# DE analysis by using presto::wilcoxauc()
res <- wilcoxauc(X = seu[["RNA"]]$data, y = seu$celltype_id, groups_use = c("Neck_WT", "Neck_Rspo3OE"))
res <- res %>% filter(group == "Neck_Rspo3OE" & (pct_in != 0 | pct_out != 0))    # removing genes with no expression in both groups
write.table(res, file = file.path(path, "Neck_Rspo3OE_vs_WT.txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
ggplot(res, aes(logFC, -log10(pval))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.25, feature,"")), colour = "red", size = 3)
ranks <- res %>% select(feature, auc) %>% deframe()

# run fgsea
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = ranks, eps=0.0, minSize = 15, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, file = file.path(path, "GSEA_H_Neck_Rspo3_vs_WT.txt"), sep ="\t", col.names = TRUE, row.names = FALSE)


####################
#### DE: Parietal Rspo3 vs WT ----

# DE analysis by using presto::wilcoxauc()
res <- wilcoxauc(X = seu[["RNA"]]$data, y = seu$celltype_id, groups_use = c("Parietal_WT", "Parietal_Rspo3OE"))
res <- res %>% filter(group == "Parietal_Rspo3OE" & (pct_in != 0 | pct_out != 0))    # removing genes with no expression in both groups
write.table(res, file = file.path(path, "Parietal_Rspo3OE_vs_WT.txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
ggplot(res, aes(logFC, -log10(pval))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(padj < 0.01, feature,"")), colour = "red", size = 3)
ranks <- res %>% select(feature, auc) %>% deframe()

# run fgsea
fgseaRes <- fgsea(pathways = collections$HALLMARKS, stats = ranks, eps=0.0, minSize = 15, maxSize = 500)
fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.25)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()
fgseaRes$leadingEdge <- fgseaRes$leadingEdge %>% lapply(paste, collapse = ",") %>% unlist(recursive = FALSE)
write.table(fgseaRes, file = file.path(path, "GSEA_H_Parietal_Rspo3_vs_WT.txt"), sep ="\t", col.names = TRUE, row.names = FALSE)
