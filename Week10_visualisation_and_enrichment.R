# Week 10: From a gene list to biology
# Input: DESeq2 results and normalised expression from Weeks 8 and 9.
# Comparison: thin biofilm (45 days) vs early biofilm (10 days).
# Author: Nicole Navarro
# Date: 2026-10-01
library(ggplot2)
library(pheatmap)
res_df <- read.csv("Week9_results/deseq2_all_genes.csv")
log_cpm <- readRDS("Week8_results/log_cpm.rds")
sample_metadata <- readRDS("Week8_results/sample_metadata.rds")
head(res_df)
dim(res_df)
str(res_df)
levels(sample_metadata$condition)
res_df$neg_log10_padj <- -log10(res_df$padj)
head(res_df[, c("gene_id", "log2FoldChange", "padj", "neg_log10_padj")])
res_df$status <- "Not significant"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up in thin"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down in thin"
table(res_df$status)
p_volcano <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, colour = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(values = c(
    "Up in thin" = "#d73027",
    "Down in thin" = "#4575b4",
    "Not significant" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey40") +
  labs(
    title = "Differential expression: thin versus early biofilm",
    x = "log2 fold change (thin / early)",
    y = "-log10 adjusted p-value",
    colour = ""
  ) +
  theme_bw()

p_volcano
top10 <- res_df[order(res_df$padj), ][1:10, ]

p_volcano +
  geom_text(
    data = top10,
    aes(label = gene_id),
    colour = "black",
    size = 2.5,
    vjust = -0.8
  )
dir.create("Week10_results", showWarnings = FALSE)
ggsave("Week10_results/volcano_plot.png", p_volcano, width = 7, height = 6, dpi = 300)
res_df$status_strict <- "Not significant"
res_df$status_strict[res_df$padj < 0.01 & res_df$log2FoldChange > 1] <- "Up in thin"
res_df$status_strict[res_df$padj < 0.01 & res_df$log2FoldChange < -1] <- "Down in thin"
table(res_df$status_strict)
table(res_df$status)
sig_ids <- res_df$gene_id[res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1]
sig_ids <- sig_ids[!is.na(sig_ids)]
length(sig_ids)

heat_matrix <- log_cpm[rownames(log_cpm) %in% sig_ids, ]
dim(heat_matrix)
annotation_col <- data.frame(Stage = sample_metadata$condition)
rownames(annotation_col) <- sample_metadata$sample_id
annotation_col

pheatmap(
  heat_matrix,
  scale = "row",
  show_rownames = FALSE,
  annotation_col = annotation_col,
  main = "Differentially expressed genes: thin vs early"
)
pheatmap(heat_matrix, scale = "none", show_rownames = FALSE, main = "Unscaled")
pheatmap(heat_matrix, scale = "row", show_rownames = FALSE, main = "Row scaled")
library(org.Sc.sgd.db)
columns(org.Sc.sgd.db)
res_df$gene_name <- mapIds(
  org.Sc.sgd.db,
  keys = res_df$gene_id,
  column = "COMMON",
  keytype = "ORF",
  multiVals = "first"
)
head(res_df[, c("gene_id", "gene_name", "log2FoldChange", "padj")], 10)

res_df$description <- mapIds(
  org.Sc.sgd.db,
  keys = res_df$gene_id,
  column = "DESCRIPTION",
  keytype = "ORF",
  multiVals = "first"
)
head(res_df$description, 3)
head(res_df[order(res_df$padj), c("gene_id", "gene_name", "log2FoldChange")], 5)
library(clusterProfiler)

gene_universe <- res_df$gene_id[!is.na(res_df$padj)]
length(gene_universe)
length(sig_ids)
ego <- enrichGO(
  gene = sig_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)
head(as.data.frame(ego)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")], 10)
dotplot(ego, showCategory = 15) +
  labs(title = "Enriched biological processes: thin vs early")
up_ids <- res_df$gene_id[res_df$padj < 0.05 & res_df$log2FoldChange > 1]
up_ids <- up_ids[!is.na(up_ids)]

ego_up <- enrichGO(
  gene = up_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pvalueCutoff = 0.05
)
head(as.data.frame(ego_up)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")], 10)
down_ids <- res_df$gene_id[res_df$padj < 0.05 & res_df$log2FoldChange < -1]
down_ids <- down_ids[!is.na(down_ids)]

ego_down <- enrichGO(
  gene = down_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pvalueCutoff = 0.05
)
head(as.data.frame(ego_down)[, c("Description", "GeneRatio", "BgRatio", "p.adjust")], 10)
themes <- "carbohydrate|carbon|glycoly|stress|cell wall"
as.data.frame(ego_up)$Description[grepl(themes, as.data.frame(ego_up)$Description, ignore.case = TRUE)]
as.data.frame(ego_down)$Description[grepl(themes, as.data.frame(ego_down)$Description, ignore.case = TRUE)]
themes <- "carbohydrate|carbon|glycoly|stress|cell wall"
as.data.frame(ego_up)$Description[grepl(themes, as.data.frame(ego_up)$Description, ignore.case = TRUE)]
as.data.frame(ego_down)$Description[grepl(themes, as.data.frame(ego_down)$Description, ignore.case = TRUE)]
write.csv(res_df, "Week10_results/annotated_results.csv", row.names = FALSE)
write.csv(as.data.frame(ego), "Week10_results/go_enrichment.csv", row.names = FALSE)
saveRDS(heat_matrix, "Week10_results/heat_matrix.rds")
top30_ids <- res_df$gene_id[order(res_df$padj)][1:30]
heat_top30 <- log_cpm[rownames(log_cpm) %in% top30_ids, ]
dim(heat_top30)

row_labels <- res_df$gene_name[match(rownames(heat_top30), res_df$gene_id)]
row_labels[is.na(row_labels)] <- rownames(heat_top30)[is.na(row_labels)]

p_heat30 <- pheatmap(
  heat_top30,
  scale = "row",
  show_rownames = TRUE,
  labels_row = row_labels,
  annotation_col = annotation_col,
  main = "Top 30 DE genes: thin vs early"
)

png("Week10_results/heatmap_top30.png", width = 7, height = 9, units = "in", res = 300)
grid::grid.draw(p_heat30$gtable)
dev.off()
sessionInfo()
head(as.data.frame(ego_up)$Description, 10)
head(as.data.frame(ego_down)$Description, 10)
as.data.frame(ego)$Description[grepl("carbohydrate|carbon|glycoly|stress|cell wall", as.data.frame(ego)$Description, ignore.case = TRUE)]
grep("cell wall", as.data.frame(ego_up)$Description, value = TRUE)
grep("cell wall", as.data.frame(ego_down)$Description, value = TRUE)
