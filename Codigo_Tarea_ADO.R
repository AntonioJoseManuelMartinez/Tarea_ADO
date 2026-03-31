
# 1) Cargar paquetes
library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(AnnotationDbi)

# 2) Leer la tabla DESeq2 de Galaxy
res <- read.delim(
  "Galaxy66-[DESeq2 result file on dataset 57-62].tabular",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 3) Poner nombres de columnas
colnames(res) <- c(
  "ENTREZID",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

# 4) Convertir tipos
res$ENTREZID <- as.character(res$ENTREZID)
res$baseMean <- as.numeric(res$baseMean)
res$log2FoldChange <- as.numeric(res$log2FoldChange)
res$lfcSE <- as.numeric(res$lfcSE)
res$stat <- as.numeric(res$stat)
res$pvalue <- as.numeric(res$pvalue)
res$padj <- as.numeric(res$padj)

# 5) Quitar NA en padj

res_clean <- res %>%
  filter(!is.na(padj))

# 6) Anotar IDs a símbolo y nombre
annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(res_clean$ENTREZID),
  columns = c("SYMBOL", "GENENAME"),
  keytype = "ENTREZID"
) %>%
  distinct(ENTREZID, .keep_all = TRUE)
res_annot <- res_clean %>%
  left_join(annot, by = "ENTREZID")

# 7) Filtrar genes diferencialmente expresados
deg <- res_annot %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

# 8) Separar genes up en control y up en ZCCHC8
up_control <- deg %>%
  filter(log2FoldChange > 1) %>%
  arrange(desc(log2FoldChange))

up_zcchc8 <- deg %>%
  filter(log2FoldChange < -1) %>%
  arrange(log2FoldChange)

# 9) Resumen rápido
cat("Genes significativos (padj < 0.05 y |log2FC| > 1):", nrow(deg), "\n")
cat("Up en control:", nrow(up_control), "\n")
cat("Up en ZCCHC8:", nrow(up_zcchc8), "\n")

# 10) Guardar tablas
write.table(res_annot, "DESeq2_results_annotated.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(deg, "DEG_padj0.05_abslog2FC1.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(up_control, "DEG_up_in_control.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(up_zcchc8, "DEG_up_in_ZCCHC8.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 11) Preparar volcano plot
volcano_df <- res_annot %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up_control",
      padj < 0.05 & log2FoldChange < -1 ~ "Up_ZCCHC8",
      TRUE ~ "Not_sig"
    ),
    negLog10Padj = -log10(padj)
  )


# 12) Elegir genes a etiquetar
top_labels <- volcano_df %>%
  filter(significance != "Not_sig") %>%
  arrange(padj) %>%
  slice_head(n = 15)


# 13) Crear volcano plot
p_volcano <- ggplot(volcano_df, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = significance), alpha = 0.7, size = 1.6) +
  scale_color_manual(values = c(
    "Up_control" = "steelblue3",
    "Up_ZCCHC8" = "firebrick3",
    "Not_sig" = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text(
    data = top_labels,
    aes(label = ifelse(is.na(SYMBOL), ENTREZID, SYMBOL)),
    size = 3,
    vjust = -0.3,
    check_overlap = TRUE
  ) +
  labs(
    title = "Volcano plot: control vs ZCCHC8",
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Grupo"
  ) +
  theme_minimal(base_size = 12)

# 14) Guardar volcano plot
ggsave("volcano_plot_control_vs_ZCCHC8.png",
       p_volcano, width = 9, height = 7, dpi = 300)

# 15) Crear listas para enriquecimiento GO BP
genes_control <- unique(up_control$ENTREZID)
genes_zcchc8 <- unique(up_zcchc8$ENTREZID)
gene_universe <- unique(res_clean$ENTREZID)

# 16) GO Biological Process para control
ego_control_bp <- enrichGO(
  gene = genes_control,
  universe = gene_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# 17) GO Biological Process para ZCCHC8
ego_zcchc8_bp <- enrichGO(
  gene = genes_zcchc8,
  universe = gene_universe,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# 18) Guardar tablas GO BP
write.table(as.data.frame(ego_control_bp), "GO_BP_up_control.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(as.data.frame(ego_zcchc8_bp), "GO_BP_up_ZCCHC8.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 19) Dotplot GO BP control

png("GO_BP_dotplot_up_control.png", width = 2200, height = 1600, res = 220)
print(dotplot(ego_control_bp, showCategory = 15) +
        ggtitle("GO BP - Up in control"))
dev.off()

ego_control_bp_s <- simplify(ego_control_bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

png("GO_BP_dotplot_up_control_simplified.png", width = 2200, height = 1600, res = 220)
print(dotplot(ego_control_bp_s, showCategory = 12) + ggtitle("GO BP - Up in control"))
dev.off()

# 20) Dotplot GO BP ZCCHC8
png("GO_BP_dotplot_up_ZCCHC8.png", width = 2200, height = 1600, res = 220)
print(dotplot(ego_zcchc8_bp, showCategory = 15) +
        ggtitle("GO BP - Up in ZCCHC8"))
dev.off()

# 21) Top 20 genes para el informe
top_control <- up_control %>%
  dplyr::select(ENTREZID, SYMBOL, GENENAME, log2FoldChange, padj) %>%
  dplyr::slice_head(n = 20)

top_zcchc8 <- up_zcchc8 %>%
  dplyr::select(ENTREZID, SYMBOL, GENENAME, log2FoldChange, padj) %>%
  dplyr::slice_head(n = 20)







