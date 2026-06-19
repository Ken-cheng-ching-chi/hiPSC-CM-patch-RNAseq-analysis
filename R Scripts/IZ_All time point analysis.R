source("load_and_metadata.R")

# ============================================================
# IZ all-time-point analysis
# Design: control time effect, test overall patch vs control effect
# Note:
# This model estimates the average patch effect across IZ Day3, 1wk, and 4wks,
# after controlling for time. It is NOT a time-point-specific comparison.
# ============================================================


# ============================================================
# 1. Subset IZ samples
# ============================================================

iz_samples <- rownames(coldata)[coldata$region == "IZ"]

counts_iz <- counts[, iz_samples]
coldata_iz <- coldata[iz_samples, ]

# Drop unused factor levels
coldata_iz$region <- droplevels(coldata_iz$region)
coldata_iz$time <- droplevels(coldata_iz$time)
coldata_iz$treatment <- droplevels(coldata_iz$treatment)

# Check sample number
table(coldata_iz$time, coldata_iz$treatment)


# ============================================================
# 2. DESeq2 analysis: IZ all time points
# Design: control time effect, test patch vs control
# ============================================================

dds_iz_all <- DESeqDataSetFromMatrix(
  countData = counts_iz,
  colData = coldata_iz,
  design = ~ time + treatment
)

# Filter low-count genes
dds_iz_all <- dds_iz_all[rowSums(counts(dds_iz_all)) >= 1, ]

dim(dds_iz_all)

# Run DESeq2
dds_iz_all <- DESeq(dds_iz_all)


# ============================================================
# 3. Extract DESeq2 result: IZ all time points patch vs control
# ============================================================

res_iz_all <- results(
  dds_iz_all,
  contrast = c("treatment", "patch", "control")
)

# Sort by adjusted p-value
res_iz_all <- res_iz_all[order(res_iz_all$padj), ]

# Check result
summary(res_iz_all)
head(res_iz_all)
class(res_iz_all)

# Convert DESeq2 result to a normal data frame
res_iz_all_df <- as.data.frame(res_iz_all)

# Add gene_id as a normal column
res_iz_all_df$gene_id <- rownames(res_iz_all_df)

# Reorder columns
res_iz_all_df <- res_iz_all_df[, c(
  "gene_id",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

View(res_iz_all_df)


# ============================================================
# 4. Define DEG thresholds: IZ all time points
# ============================================================

# More conservative DEG threshold
# padj < 0.05 and |log2FC| > 1
padj_cutoff <- 0.05
padj_log2fc_cutoff <- 1

# Original paper-like DEG threshold
# raw p-value < 0.05 and fold-change > 1.5
# Equivalent in DESeq2:
# pvalue < 0.05 and |log2FoldChange| > log2(1.5)
pvalue_cutoff <- 0.05
pvalue_log2fc_cutoff <- log2(1.5)


# ============================================================
# 5. Count and export DEG results: padj threshold
# ============================================================

sig_deg_iz_all_padj <- res_iz_all_df[
  !is.na(res_iz_all_df$padj) &
    res_iz_all_df$padj < padj_cutoff &
    abs(res_iz_all_df$log2FoldChange) > padj_log2fc_cutoff,
]

sig_deg_iz_all_padj <- sig_deg_iz_all_padj[
  order(sig_deg_iz_all_padj$padj),
]

sig_up_iz_all_padj <- sig_deg_iz_all_padj[
  sig_deg_iz_all_padj$log2FoldChange > 0,
]

sig_down_iz_all_padj <- sig_deg_iz_all_padj[
  sig_deg_iz_all_padj$log2FoldChange < 0,
]

nrow(sig_deg_iz_all_padj)
nrow(sig_up_iz_all_padj)
nrow(sig_down_iz_all_padj)


# ============================================================
# 6. Count and export DEG results: pvalue threshold
# ============================================================

sig_deg_iz_all_pvalue <- res_iz_all_df[
  !is.na(res_iz_all_df$pvalue) &
    res_iz_all_df$pvalue < pvalue_cutoff &
    abs(res_iz_all_df$log2FoldChange) > pvalue_log2fc_cutoff,
]

sig_deg_iz_all_pvalue <- sig_deg_iz_all_pvalue[
  order(sig_deg_iz_all_pvalue$pvalue),
]

sig_up_iz_all_pvalue <- sig_deg_iz_all_pvalue[
  sig_deg_iz_all_pvalue$log2FoldChange > 0,
]

sig_down_iz_all_pvalue <- sig_deg_iz_all_pvalue[
  sig_deg_iz_all_pvalue$log2FoldChange < 0,
]

nrow(sig_deg_iz_all_pvalue)
nrow(sig_up_iz_all_pvalue)
nrow(sig_down_iz_all_pvalue)


# ============================================================
# 7. Export DESeq2 tables and summary tables
# ============================================================

# Export all DESeq2 results
write.csv(
  res_iz_all_df,
  "DESeq2_IZ_all_time_patch_vs_control_all_genes.csv",
  row.names = FALSE
)

# Export padj DEGs
write.csv(
  sig_deg_iz_all_padj,
  "DESeq2_IZ_all_time_patch_vs_control_DEGs_padj0.05_log2FC1.csv",
  row.names = FALSE
)

# Export pvalue DEGs
write.csv(
  sig_deg_iz_all_pvalue,
  "DESeq2_IZ_all_time_patch_vs_control_DEGs_pvalue0.05_FC1.5.csv",
  row.names = FALSE
)

# Summary table: padj threshold
deg_summary_iz_all_padj <- data.frame(
  comparison = "IZ all time points patch vs control",
  design = "~ time + treatment",
  filter_threshold = "rowSums >= 1",
  DEG_threshold = "padj < 0.05 and |log2FC| > 1",
  total_tested_genes = nrow(res_iz_all_df),
  significant_DEGs = nrow(sig_deg_iz_all_padj),
  upregulated_in_patch = nrow(sig_up_iz_all_padj),
  downregulated_in_patch = nrow(sig_down_iz_all_padj)
)

# Summary table: pvalue threshold
deg_summary_iz_all_pvalue <- data.frame(
  comparison = "IZ all time points patch vs control",
  design = "~ time + treatment",
  filter_threshold = "rowSums >= 1",
  DEG_threshold = "raw pvalue < 0.05 and |log2FC| > log2(1.5)",
  total_tested_genes = nrow(res_iz_all_df),
  significant_DEGs = nrow(sig_deg_iz_all_pvalue),
  upregulated_in_patch = nrow(sig_up_iz_all_pvalue),
  downregulated_in_patch = nrow(sig_down_iz_all_pvalue)
)

# Combined summary table
deg_summary_iz_all <- rbind(
  deg_summary_iz_all_padj,
  deg_summary_iz_all_pvalue
)

View(deg_summary_iz_all)

write.csv(
  deg_summary_iz_all,
  "DESeq2_IZ_all_time_patch_vs_control_summary_padj_and_pvalue.csv",
  row.names = FALSE
)


# ============================================================
# 8. PCA plot: IZ all time points
# ============================================================

library(ggplot2)

vsd_iz_all <- vst(dds_iz_all, blind = FALSE)

# Treatment-based PCA
pca_iz_all_treatment <- plotPCA(vsd_iz_all, intgroup = "treatment") +
  ggtitle("IZ all time points: PCA by treatment")

pca_iz_all_treatment

ggsave(
  filename = "PCA_IZ_all_time_patch_vs_control_by_treatment.bmp",
  plot = pca_iz_all_treatment,
  width = 7,
  height = 5,
  dpi = 600
)

# Time + treatment PCA
pca_iz_all_time_treatment <- plotPCA(vsd_iz_all, intgroup = c("time", "treatment")) +
  ggtitle("IZ all time points: PCA by time and treatment")

pca_iz_all_time_treatment

ggsave(
  filename = "PCA_IZ_all_time_patch_vs_control_by_time_treatment.bmp",
  plot = pca_iz_all_time_treatment,
  width = 7,
  height = 5,
  dpi = 600
)


# ============================================================
# 9. Volcano plot: IZ all time points using padj threshold
# ============================================================

volcano_iz_all_padj_df <- res_iz_all_df

volcano_iz_all_padj_df$significant <- ifelse(
  !is.na(volcano_iz_all_padj_df$padj) &
    volcano_iz_all_padj_df$padj < padj_cutoff &
    abs(volcano_iz_all_padj_df$log2FoldChange) > padj_log2fc_cutoff,
  "padj < 0.05 & |log2FC| > 1",
  "Not significant"
)

volcano_iz_all_padj_df_clean <- volcano_iz_all_padj_df[
  !is.na(volcano_iz_all_padj_df$pvalue) &
    !is.na(volcano_iz_all_padj_df$padj) &
    !is.na(volcano_iz_all_padj_df$log2FoldChange) &
    is.finite(volcano_iz_all_padj_df$pvalue) &
    is.finite(volcano_iz_all_padj_df$padj) &
    is.finite(volcano_iz_all_padj_df$log2FoldChange),
]

volcano_iz_all_padj_df_clean$significant <- factor(
  volcano_iz_all_padj_df_clean$significant,
  levels = c("padj < 0.05 & |log2FC| > 1", "Not significant")
)

volcano_iz_all_padj_plot <- ggplot(
  volcano_iz_all_padj_df_clean,
  aes(x = log2FoldChange, y = -log10(pvalue), color = significant)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "padj < 0.05 & |log2FC| > 1" = "#F8766D",
      "Not significant" = "darkturquoise"
    ),
    drop = FALSE
  ) +
  theme_bw() +
  labs(
    title = "IZ all time points: patch vs control",
    subtitle = "DESeq2 design: ~ time + treatment; padj threshold",
    x = "log2 fold change",
    y = "-log10 p-value",
    color = "significant"
  )

volcano_iz_all_padj_plot

ggsave(
  filename = "Volcano_IZ_all_time_patch_vs_control_padj0.05_log2FC1.bmp",
  plot = volcano_iz_all_padj_plot,
  width = 7,
  height = 5,
  dpi = 600
)


# ============================================================
# 10. Volcano plot: IZ all time points using pvalue threshold
# ============================================================

volcano_iz_all_pvalue_df <- res_iz_all_df

volcano_iz_all_pvalue_df$significant <- ifelse(
  !is.na(volcano_iz_all_pvalue_df$pvalue) &
    volcano_iz_all_pvalue_df$pvalue < pvalue_cutoff &
    abs(volcano_iz_all_pvalue_df$log2FoldChange) > pvalue_log2fc_cutoff,
  "raw p < 0.05 & FC > 1.5",
  "Not significant"
)

volcano_iz_all_pvalue_df_clean <- volcano_iz_all_pvalue_df[
  !is.na(volcano_iz_all_pvalue_df$pvalue) &
    !is.na(volcano_iz_all_pvalue_df$log2FoldChange) &
    is.finite(volcano_iz_all_pvalue_df$pvalue) &
    is.finite(volcano_iz_all_pvalue_df$log2FoldChange),
]

volcano_iz_all_pvalue_df_clean$significant <- factor(
  volcano_iz_all_pvalue_df_clean$significant,
  levels = c("raw p < 0.05 & FC > 1.5", "Not significant")
)

volcano_iz_all_pvalue_plot <- ggplot(
  volcano_iz_all_pvalue_df_clean,
  aes(x = log2FoldChange, y = -log10(pvalue), color = significant)
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "raw p < 0.05 & FC > 1.5" = "#F8766D",
      "Not significant" = "darkturquoise"
    ),
    drop = FALSE
  ) +
  theme_bw() +
  labs(
    title = "IZ all time points: patch vs control",
    subtitle = "DESeq2 design: ~ time + treatment; raw p-value threshold",
    x = "log2 fold change",
    y = "-log10 p-value",
    color = "significant"
  )

volcano_iz_all_pvalue_plot

ggsave(
  filename = "Volcano_IZ_all_time_patch_vs_control_pvalue0.05_FC1.5.bmp",
  plot = volcano_iz_all_pvalue_plot,
  width = 7,
  height = 5,
  dpi = 600
)


# ============================================================
# 11. Optional: check software versions
# ============================================================

packageVersion("DESeq2")
R.version.string


# ============================================================
# 12. GO / pathway enrichment analysis: IZ all time points
# ============================================================

# Required packages
# If not installed, run these lines once:
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Rn.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Rn.eg.db)
library(enrichplot)
library(ggplot2)


# ============================================================
# 12-1. Prepare gene ID conversion
# ============================================================

# Remove possible Ensembl version suffix if present
res_iz_all_df$gene_id_clean <- gsub("\\..*$", "", res_iz_all_df$gene_id)

sig_deg_iz_all_pvalue$gene_id_clean <- gsub("\\..*$", "", sig_deg_iz_all_pvalue$gene_id)
sig_deg_iz_all_padj$gene_id_clean <- gsub("\\..*$", "", sig_deg_iz_all_padj$gene_id)

# Background genes: all tested genes
background_genes <- unique(res_iz_all_df$gene_id_clean)

background_map <- bitr(
  background_genes,
  fromType = "ENSEMBL",
  toType = c("ENTREZID", "SYMBOL"),
  OrgDb = org.Rn.eg.db
)

background_entrez <- unique(background_map$ENTREZID)

write.csv(
  background_map,
  "IZ_all_time_background_ENSEMBL_to_ENTREZ_SYMBOL.csv",
  row.names = FALSE
)


# ============================================================
# 12-2. Helper function: convert DEG table to Entrez IDs
# ============================================================

convert_deg_to_entrez <- function(deg_df) {
  
  if (is.null(deg_df) || nrow(deg_df) == 0) {
    message("DEG table is empty. Skip ID conversion.")
    return(list(
      map = data.frame(),
      entrez = character(0)
    ))
  }
  
  deg_genes <- unique(deg_df$gene_id_clean)
  deg_genes <- deg_genes[!is.na(deg_genes) & deg_genes != ""]
  
  if (length(deg_genes) == 0) {
    message("No valid gene_id_clean found. Skip ID conversion.")
    return(list(
      map = data.frame(),
      entrez = character(0)
    ))
  }
  
  valid_ensembl <- keys(org.Rn.eg.db, keytype = "ENSEMBL")
  deg_genes_valid <- intersect(deg_genes, valid_ensembl)
  
  message("Input DEG genes: ", length(deg_genes))
  message("Valid ENSEMBL genes in org.Rn.eg.db: ", length(deg_genes_valid))
  
  if (length(deg_genes_valid) == 0) {
    message("None of the DEG gene IDs are valid ENSEMBL IDs in org.Rn.eg.db. Skip ID conversion.")
    return(list(
      map = data.frame(),
      entrez = character(0)
    ))
  }
  
  deg_map <- bitr(
    deg_genes_valid,
    fromType = "ENSEMBL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Rn.eg.db
  )
  
  deg_entrez <- unique(deg_map$ENTREZID)
  
  return(list(
    map = deg_map,
    entrez = deg_entrez
  ))
}


# ============================================================
# 12-3. Helper function: run GO and KEGG enrichment
# ============================================================

run_go_kegg <- function(gene_entrez, output_prefix) {
  
  if (length(gene_entrez) < 10) {
    message(paste0(output_prefix, ": gene number < 10, skip enrichment."))
    return(NULL)
  }
  
  # ----------------------------
  # GO Biological Process
  # ----------------------------
  ego_bp <- enrichGO(
    gene = gene_entrez,
    universe = background_entrez,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  write.csv(
    as.data.frame(ego_bp),
    paste0(output_prefix, "_GO_BP.csv"),
    row.names = FALSE
  )
  
  if (nrow(as.data.frame(ego_bp)) > 0) {
    p_bp <- dotplot(ego_bp, showCategory = 10) +
      ggtitle(paste0(output_prefix, " GO Biological Process"))
    
    ggsave(
      filename = paste0(output_prefix, "_GO_BP_dotplot.bmp"),
      plot = p_bp,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
  
  
  # ----------------------------
  # GO Cellular Component
  # ----------------------------
  ego_cc <- enrichGO(
    gene = gene_entrez,
    universe = background_entrez,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  write.csv(
    as.data.frame(ego_cc),
    paste0(output_prefix, "_GO_CC.csv"),
    row.names = FALSE
  )
  
  if (nrow(as.data.frame(ego_cc)) > 0) {
    p_cc <- dotplot(ego_cc, showCategory = 10) +
      ggtitle(paste0(output_prefix, " GO Cellular Component"))
    
    ggsave(
      filename = paste0(output_prefix, "_GO_CC_dotplot.bmp"),
      plot = p_cc,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
  
  
  # ----------------------------
  # GO Molecular Function
  # ----------------------------
  ego_mf <- enrichGO(
    gene = gene_entrez,
    universe = background_entrez,
    OrgDb = org.Rn.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  write.csv(
    as.data.frame(ego_mf),
    paste0(output_prefix, "_GO_MF.csv"),
    row.names = FALSE
  )
  
  if (nrow(as.data.frame(ego_mf)) > 0) {
    p_mf <- dotplot(ego_mf, showCategory = 10) +
      ggtitle(paste0(output_prefix, " GO Molecular Function"))
    
    ggsave(
      filename = paste0(output_prefix, "_GO_MF_dotplot.bmp"),
      plot = p_mf,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
  
  
  # ----------------------------
  # KEGG pathway
  # ----------------------------
  ekegg <- tryCatch(
    {
      enrichKEGG(
        gene = gene_entrez,
        universe = background_entrez,
        organism = "rno",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
    },
    error = function(e) {
      message("KEGG enrichment failed. This may be due to internet/API issues.")
      return(NULL)
    }
  )
  
  if (!is.null(ekegg)) {
    write.csv(
      as.data.frame(ekegg),
      paste0(output_prefix, "_KEGG.csv"),
      row.names = FALSE
    )
    
    if (nrow(as.data.frame(ekegg)) > 0) {
      p_kegg <- dotplot(ekegg, showCategory = 10) +
        ggtitle(paste0(output_prefix, " KEGG pathway"))
      
      ggsave(
        filename = paste0(output_prefix, "_KEGG_dotplot.bmp"),
        plot = p_kegg,
        width = 10,
        height = 6,
        dpi = 300
      )
    }
  }
  
  return(list(
    GO_BP = ego_bp,
    GO_CC = ego_cc,
    GO_MF = ego_mf,
    KEGG = ekegg
  ))
}


# ============================================================
# 12-4. IZ all-time pvalue DEG list enrichment
# raw pvalue < 0.05 and FC > 1.5
# ============================================================

iz_all_pvalue_conv <- convert_deg_to_entrez(sig_deg_iz_all_pvalue)

write.csv(
  iz_all_pvalue_conv$map,
  "IZ_all_time_pvalue_DEG_ENSEMBL_to_ENTREZ_SYMBOL.csv",
  row.names = FALSE
)

length(iz_all_pvalue_conv$entrez)

enrich_iz_all_pvalue <- run_go_kegg(
  gene_entrez = iz_all_pvalue_conv$entrez,
  output_prefix = "IZ_all_time_pvalue_DEG"
)


# ============================================================
# 12-5. IZ all-time padj DEG list enrichment
# padj < 0.05 and |log2FC| > 1
# ============================================================

iz_all_padj_conv <- convert_deg_to_entrez(sig_deg_iz_all_padj)

write.csv(
  iz_all_padj_conv$map,
  "IZ_all_time_padj_DEG_ENSEMBL_to_ENTREZ_SYMBOL.csv",
  row.names = FALSE
)

length(iz_all_padj_conv$entrez)

enrich_iz_all_padj <- run_go_kegg(
  gene_entrez = iz_all_padj_conv$entrez,
  output_prefix = "IZ_all_time_padj_DEG"
)


# ============================================================
# 12-6. Quick check: top GO BP results
# ============================================================

# pvalue DEG list top GO BP
if (!is.null(enrich_iz_all_pvalue)) {
  head(as.data.frame(enrich_iz_all_pvalue$GO_BP), 20)
}

# padj DEG list top GO BP
if (!is.null(enrich_iz_all_padj)) {
  head(as.data.frame(enrich_iz_all_padj$GO_BP), 20)
}


# ============================================================
# 13. Check TGF-beta / fibrosis / ECM-related genes in DESeq2 results
# ============================================================

library(org.Rn.eg.db)
library(AnnotationDbi)

# Add SYMBOL to DESeq2 result
add_symbol <- function(res_df) {
  gene_id_clean <- gsub("\\..*$", "", res_df$gene_id)
  
  symbol <- mapIds(
    org.Rn.eg.db,
    keys = gene_id_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  res_df$SYMBOL <- symbol
  return(res_df)
}

res_iz_all_symbol <- add_symbol(res_iz_all_df)

# TGF-beta / fibrosis / ECM-related genes to inspect
tgfb_related_genes <- c(
  "Tgfb1", "Tgfb2", "Tgfb3",
  "Tgfbr1", "Tgfbr2",
  "Smad2", "Smad3", "Smad4", "Smad7",
  "Acta2", "Cxcl12",
  "Col1a1", "Col1a2", "Col3a1",
  "Fn1", "Postn",
  "Mmp2", "Mmp3", "Mmp7", "Mmp9", "Mmp12",
  "Timp1", "Timp2",
  "Lox", "Loxl2",
  "Ctgf", "Serpine1"
)

# Check IZ all time points
tgfb_check_iz_all <- res_iz_all_symbol[
  res_iz_all_symbol$SYMBOL %in% tgfb_related_genes,
  c("gene_id", "SYMBOL", "baseMean", "log2FoldChange", "pvalue", "padj")
]

tgfb_check_iz_all <- tgfb_check_iz_all[
  order(tgfb_check_iz_all$pvalue),
]

View(tgfb_check_iz_all)

write.csv(
  tgfb_check_iz_all,
  "IZ_all_time_TGF_beta_ECM_related_gene_check.csv",
  row.names = FALSE
)


# ============================================================
# 14. Search TGF-beta / ECM / wound-related terms in GO results
# ============================================================

search_terms <- function(enrich_result, keyword) {
  if (is.null(enrich_result)) {
    return(data.frame())
  }
  
  df <- as.data.frame(enrich_result)
  
  if (nrow(df) == 0) {
    return(df)
  }
  
  df[grepl(keyword, df$Description, ignore.case = TRUE), ]
}

# IZ all-time pvalue GO BP
if (!is.null(enrich_iz_all_pvalue)) {
  search_terms(enrich_iz_all_pvalue$GO_BP, "transforming|TGF|SMAD|fibrosis|collagen|matrix|wound")
}

# IZ all-time padj GO BP
if (!is.null(enrich_iz_all_padj)) {
  search_terms(enrich_iz_all_padj$GO_BP, "transforming|TGF|SMAD|fibrosis|collagen|matrix|wound")
}
