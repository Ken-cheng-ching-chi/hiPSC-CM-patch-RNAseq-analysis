library(dplyr)
library(purrr)
library(tidyr)
library(openxlsx)
library(org.Rn.eg.db)
library(AnnotationDbi)

# ============================================================
# Representative ECM / fibrosis / TGF-beta-related genes
# ============================================================

genes_interest <- c(
  "Col1a1", "Col1a2", "Col3a1", "Fn1", "Lox", "Postn",
  "Mmp2", "Mmp9", "Tgfb1", "Serpine1", "Tgfbr1", "Tgfbr2",
  "Smad2", "Smad3", "Timp1"
)

gene_category <- data.frame(
  SYMBOL = genes_interest,
  Category = c(
    "Collagen / ECM",
    "Collagen / ECM",
    "Collagen / ECM",
    "ECM glycoprotein",
    "Collagen crosslinking",
    "Fibrosis / repair",
    "Matrix remodeling",
    "Matrix remodeling",
    "TGF-beta signaling",
    "TGF-beta / fibrosis",
    "TGF-beta receptor",
    "TGF-beta receptor",
    "TGF-beta signaling",
    "TGF-beta signaling",
    "Matrix remodeling"
  )
)

# ============================================================
# Add SYMBOL if the dataframe does not have SYMBOL column
# ============================================================

add_symbol_if_missing <- function(res_df) {
  
  if ("SYMBOL" %in% colnames(res_df)) {
    return(res_df)
  }
  
  if (!"gene_id" %in% colnames(res_df)) {
    stop("This dataframe has no SYMBOL column and no gene_id column.")
  }
  
  gene_id_clean <- gsub("\\..*$", "", res_df$gene_id)
  
  symbol <- mapIds(
    org.Rn.eg.db,
    keys = gene_id_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  res_df$SYMBOL <- as.character(symbol)
  
  return(res_df)
}

# ============================================================
# Prepare each comparison result
# ============================================================

prep_result <- function(df, label) {
  
  df <- add_symbol_if_missing(df)
  
  df %>%
    filter(SYMBOL %in% genes_interest) %>%
    select(SYMBOL, log2FoldChange, pvalue, padj) %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    rename(
      !!paste0(label, "_log2FC") := log2FoldChange,
      !!paste0(label, "_pvalue") := pvalue,
      !!paste0(label, "_padj") := padj
    )
}

# ============================================================
# Create tables for each comparison
# ============================================================

tbl_iz_day3 <- prep_result(res_iz_day3_df, "IZ_Day3")
tbl_iz_1wk  <- prep_result(res_iz_1wk_df,  "IZ_1wk")
tbl_iz_4wks <- prep_result(res_iz_4wks_df, "IZ_4wks")
tbl_bz_1wk  <- prep_result(res_bz_1wk_df,  "BZ_1wk")

# ============================================================
# Merge all comparisons
# ============================================================

rep_gene_table <- reduce(
  list(tbl_iz_day3, tbl_iz_1wk, tbl_iz_4wks, tbl_bz_1wk),
  full_join,
  by = "SYMBOL"
) %>%
  left_join(gene_category, by = "SYMBOL") %>%
  select(SYMBOL, Category, everything()) %>%
  mutate(
    SYMBOL = factor(SYMBOL, levels = genes_interest)
  ) %>%
  arrange(SYMBOL) %>%
  mutate(SYMBOL = as.character(SYMBOL))

# View result
View(rep_gene_table)

# Export Excel
write.xlsx(
  rep_gene_table,
  "Representative_ECM_fibrosis_TGFb_genes_all_comparisons.xlsx",
  rowNames = FALSE
)

rep_gene_table








#representative ECM / fibrosis / TGF-β-related genes

rep_gene_table <- rep_gene_table %>%
  mutate(
    Group = case_when(
      SYMBOL %in% c("Col1a1", "Col1a2", "Col3a1", "Fn1") ~ "ECM / collagen",
      SYMBOL %in% c("Lox", "Postn") ~ "Fibrosis / repair",
      SYMBOL %in% c("Mmp2", "Mmp9", "Timp1") ~ "Matrix remodeling",
      SYMBOL %in% c("Tgfb1", "Tgfbr1", "Tgfbr2", "Smad2", "Smad3", "Serpine1") ~ "TGF-beta signaling",
      TRUE ~ "Other"
    )
  )


plot_df <- rep_gene_table %>%
  pivot_longer(
    cols = matches("_(log2FC|pvalue|padj)$"),
    names_to = c("Comparison", ".value"),
    names_pattern = "(.*)_(log2FC|pvalue|padj)"
  ) %>%
  mutate(
    sig_label = case_when(
      !is.na(padj) & padj < 0.05 ~ "**",
      !is.na(pvalue) & pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    SYMBOL = factor(SYMBOL, levels = rev(genes_interest)),
    Group = factor(
      Group,
      levels = c("ECM / collagen", "Fibrosis / repair", "Matrix remodeling", "TGF-beta signaling")
    )
  )

ggplot(plot_df, aes(x = Comparison, y = SYMBOL, fill = log2FC)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sig_label), size = 6) +
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "log2FC"
  ) +
  labs(
    title = "Representative ECM / fibrosis / TGF-beta-related genes",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0, face = "bold"),
    strip.placement = "outside"
  )


ggsave(
  "Representative_ECM_fibrosis_TGFb_genes_heatmap.png",
  width = 10,
  height = 7,
  dpi = 300
)
