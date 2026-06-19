# Metascape BZ 1 week p-value 

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

meta_enrich <- read_excel("metascape_BZ 1week_p-value.xlsx", sheet = "Enrichment")

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  slice_head(n = 20)

if (median(meta_top20$`Log(q-value)`, na.rm = TRUE) < 0) {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = -`Log(q-value)`)
} else {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = `Log(q-value)`)
}

meta_top20 <- meta_top20 %>%
  arrange(desc(minus_log10_q)) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )

highlight_terms <- c(
  "extracellular matrix organization",
  "ECM-receptor interaction",
  "tissue morphogenesis",
  "heart development"
)


meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  mutate(
    minus_log10_q = -`Log(q-value)`,
    Description_clean = tolower(Description)
  ) %>%
  arrange(desc(minus_log10_q)) %>%
  distinct(Description_clean, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )



meta_top20 <- meta_top20 %>%
  mutate(
    highlight = ifelse(Description %in% highlight_terms, "highlight", "normal")
  )

y_min <- 0.5
y_max <- nrow(meta_top20) + 0.5
x_max <- ceiling(max(meta_top20$minus_log10_q, na.rm = TRUE)) + 1

p_bar <- ggplot(meta_top20) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = minus_log10_q,
      ymin = y_pos - 0.35,
      ymax = y_pos + 0.35,
      fill = minus_log10_q
    ),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_gradient(
    low = "#FFD966",
    high = "#D95F02"
  ) +
  scale_x_continuous(
    breaks = seq(0, x_max, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  labs(
    title = "BZ 1wk Top enriched terms",
    x = expression(-log[10]("q-value")),
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 26),
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_text <- ggplot(meta_top20) +
  geom_text(
    aes(x = 0, y = y_pos, label = Description_display, color = highlight),
    hjust = 0,
    size = 5.4
  ) +
  scale_color_manual(
    values = c(
      "highlight" = "red",
      "normal" = "black"
    )
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(54, 5.5, 5.5, 5.5)
  )

p_final <- p_bar + p_text + plot_layout(widths = c(2.2, 1.45))

p_final

ggsave(
  "Metascape_BZ_1wk_qvalue_paper_style_fixed_horizontal.png",
  plot = p_final,
  width = 15,
  height = 6,
  dpi = 300
)

ggsave(
  "Metascape_BZ_1wk_qvalue_paper_style_fixed_horizontal.pdf",
  plot = p_final,
  width = 15,
  height = 6
)





library(readxl)
library(dplyr)

meta_enrich <- read_excel("metascape_BZ 1week_p-value.xlsx", sheet = "Enrichment")

ecm_collagen_terms <- meta_enrich %>%
  filter(
    grepl(
      "collagen|fiber|fibril|matrix|ECM|matrisome|adhesion",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  select(
    GroupID,
    Category,
    Description,
    `LogP`,
    `Log(q-value)`,
    everything()
  ) %>%
  arrange(`Log(q-value)`)

View(ecm_collagen_terms)




# Metascape IZ 1 week p-value 
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

meta_enrich <- read_excel("metascape_IZ 1week_p-value.xlsx", sheet = "Enrichment")

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  slice_head(n = 20)

if (median(meta_top20$`Log(q-value)`, na.rm = TRUE) < 0) {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = -`Log(q-value)`)
} else {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = `Log(q-value)`)
}

meta_top20 <- meta_top20 %>%
  arrange(desc(minus_log10_q)) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )

highlight_terms <- c(
  "regulation of response to wounding",
  "positive regulation of cell migration",
  "cytokine-mediated signaling pathway",
  "positive regulation of macrophage activation",
  "regulation of phagocytosis"
)


meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  mutate(
    minus_log10_q = -`Log(q-value)`,
    Description_clean = tolower(Description)
  ) %>%
  arrange(desc(minus_log10_q)) %>%
  distinct(Description_clean, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )



meta_top20 <- meta_top20 %>%
  mutate(
    highlight = ifelse(Description %in% highlight_terms, "highlight", "normal")
  )

y_min <- 0.5
y_max <- nrow(meta_top20) + 0.5
x_max <- ceiling(max(meta_top20$minus_log10_q, na.rm = TRUE)) + 1

p_bar <- ggplot(meta_top20) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = minus_log10_q,
      ymin = y_pos - 0.35,
      ymax = y_pos + 0.35,
      fill = minus_log10_q
    ),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_gradient(
    low = "#FFD966",
    high = "#D95F02"
  ) +
  scale_x_continuous(
    breaks = seq(0, x_max, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  labs(
    title = "IZ 1wk Top enriched terms",
    x = expression(-log[10]("q-value")),
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 26),
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_text <- ggplot(meta_top20) +
  geom_text(
    aes(x = 0, y = y_pos, label = Description_display, color = highlight),
    hjust = 0,
    size = 5.4
  ) +
  scale_color_manual(
    values = c(
      "highlight" = "red",
      "normal" = "black"
    )
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(54, 5.5, 5.5, 5.5)
  )

p_final <- p_bar + p_text + plot_layout(widths = c(2.2, 1.45))

p_final

ggsave(
  "Metascape_IZ_1wk_qvalue_paper_style_fixed_horizontal.png",
  plot = p_final,
  width = 15,
  height = 6,
  dpi = 300
)

ggsave(
  "Metascape_IZ_1wk_qvalue_paper_style_fixed_horizontal.pdf",
  plot = p_final,
  width = 15,
  height = 6
)



library(readxl)
library(dplyr)

meta_enrich <- read_excel("metascape_IZ 1week_p-value.xlsx", sheet = "Enrichment")

ecm_collagen_terms <- meta_enrich %>%
  filter(
    grepl(
      "collagen|fiber|fibril|matrix|ECM|matrisome|adhesion",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  select(
    GroupID,
    Category,
    Description,
    `LogP`,
    `Log(q-value)`,
    everything()
  ) %>%
  arrange(`Log(q-value)`)

View(ecm_collagen_terms)




# Metascape IZ Day3 p-value 
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

meta_enrich <- read_excel("metascape_IZ Day3_p-value.xlsx", sheet = "Enrichment")

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  slice_head(n = 20)

if (median(meta_top20$`Log(q-value)`, na.rm = TRUE) < 0) {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = -`Log(q-value)`)
} else {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = `Log(q-value)`)
}

meta_top20 <- meta_top20 %>%
  arrange(desc(minus_log10_q)) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )

highlight_terms <- c(
  "positive regulation of cytokine production",
  "TNF signaling pathway",
  "angiogenesis",
  "extracellular matrix organization"
)


meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  mutate(
    minus_log10_q = -`Log(q-value)`,
    Description_clean = tolower(Description)
  ) %>%
  arrange(desc(minus_log10_q)) %>%
  distinct(Description_clean, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )



meta_top20 <- meta_top20 %>%
  mutate(
    highlight = ifelse(Description %in% highlight_terms, "highlight", "normal")
  )

y_min <- 0.5
y_max <- nrow(meta_top20) + 0.5
x_max <- ceiling(max(meta_top20$minus_log10_q, na.rm = TRUE)) + 1

p_bar <- ggplot(meta_top20) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = minus_log10_q,
      ymin = y_pos - 0.35,
      ymax = y_pos + 0.35,
      fill = minus_log10_q
    ),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_gradient(
    low = "#FFD966",
    high = "#D95F02"
  ) +
  scale_x_continuous(
    breaks = seq(0, x_max, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  labs(
    title = "IZ Day3 Top enriched terms",
    x = expression(-log[10]("q-value")),
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 26),
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_text <- ggplot(meta_top20) +
  geom_text(
    aes(x = 0, y = y_pos, label = Description_display, color = highlight),
    hjust = 0,
    size = 5.4
  ) +
  scale_color_manual(
    values = c(
      "highlight" = "red",
      "normal" = "black"
    )
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(54, 5.5, 5.5, 5.5)
  )

p_final <- p_bar + p_text + plot_layout(widths = c(2.2, 1.45))

p_final

ggsave(
  "Metascape_IZ_Day3_qvalue_paper_style_fixed_horizontal.png",
  plot = p_final,
  width = 15,
  height = 6,
  dpi = 300
)

ggsave(
  "Metascape_IZ_Day3_qvalue_paper_style_fixed_horizontal.pdf",
  plot = p_final,
  width = 15,
  height = 6
)



library(readxl)
library(dplyr)

meta_enrich <- read_excel("metascape_IZ Day3_p-value.xlsx", sheet = "Enrichment")

ecm_collagen_terms <- meta_enrich %>%
  filter(
    grepl(
      "collagen|fiber|fibril|matrix|ECM|matrisome|adhesion",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  select(
    GroupID,
    Category,
    Description,
    `LogP`,
    `Log(q-value)`,
    everything()
  ) %>%
  arrange(`Log(q-value)`)

View(ecm_collagen_terms)











# Metascape IZ 4wks p-value 
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

meta_enrich <- read_excel("metascape_IZ 4wks_p-value.xlsx", sheet = "Enrichment")

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  slice_head(n = 20)

if (median(meta_top20$`Log(q-value)`, na.rm = TRUE) < 0) {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = -`Log(q-value)`)
} else {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = `Log(q-value)`)
}

meta_top20 <- meta_top20 %>%
  arrange(desc(minus_log10_q)) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )

highlight_terms <- c(
  "neutrophil chemotaxis",
  "regulation of monocyte chemotactic protein-1 production"
)


meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  mutate(
    minus_log10_q = -`Log(q-value)`,
    Description_clean = tolower(Description)
  ) %>%
  arrange(desc(minus_log10_q)) %>%
  distinct(Description_clean, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )



meta_top20 <- meta_top20 %>%
  mutate(
    highlight = ifelse(Description %in% highlight_terms, "highlight", "normal")
  )

y_min <- 0.5
y_max <- nrow(meta_top20) + 0.5
x_max <- ceiling(max(meta_top20$minus_log10_q, na.rm = TRUE)) + 1

p_bar <- ggplot(meta_top20) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = minus_log10_q,
      ymin = y_pos - 0.35,
      ymax = y_pos + 0.35,
      fill = minus_log10_q
    ),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_gradient(
    low = "#FFD966",
    high = "#D95F02"
  ) +
  scale_x_continuous(
    breaks = seq(0, x_max, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  labs(
    title = "IZ 4wks Top enriched terms",
    x = expression(-log[10]("q-value")),
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 26),
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_text <- ggplot(meta_top20) +
  geom_text(
    aes(x = 0, y = y_pos, label = Description_display, color = highlight),
    hjust = 0,
    size = 5.4
  ) +
  scale_color_manual(
    values = c(
      "highlight" = "red",
      "normal" = "black"
    )
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(54, 5.5, 5.5, 5.5)
  )

p_final <- p_bar + p_text + plot_layout(widths = c(2.2, 1.45))

p_final

ggsave(
  "Metascape_IZ_4wks_qvalue_paper_style_fixed_horizontal.png",
  plot = p_final,
  width = 15,
  height = 6,
  dpi = 300
)

ggsave(
  "Metascape_IZ_4wks_qvalue_paper_style_fixed_horizontal.pdf",
  plot = p_final,
  width = 15,
  height = 6
)



library(readxl)
library(dplyr)

meta_enrich <- read_excel("metascape_IZ 4wks_p-value.xlsx", sheet = "Enrichment")

ecm_collagen_terms <- meta_enrich %>%
  filter(
    grepl(
      "collagen|fiber|fibril|matrix|ECM|matrisome|adhesion",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  select(
    GroupID,
    Category,
    Description,
    `LogP`,
    `Log(q-value)`,
    everything()
  ) %>%
  arrange(`Log(q-value)`)

View(ecm_collagen_terms)











# Metascape IZ all-time-point p-value 
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

meta_enrich <- read_excel("metascape_IZ all-time-point_p-value.xlsx", sheet = "Enrichment")

capitalize_first <- function(x) {
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}

meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  slice_head(n = 20)

if (median(meta_top20$`Log(q-value)`, na.rm = TRUE) < 0) {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = -`Log(q-value)`)
} else {
  meta_top20 <- meta_top20 %>%
    mutate(minus_log10_q = `Log(q-value)`)
}

meta_top20 <- meta_top20 %>%
  arrange(desc(minus_log10_q)) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )

highlight_terms <- c(
  "regulation of response to wounding",
  "positive regulation of cell migration",
  "cytokine-mediated signaling pathway",
  "positive regulation of macrophage activation",
  "regulation of phagocytosis"
)


meta_top20 <- meta_enrich %>%
  filter(grepl("Summary", GroupID)) %>%
  mutate(
    minus_log10_q = -`Log(q-value)`,
    Description_clean = tolower(Description)
  ) %>%
  arrange(desc(minus_log10_q)) %>%
  distinct(Description_clean, .keep_all = TRUE) %>%
  slice_head(n = 20) %>%
  mutate(
    y_pos = rev(seq_len(n())),
    Description_display = capitalize_first(Description)
  )



meta_top20 <- meta_top20 %>%
  mutate(
    highlight = ifelse(Description %in% highlight_terms, "highlight", "normal")
  )

y_min <- 0.5
y_max <- nrow(meta_top20) + 0.5
x_max <- ceiling(max(meta_top20$minus_log10_q, na.rm = TRUE)) + 1

p_bar <- ggplot(meta_top20) +
  geom_rect(
    aes(
      xmin = 0,
      xmax = minus_log10_q,
      ymin = y_pos - 0.35,
      ymax = y_pos + 0.35,
      fill = minus_log10_q
    ),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_gradient(
    low = "#FFD966",
    high = "#D95F02"
  ) +
  scale_x_continuous(
    breaks = seq(0, x_max, by = 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  labs(
    title = "IZ all-time-point Top enriched terms",
    x = expression(-log[10]("q-value")),
    y = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 26),
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_text <- ggplot(meta_top20) +
  geom_text(
    aes(x = 0, y = y_pos, label = Description_display, color = highlight),
    hjust = 0,
    size = 5.4
  ) +
  scale_color_manual(
    values = c(
      "highlight" = "red",
      "normal" = "black"
    )
  ) +
  scale_x_continuous(
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(y_min, y_max),
    clip = "off"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(54, 5.5, 5.5, 5.5)
  )

p_final <- p_bar + p_text + plot_layout(widths = c(2.2, 1.45))

p_final

ggsave(
  "Metascape_IZ_all-time-point_qvalue_paper_style_fixed_horizontal.png",
  plot = p_final,
  width = 15,
  height = 6,
  dpi = 300
)

ggsave(
  "Metascape_IZ_all-time-point_qvalue_paper_style_fixed_horizontal.pdf",
  plot = p_final,
  width = 15,
  height = 6
)



library(readxl)
library(dplyr)

meta_enrich <- read_excel("metascape_IZ all-time-point_p-value.xlsx", sheet = "Enrichment")

ecm_collagen_terms <- meta_enrich %>%
  filter(
    grepl(
      "collagen|fiber|fibril|matrix|ECM|matrisome|adhesion",
      Description,
      ignore.case = TRUE
    )
  ) %>%
  select(
    GroupID,
    Category,
    Description,
    `LogP`,
    `Log(q-value)`,
    everything()
  ) %>%
  arrange(`Log(q-value)`)

View(ecm_collagen_terms)

