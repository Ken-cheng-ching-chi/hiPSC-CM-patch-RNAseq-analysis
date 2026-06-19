# 1. Set working directory
getwd()
setwd("C:/Users/asd33/OneDrive/桌面/bioinformatics/Final project/R scripts")

list.files()


counts <- read.delim(
  "GSE318811_raw_counts_All_Samples.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)


# View count matrix as a table in RStudio
View(counts)

# Check column names
colnames(counts)

# Check number of genes and samples
dim(counts)

# Check first few rows and columns
head(counts[, 1:6])


search()

library("DESeq2")


# ============================================================
# 2. Create sample metadata
# ============================================================

sample_names <- colnames(counts)

coldata <- data.frame(
  sample = sample_names,
  row.names = sample_names
)

parts <- strsplit(sample_names, "_")

coldata$region <- sapply(parts, function(x) x[2])
coldata$treatment <- sapply(parts, function(x) x[3])
coldata$time <- sapply(parts, function(x) x[4])
coldata$replicate <- sapply(parts, function(x) x[5])

coldata$region <- factor(coldata$region, levels = c("IZ", "BZ", "RZ"))
coldata$treatment <- factor(coldata$treatment, levels = c("control", "patch"))
coldata$time <- factor(coldata$time, levels = c("Day3", "1wk", "4wks"))

# Check metadata
head(coldata)
View(coldata)

# Check sample distribution
table(coldata$region, coldata$time, coldata$treatment)