if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

if (!require("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!require("magrittr", quietly = TRUE))
  install.packages("magrittr")

if (!require("readr", quietly = TRUE))
  install.packages("readr")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")

# Install the Human package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

# Attach the library
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

expression_df <- readr::read_tsv("GSE183516_Normalized_counts_AT2_pAT2_samples.txt")

# Tuck away the Gene ID column as row names
expression_df = tibble::column_to_rownames(expression_df,var="...1")

expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

mapped_list <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "list"
)



