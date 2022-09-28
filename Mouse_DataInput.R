#### Mouse Gene ####
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


# Install the Mouse package
if (!("org.Mm.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}

# Attach the library
library(org.Mm.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

expression_df <- readr::read_tsv("GSE183548_Normalized_counts_tumour_samples.txt")

# Tuck away the Gene ID column as row names
expression_df = tibble::column_to_rownames(expression_df,var="...1")

expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")


expression_df_data <- expression_df
for( i in 1:nrow(expression_df_data)){
  expression_df_data[i, 1] <- expression_df[i, 3] 
  expression_df_data[i, 2] <- var(as.vector(t(expression_df_data[i,4:ncol(expression_df)])))
  }

geneVariation <- expression_df_data[,1:2]
print(geneVariation$Row.names)
print(mean(geneVariation$Row.names))
print(typeof(geneVariation$Row.names))
geneVariation[,2] <- as.numeric(unlist(geneVariation$Row.names))
genevariation_avg <- mean(geneVariation$Row.names) 

library(ggplot2)
ggplot(geneVariation, aes(x = Gene, y = Row.names)) + 
  geom_point()+
  scale_y_continuous(trans = 'log2') +
  ylab("Average Sample Variance")
