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

# Github
if(!require(devtools)) install.packages("devtools")
devtools::install_github("sinhrks/ggfortify")

# Install the Human package
if (!("org.Hs.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}

# Attach the library
library(org.Hs.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

human_cancer <- readr::read_tsv("GSE183516_Normalized_counts_H23_H358_A549_cell_lines.txt")
# Tuck away the Gene ID column as row names
human_cancer = tibble::column_to_rownames(human_cancer,var="...1")

human_cancer <- human_cancer %>%
  tibble::rownames_to_column("Gene")

mapped_list <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = human_cancer$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "filter"
)

# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Entrez") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Entrez)

multi_mapped <- mapped_df %>%
  # Let's count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(Ensembl, name = "entrez_id_count") %>%
  # Arrange by the genes with the highest number of Entrez IDs mapped
  dplyr::arrange(desc(entrez_id_count))

# Let's look at the first 6 rows of our `multi_mapped` object
head(multi_mapped)

human_cancer_data <- human_cancer
for( i in 1:nrow(human_cancer_data)){
  human_cancer_data[i, 1] <- mapped_df[i, 2] 
}
view(human_cancer_data)

human_cancer_data_var <- human_cancer_data
for( i in 1:nrow(human_cancer_data_var)){
  human_cancer_data_var[i, 2] <- var(as.vector(t(human_cancer_data[i,2:ncol(human_cancer_data_var)])))
  #print(as.vector(t(expression_df_data[i,2:ncol(df_data_var)])))
}
human_geneVariation <- human_cancer_data[,1:2]
#print(as.vector(t(expression_df_data[1,2:ncol(df_data_var)])))
view(human_geneVariation)
library(ggplot2)
p_human <-ggplot(human_geneVariation, aes(x = Gene, y = MUG207A1)) + 
  geom_point()+
  scale_y_continuous(trans = 'log2') +
  ylab("Average Sample Variance")
slotNames("Matrix")
p_human + scale_y_continuous(trans = 'log2')
view(p_human)



#MOUSE
# Install the Mouse package
if (!("org.Mm.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}

# Attach the library
library(org.Mm.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

mouse_cancer <- readr::read_tsv("GSE183548_Normalized_counts_3LL_samples.txt")

# Tuck away the Gene ID column as row names
mouse_cancer = tibble::column_to_rownames(mouse_cancer,var="...1")

mouse_cancer <- mouse_cancer %>%
  tibble::rownames_to_column("Gene")


mouse_cancer_data <- mouse_cancer
for( i in 1:nrow(mouse_cancer_data)){
  mouse_cancer_data[i, 1] <- mouse_cancer[i, 3] 
  mouse_cancer_data[i, 2] <- var(as.vector(t(mouse_cancer_data[i,4:ncol(mouse_cancer)])))
}
view(mouse_cancer_data)
mouse_geneVariation <- mouse_cancer_data[,1:2]
mouse_geneVariation[,2] <- as.numeric(unlist(mouse_geneVariation$Row.names))
mouse_genevariation_avg <- mean(mouse_geneVariation$Row.names) 

library(ggplot2)
p_mouse <- ggplot(mouse_geneVariation, aes(x = Gene, y = Row.names)) + 
  geom_point() +
  scale_y_continuous(trans = 'log2') +
  ylab("Average Sample Variance")
slotNames("Matrix")

# Log base 10 scale + log ticks (on left and bottom side)
p_mouse + scale_y_continuous(trans = 'log2')
view(p_mouse)
