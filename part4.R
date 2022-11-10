#Get statistically significant genes from mice
sig_m <- m_deseq_df[which(apply(m_deseq_df, 1, any)), ]

#Get statistically significant genes from human
sig_h <- h_deseq_df[which(apply(h_deseq_df, 1, any)), ]

library(devtools)
install_github("jokergoo/ComplexHeatmap")
