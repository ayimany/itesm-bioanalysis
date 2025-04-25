#Ejemplo de lectura de archivo fasta con Biostrings

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Biostrings")

library(Biostrings)

seqs <- readDNAStringSet("covid19_mex.fasta")
seqs
seqs <- as.character(seqs)
seqs