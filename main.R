# Libraries
library(ade4)
library(ape)
library(adegenet)
library(Biostrings)
library(ggplot2)
library(tidyr)

# Handmade methods for performing operations on data
decomposeDnaString <- function(str) {
  nstr <- toupper(str)
  char_vec <- strsplit(nstr, "")[[1]]
  
  count   <- 0;
  a_count <- 0;
  c_count <- 0;
  g_count <- 0;
  t_count <- 0;
  n_count <- 0;
  
  for (c in char_vec) {
    count <- count + 1
    if (c == 'A') {
      a_count <- a_count + 1;
    } else if (c == 'C') {
      c_count <- c_count + 1;
    } else if (c == 'G') {
      g_count <- g_count + 1;
    } else if (c == 'T') {
      t_count <- t_count + 1;
    } else {
      n_count <- n_count + 1;
    }
  }
  
  return(list(
    'seq' = nstr,
    'qes' = rev(nstr),
    'len' = count,
    'N_a' = a_count,
    'N_c' = c_count,
    'N_g' = g_count,
    'N_t' = t_count,
    'N_n' = n_count,
    'P_a' = a_count / count,
    'P_c' = c_count / count,
    'P_g' = g_count / count,
    'P_t' = t_count / count,
    'P_n' = n_count / count
));
}


# Obtain the genome sequences from FASTA files
genome_1 = readDNAStringSet('./resources/GCA_011537295.1/GCA_011537295.1_ASM1153729v1_genomic.fna');
genome_2 = readDNAStringSet('./resources/GCA_011537695.1/GCA_011537695.1_ASM1153769v1_genomic.fna');
genome_3 = readDNAStringSet('./resources/GCA_011537705.1/GCA_011537705.1_ASM1153770v1_genomic.fna');
genome_4 = readDNAStringSet('./resources/GCA_011545165.1/GCA_011545165.1_ASM1154516v1_genomic.fna');
genome_5 = readDNAStringSet('./resources/GCA_011545245.1/GCA_011545245.1_ASM1154524v1_genomic.fna');
genome_6 = readDNAStringSet('./resources/GCA_011545275.1/GCA_011545275.1_ASM1154527v1_genomic.fna');
genome_7 = readDNAStringSet('./resources/GCA_011741995.1/GCA_011741995.1_ASM1174199v1_genomic.fna');
genome_8 = readDNAStringSet('./resources/GCA_011742015.1/GCA_011742015.1_ASM1174201v1_genomic.fna');

alpha_variant = readDNAStringSet("./resources/variants/b117.fasta")
beta_variant = readDNAStringSet("./resources/variants/b1351.fasta")
gamma_variant = readDNAStringSet("./resources/variants/p1.fasta")
delta_variant = readDNAStringSet("./resources/variants/b16172.fasta")
omicron_variant = readDNAStringSet("./resources/variants/b11529.fasta")

# WARNING: delta_variant seems not to work

# Store the sequences in a vector
genomes <- list(
  genome_1,
  alpha_variant,
  beta_variant, 
  gamma_variant, 
  omicron_variant
);

# Convert to character strings
genome_strs <- lapply(genomes, as.character);

# Convert to lengths, frequencies and percentages
genome_freq <- lapply(genome_strs, decomposeDnaString);

# Begin the frequency plot
genome_df <- do.call(rbind, lapply(genome_freq, as.data.frame.list))
genome_df$group <- factor(c("Base", "Alpha", "Beta", "Gamma", "Omicron"))
genome_df_pl <- pivot_longer(genome_df, cols = c(N_a, N_c, N_g, N_t), names_to = "property", values_to = "value")

ggplot(genome_df_pl, aes(x = group, y = value, fill = property)) +
  scale_fill_discrete(
    name = "Base",
    labels = c(
      "N_a" = "A Bases",
      "N_c" = "C Bases",
      "N_g" = "G Bases",
      "N_t" = "T Bases"
  )) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Genome", y = "Value", fill = "Property") +
  ggtitle("Grouped Bar Plot of Vector Properties") +
  theme_minimal()

