#!/blast_lowest_evalue_fasta.R

#The script takes the result of the BLAST analysis and it identifies for the target with the lowest evalue.
#Then, it creates a fasta file with the sequence of the query and the target.

#The script has the following input:
#1) Path to the BLAST results (in this case for AP1 and NZ131). The input file has the following columns but it is without #headers: qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, stitle.
#2) Path to the proteome fasta file (in this case for AP1 and NZ131).

#The script output files are:
#1) A text file with the same structure as the BLAST result (for AP1 and NZ131) but it only showcases the target with the #lowest evalue.
#2) The output fasta file has the protein sequence for both query and target for both AP1 and NZ131 BLAST analysis.

#Load Library
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(Biostrings) #Efficient manipulation of biological strings

rm(list = ls()) #Remove all saved objects

###### USER's INPUT ######
ap1_hom_res_path <- "7_Homology_Analysis/Sp_selected_AP1_vs_NZ131_db_blast_res.txt"
ap1_proteome_path <- "Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.pep.all.fa"

nz131_hom_res_path <- "7_Homology_Analysis/Sp_selected_NZ131_vs_AP1_db_blast_res.txt"
nz131_proteome_path <- "Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.pep.all.fa"
##### ####### ####### ######

######Proteome
ap1_proteome_fa <- readDNAStringSet(ap1_proteome_path)
nz131_proteome_fa <- readDNAStringSet(nz131_proteome_path)
ap1_nz131_proteome_fa <- c(ap1_proteome_fa, nz131_proteome_fa)

###### S. pyogenes AP1 Homology Analysis
ap1_hom_res <- read.table(ap1_hom_res_path, header=FALSE, fill = TRUE, sep = "\t")
colnames(ap1_hom_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle") 
ap1_hom_res_final <- ap1_hom_res %>% group_by(qseqid) %>% filter(evalue==min(evalue))
write.csv(ap1_hom_res_final, file = "7_Homology_Analysis/AP1_hom_blast_res_final.csv")

###### S. pyogenes NZ131 Homology Analysis
nz131_hom_res <- read.table(nz131_hom_res_path, header=FALSE, fill = TRUE, sep = "\t")
colnames(nz131_hom_res) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle") 
nz131_hom_res_final <- nz131_hom_res %>% group_by(qseqid) %>% filter(evalue==min(evalue))
write.csv(nz131_hom_res_final, file = "7_Homology_Analysis/NZ131_hom_blast_res_final.csv")

###### Get the sequence for both query and target from both homology analyses
ap1_queryTarget <- c(ap1_hom_res_final$qseqid, ap1_hom_res_final$stitle)
nz131_queryTarget <- c(nz131_hom_res_final$qseqid, nz131_hom_res_final$stitle)
ap1_nz131_queryTarget <- unique (c(ap1_queryTarget, nz131_queryTarget))

selIndex <- c()
for (qT in ap1_nz131_queryTarget) {selIndex <- append(selIndex, grep(qT, names(ap1_nz131_proteome_fa)))}
ap1_nz131_queryTarget_fa <- ap1_nz131_proteome_fa[selIndex]
writeXStringSet(ap1_nz131_queryTarget_fa, "7_Homology_Analysis/ap1_nz131_queryTarget.fa", append=FALSE, compress=FALSE, compression_level=NA, format="fasta")
