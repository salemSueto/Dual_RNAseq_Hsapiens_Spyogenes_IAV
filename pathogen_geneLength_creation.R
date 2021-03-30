#!/pathogen_geneLength_creation.R
#The script takes, as inputs, the annotation files of the three pathogens and creates a text file called "geneLength.txt".
#The geneLength.txt file has the following columns: ID (geneID), Length (gene's length), Reference (name of the annotation), Type (Gene/Transcript).

#Insert the necessary information in the USER's INPUT section:
#annoPath -> vector list of the pathogens' annotation files.

#Load Library
library(tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.

rm(list = ls()) #Remove all saved objects

###### USER's INPUT ######
annoPath <- c("Reference/Spyogenes_AP1/Streptococcus_pyogenes_gca_000993765.ASM99376v1.37.gtf",
             "Reference/Spyogenes_NZ131/Streptococcus_pyogenes_nz131.ASM1812v1.gtf",
             "Reference/IAV/GCF_000865085.1_ViralMultiSegProj15622_genomic.gff")
##########################


geneLength_df <- as.data.frame(matrix(NA, nrow = 0, ncol = 4))
colnames(geneLength_df) <- c("ID", "Length", "Reference", "Type")

for (iAnnoPath in annoPath)
{
  anno <- rtracklayer::import(iAnnoPath) 
  
  if (grepl( "gtf", iAnnoPath, fixed = TRUE) == TRUE)
  {
    gtf_df <- as.data.frame(anno) %>% filter(type=="exon" & gene_id != "") %>% mutate (Length = end - start + 1)
    gtf_df$ID <- gsub("gene:*","", gtf_df$gene_id)
    
    filter_gtf_df <- gtf_df %>% select(ID, Length)
    filter_gtf_df$Reference <- tail(str_split(iAnnoPath, pattern = "/")[[1]], n=1)
    filter_gtf_df$Type <- "Gene"
    geneLength_df <- rbind(geneLength_df, filter_gtf_df)
  }
  
  if (grepl( "gff", iAnnoPath, fixed = TRUE) == TRUE)
  {
    gff_df <- as.data.frame(anno) %>% filter(type=="gene") %>% mutate (Length = end - start + 1)
    
    filter_gff_df <- gff_df %>% select(ID, Length)
    filter_gff_df$Reference <- tail(str_split(iAnnoPath, pattern = "/")[[1]], n=1)
    filter_gff_df$Type <- "Gene"
    geneLength_df <- rbind(geneLength_df, filter_gff_df)
  }
}

write.table (geneLength_df, file="geneLength.txt", sep="\t", row.names=FALSE)
