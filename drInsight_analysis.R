#!/drInsight_analysis.R

#The script takes the result of DEGs identified by DESeq2 tool and find drugs which can revert the expression showcased by the input.

#The present script has the following input data:
#1) Path to the directory where DESeq2 output files reside saved in the dirPath object.
#The DESeq2 file has the following columns: GeneID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj.
#2) Specify a pattern in order to select the right files from the directory saved in the "patternFiles" object.
#3) Path to the rank matrix file obtained from the Broad Institute saved in the "Drug_rankMatrix" object.

#The present script has the following output files:
#1) A DrInsight result file table for each DESeq2 file found (format: drInsight_inputFileName_res.csv)
#The output file has the following columns: drug, pval, FDR, pval_pass, FDR_pass.

#Load Library
library(DrInsight)  #Drug Repurposing: Integration and Systematic Investigation of Genomic High Throughput Data
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(gprofiler2) #Gene list functional enrichment analysis and namespace conversion

rm(list=ls(all=TRUE)) #Remove all object saved in the system

###### USER's Input Data
dirPath <- "4_DEG_Analysis/Hsapiens"   #Path to the directory where the DESeq2 output files are saved.
patternFiles <- "DESeq2" #Specify the pattern in order to select all the input data
Drug_rankMatrix <- "rankMatrix.txt" #Drug CMap downloaded from the Broad Institute
#######################

###### Upload the rankMatrix from the Broad Institute
cmap.ref.profiles = get.cmap.ref(cmap.data.path = Drug_rankMatrix, probe.to.genes = probe.to.genes, drug.info = drug.info)

###### UPLOAD INPUT FILES
list_of_files <- list()
files <- list.files(path=dirPath, pattern=patternFiles)

if (dirPath==".")
{
  for (i in files) 
  {list_of_files[[i]] <- read.csv2 (i, header=TRUE, sep = ",")}
} else {
  for (i in files) 
  {list_of_files[[i]] <- read.csv2 (paste(dirPath, "/", i, sep = ""), header=FALSE, sep = ",")}
}

##########
for (i in 1:length(list_of_files))
{
  nowInput <- list_of_files[[i]] [,c(1,3)]
  colnames(nowInput) <- c("ID", "score")
  
  #Get Table with gene symbol and score (Log2FC)
  idConvert <- gconvert(nowInput[,1], organism = "hsapiens", target = "ENSG", numeric_ns = "", mthreshold = Inf, filter_na = TRUE) %>% arrange(input)
  nowInput <- nowInput %>% filter(ID %in% idConvert$input) %>% arrange(ID)
  nowInput <- cbind(nowInput, idConvert) %>% select(name, score)
  colnames(nowInput) <- c("geneSymbol", "score")
  
  # DrInsight Drug Identification
  drug.ident.res <- drug.ident(query.data=nowInput, cmap.ref.profiles=cmap.ref.profiles, CEG.threshold=0.05, repurposing.unit="treatment", connectivity="negative")
  drug.pvals <- drug.ident.res$drug.pvals %>% filter(pval < 0.05)
  drug.pvals$pval_pass <- 1
  drug.pvals$FDR_pass <- ifelse(drug.pvals$FDR < 0.05, 1, 0)
  
  #Save result
  write.csv(drug.pvals, file = paste("5_Drug_DrInsight_Analysis/Hsapiens/","drInsight_", names(list_of_files)[i], "_res.csv", sep = ""))
}
