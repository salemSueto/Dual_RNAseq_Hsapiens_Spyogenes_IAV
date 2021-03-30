#!/count_to_tpm.R

#The script takes a list of count tables and calculate the TPM values 
#based on the gene's length and the mean of the fragment length sequenced.

#The script has the following input:
#1) Path for the selected count table files saved in the "pathFiles" vector object.
#The count tables need to have the first column as ID followed by the samples. 
#2) Path to a text file which shows the gene's ID and its length in two columns called: ID and Length.
#3) The mean of the Fragment Length sequenced for the selected samples.

#The script outputs one file for each count table in the analysis called "star_htseq_rev_tpm.csv" 
#and it is saved in the directory the input count table is present. 
#The output has the following columns: ID, sample1, sampleN, Length, sample1_tpm, sampleN_tpm, TPM_Average.

#Load Library
library (tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.

rm(list = ls()) #Remove all saved objects

########### User's Input ###########

#Path to the count table files
pathFiles <- c("3_Quantification/Spyogenes_AP1/htseq_all_sample_strand_rev_count.txt", 
               "3_Quantification/Spyogenes_NZ131/htseq_all_sample_strand_rev_count.txt",
               "3_Quantification/IAV/htseq_all_sample_strand_rev_count.txt")

#The file needs to have atleast two columns called "ID" (list of gene/transcript ID) and "Length" (length of the gene/transcript).
geneLength <- "geneLength.txt"

# Mean of the Fragment Length sequenced
mean_read_seq <- 75  
################################################


###### Get the IDs' length file
ID_length <- read.table(geneLength, header = TRUE, sep = "\t")

###### TPM calculation function (https://gist.github.com/slowkow/c6ab0348747f86e2748b)
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))
  
  # Exclude genes with length less than the mean fragment length (NOT IN THIS CASE)
  idx <- apply(effLen, 1, function(x) min(x) > 1) 
  
  ######
  idx_true <- ID_length$ID [which(idx == "TRUE")]   #IDs where it is true
  idx_false <- ID_length$ID [which(idx == "FALSE")] #IDs where it is false
  ######
  
  counts <- counts[idx,]
  effLen <- effLen[idx,]
  featureLength <- featureLength[idx]
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  #Prepare TPM output
  colnames(tpm) <- colnames(counts) # Copy the column names from the original matrix.
  rownames(tpm) <- idx_true
  
  idx_false_df <- matrix(data = 0, nrow = length(idx_false), ncol = ncol(counts)) 
  colnames(idx_false_df) <- colnames(counts)
  rownames(idx_false_df) <- idx_false
  
  final_tpm <- as.data.frame(rbind(tpm, idx_false_df)) %>% rownames_to_column(var = "ID") %>% arrange_at(1) %>% column_to_rownames("ID")
  return(final_tpm)
}


###### Calculate TPM values for each count table
for (iFile in pathFiles) 
{
  iFile_name <- str_split(iFile, pattern = "/")[[1]]
  iCT <- read.table (iFile, header=TRUE, sep = "") %>% arrange(ID) %>% column_to_rownames(var = "ID")
  filterInfo <- ID_length %>% filter(ID %in% rownames(iCT)) %>% arrange(ID)
  
  #Calculate TPM from the count table
  meanFragmentLength <- rep(mean_read_seq, ncol(iCT))
  iCT_tpm <- counts_to_tpm (iCT, filterInfo$Length, meanFragmentLength)
  colnames(iCT_tpm) <- paste (colnames(iCT_tpm), "_tpm", sep = "")
  
  iCT$Length <- filterInfo$Length
  iCT_final <- cbind(iCT, iCT_tpm)
  
  tpm_col_index <- grep(pattern = "tpm", colnames(iCT_final))
  iCT_final$TPM_Average <- rowSums(iCT_final[tpm_col_index])/length(tpm_col_index)
  iCT_final <- iCT_final %>% arrange(colnames(iCT_final)[ncol(iCT_final)])
  iCT_final <- iCT_final %>% arrange(desc(TPM_Average))
  
  write.csv (iCT_final, file = paste("6_TPM_Analysis", "/", iFile_name[2], "/", "star_htseq_rev_tpm.csv", sep = ""))   #Save the TPM result
}
