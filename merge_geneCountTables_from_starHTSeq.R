#!/merge_geneCountTables_from_starHTSeq.R

#The script merges the files generated after the mapping-quantification step with STAR and HTSeq-count with the flag "--quantMode".
#Each input file has the information for one single sample under study and it is divided in two parts.
#The upper section has general information of the quantification result.
#The below section has four columns without headers.
#The first column has the gene's ID. Then, counts for unstranded reads. Then, counts for stranded reads. Finally, counts for reversed stranded reads.
#The script outputs two type of files (count table and summary table) for the unstranded, strand, and reverse strand.
#And, it also outputs one pdf file of barplots for the three quantification types.
#These output files are created for each directory these input files are found.

#Insert the necessary information in the USER's INPUT section:
#pathFiles -> vector list with the directories where the STAR-HTSeq-count files are saved.
#patternFiles -> pattern used to select the right files.

#Load Library
library(tidyverse) #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(gridExtra)  #Provides the ability to arrange multiple grid-based plots on a page, and draw tables.

rm(list = ls()) #Remove all saved objects

###### USER's INPUT ######
pathFiles <- c("3_Quantification/Hsapiens", "3_Quantification/Spyogenes_AP1", "3_Quantification/Spyogenes_NZ131") #List of directories where input files reside
patternFiles <- "ReadsPerGene" #Pattern for the right selection of files
#########################

###### Loop through each Directory to find the count files and merge them
for (iDr in pathFiles) 
{
  #Get filenames
  files <- list.files(path=iDr, pattern=patternFiles)
  
  # Save all files
  list_of_files <- list()
  for (i in files) 
  {list_of_files[[i]] <- read.table(paste(iDr, "/", i, sep = ""), header=FALSE)}
  
  summary_plot_df <- as.data.frame (matrix(nrow = 0, ncol = 4))
  colnames(summary_plot_df) <- c("ID", "variable", "value", "MapType")
  selectColID <- c("unstranded", "strand_yes", "strand_rev")
  
  #Merge the files & separate the count table from the summary
  for (selectCol in c(2,3,4)) #2: unstranded - #3: 1st read strand aligned with RNA (htseq-count -s yes) - #4: 1st read strand aligned with RNA (htseq-count -s reverse)
  {
    # Select the column for each files
    merge_htseq_df <- as.data.frame (matrix(data = NA, nrow = nrow(list_of_files[[1]])-4, ncol = length(list_of_files)))
    summary_htseq_df <- as.data.frame (matrix(data = NA, nrow = 4, ncol = length(list_of_files)))
    
    colnames(merge_htseq_df) <- names(list_of_files)
    colnames(summary_htseq_df) <- names(list_of_files)
    
    rownames(merge_htseq_df) <- list_of_files[[1]]$V1 [-(1:4)]
    rownames(summary_htseq_df) <- list_of_files[[1]]$V1 [1:4]
    
    for (i in 1:length(list_of_files)) 
    {
      sFile <- list_of_files[[i]] %>% select (contains(as.character(selectCol)))
      summary_htseq_df[,i] <- sFile[1:4, 1]
      merge_htseq_df [,i] <- tail(sFile, -4)
    }
    
    merge_htseq_colsums <- colSums(merge_htseq_df)
    merge_htseq_df <- merge_htseq_df %>% rownames_to_column(var = "ID")
    
    summary_htseq_df <- rbind (summary_htseq_df, N_mapped_Quantified = merge_htseq_colsums) %>% rownames_to_column(var = "ID")
    summary_htseq_df$MapType <- selectColID[selectCol-1]
    
    temp_summary <- reshape2::melt(summary_htseq_df,id.vars=c("ID"))
    temp_summary$MapType <- selectColID[selectCol-1]
    summary_plot_df <- rbind(summary_plot_df, temp_summary)
    
    fileName <- paste("htseq_all_sample_", selectColID[selectCol-1], sep = "")
    write.table (merge_htseq_df, paste (iDr, "/", fileName, "_count.txt", sep = ""), sep="\t", row.names=FALSE)
    write.table (summary_htseq_df,paste(iDr, "/", fileName, "_summary.txt", sep = ""), sep="\t", row.names=FALSE)
    
    ####PLOT SUMMARY INFO #####
    summary_plot_df <- summary_plot_df %>% filter(value!=MapType)
    plots_list <- list()
    
    for (i in unique(summary_plot_df$ID))
    {
      plots_list[[paste("Variable:", i)]] <- summary_plot_df %>%
        filter(ID==i) %>%
        ggplot(aes(x=variable, y=value, fill=MapType)) + 
        geom_bar(position = "dodge", stat="identity") +
        labs(title=paste("Variable:", i)) + 
        theme(axis.text.x = element_text(angle=65, hjust=1.0))
    }
    
    pdf(paste(iDr, "/", "htseq_summary_plots.pdf", sep = ""), onefile=TRUE)
    for (i in seq(length(plots_list))) {grid.arrange (plots_list[[i]])}
    dev.off()
  }
}
