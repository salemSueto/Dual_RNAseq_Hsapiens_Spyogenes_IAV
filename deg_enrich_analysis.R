#!/deg_enrich_analysis.R
#The script takes several inputs from the user.
#The path to the count table files saved in the vector object called "path_RawCount".
#The raw table files need to have the first column as ID followed by the samples with their gene count.
#The "path_samplesMetadata" object takes the path to the metadata about the samples.
#The samples's metadata have two columns: ID column with the sample name and Group column with the sample's group membership.
#The "groupComparisons" object vector has the list of group comparisosns to perform saved as Group1 - Group2.
#The DESeq2 input section has the value to identify DEGs. 
#The "cutOff" remove genes which row sum is less than its value. The "padj_FDR_value" is the adjusted p-value.
#The "gProfiler_supported_org" need the scientific name of the organism to perform enrichment analysis.
#The "REVIGO's Input format" has all the inputs need to perform at the REVIGO website.

#The following output files are created:
#1) DESeq2 output file (only DEGs) for each group comparison (format: Group1_vs_Group2_DESeq2_cutOff_adjPvalue.csv).
#2) One single file for the enrichment analysis for all group comparisons (format: OrganismName_gProfiler_enrich.csv).
#3) REVIGO results for each GO class (format: OrganismName_revigo_GO/class.csv).
#4) Tables obtained after the cluster generation from REVIGO's GO semantic space (format: OrganismName_revigoCluster_GO/class.csv) 
#5) Table used to generate the heatmap for the three GO class (format: OrganismName_revigoClusterHeatmap_GO/class.csv)
#6) Table used to generate the heatmap for KEGG (format: OrganismName_heatmapKEGG.csv)
#7) A PDF file with heatmap plots (format: OrganismName_heatmapPlots.pdf).

#Load Library
library(tidyverse)  #Collection of R packages designed for data science: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats.
library(DESeq2) #Differential gene expression analysis based on the negative binomial distribution
library(gprofiler2) #Gene list functional enrichment analysis and namespace conversion
library(rvest)  #Wrappers around the 'xml2' and 'httr' packages to make it easy to download, then manipulate, HTML and XML.
library(tidytext)   #Text mining tasks, and plot generation are easier and more effective with tools already in wide use: dplyr, broom, tidyr, and ggplot2.
library(gridExtra)  #Provides the ability to arrange multiple grid-based plots on a page, and draw tables.
library(stringi)    #A multitude of language processing tools: pattern searching, concatenation, sorting, padding, wrapping, and many more.
library(scales) #Graphical scales map data to aesthetics, and provide methods for automatically determining breaks and labels for axes and legends.

rm(list=ls(all=TRUE))

##### User's INPUT #####
path_RawCount <- c("3_Quantification/Hsapiens/htseq_all_sample_strand_rev_count.txt")   #Path to the count table files
path_samplesMetadata <- "samples_metadata.txt"                                          #Path to the metadata file (#ID #Group)
groupComparisons <- c("Inf_AP1 - Control", "Inf_NZ131 - Control", "Inf_IAV - Control")  #Specify the group comparisons as: Group1 - Group2

##### DESeq2 input
cutOff <- 0
padj_FDR_value <- 0.0005

##### Enrichment Analysis -> Supported Organims for g:Profiler tool (e.i. "Homo sapiens") -> (https://biit.cs.ut.ee/gprofiler/page/organism-list)
gProfiler_supported_org <- c("Homo sapiens") #Scientific Name 

# REVIGO's Input format
baseurl <- "http://revigo.irb.hr/"  #REVIGO website link (Do not change)
cutoff <- "0.40" #Allowed values: "0.90" "0.70" "0.50" "0.40" 
isPValue <- "yes" #Allowed values: "yes"  "no"
whatIsBetter <- "higher" #Allowed values: "higher" "lower" "absolute" "abs_log"
goSizes <- "0" #("whole UniProt")
measure <- "SIMREL" #Allowed values: "RESNIK" "LIN" "SIMREL" "JIANG"
######


###### Loop through each count table to identify the DEGs ######
for (iCT in path_RawCount) 
{
  deg_list <- list() #List of the DEGs identified
  
  # Get the count table & samples' metadata files
  rawCount_df <- read.table(iCT, header=TRUE)
  samplesMetadata_df <- read.table(path_samplesMetadata, header=TRUE)
  
  # Loop through the different group comparisons
  for (iGC in groupComparisons) 
  {
    #Check that there are two saved groups
    if (length (str_split(iGC, pattern = " - ")[[1]]) == 2)
    {
      groups <- str_split(iGC, pattern = " - ")[[1]]
      
      #Group1
      gr1ColNames <- samplesMetadata_df %>% filter(Group==groups[1])
      group1 <- grep(paste(gr1ColNames$ID, collapse="|"), colnames(rawCount_df))
      group1_df <- rawCount_df[, group1]
      
      #Group2
      gr1ColNames <- samplesMetadata_df %>% filter(Group==groups[2])
      group2 <- grep(paste(gr1ColNames$ID, collapse="|"), colnames(rawCount_df))
      group2_df <- rawCount_df[, group2]
      
      #Replicates
      replicate <- c(ncol(group1_df), ncol(group2_df))
      
      #Check that each group has more than one replicate
      if (all(replicate > 1)) 
      {
        cat("The Differential Expressed Analisis will be started!\n")
        
        #Fix the Count Table -> nowTable
        nowTable <- cbind(group1_df, group2_df)
        nowTable <- nowTable %>% mutate_if (is.character,as.numeric)
        rownames(nowTable) <- rawCount_df[,1]
        
        #Loop through each instance of cutoff value
        for (c in cutOff) 
        {
          print(paste(iGC, " * ", "CutOff: ", c, sep = ""))
          
          colSelect <- append(c(1:replicate[1]), c((replicate[1]+1):(sum(replicate))))
          colSelect <- colnames(nowTable)[colSelect]
          
          compaDF <- nowTable %>% 
            dplyr::select(colSelect) %>% 
            rownames_to_column() %>%
            mutate(Sum=rowSums(dplyr::select(., -starts_with("rowname")))) %>% 
            filter(Sum >= c) %>% 
            dplyr::select(-Sum) %>%
            column_to_rownames()
          
          # Group Factors
          sampleGroups <- append (rep(1, replicate[1]), rep(2, replicate[2]))
          myfactors_D_E <- factor(sampleGroups)             #DESeq2 factors
          
          #------- DESeq2 ------#
          matrixCountTable <- round (as.matrix(compaDF))
          coldata <- data.frame(row.names=colnames(matrixCountTable), myfactors_D_E)
          dds <- DESeqDataSetFromMatrix(countData=matrixCountTable, colData=coldata, design=~myfactors_D_E)
          res <- as.data.frame(results(DESeq(dds)))
          
          #Adjusted P-values
          for (pv in padj_FDR_value) 
          {
            #DESeq2
            deseq2_deg <- res %>% filter(padj < pv)
            #deseq2_deg$GroupComparison <- iGC
            write.csv (deseq2_deg, file = paste("4_DEG_Analysis/Hsapiens/", groups[1], "_vs_", groups[2], "_DESeq2_", c, "_", pv, ".csv", sep = ""))
            deg_list[[paste (groups[1], groups[2], "DESeq2", c, pv, sep = "_")]] <- rownames(deseq2_deg)
          }
        }
      } else {cat(paste(gc, " * ", "Replicates -> ", "The samples' replicates need to be above 1\n"))}
    } else {cat(paste(gc, "-> ", "The present comparison needs exactly two groups for the comparison 1\n", sep = ""))}
  }
}


###### g:Profiler Supported Organism Enrichment Analysis ######
enrichPlot <- list()

for (iSpecies in gProfiler_supported_org) 
{
  ##### Enrichment Analysis -> g:Profiler #####
  gProfiler_Res <- NULL
  gPr_iSpecies <- paste(substr (tolower(iSpecies), 1, 1), str_split (tolower(iSpecies), pattern = " ")[[1]][2], sep = "")
  
  tryCatch (
    {
      gProfiler_Res <- gost(deg_list, organism = gPr_iSpecies, ordered_query = FALSE,
                            multi_query = FALSE, significant = TRUE, exclude_iea = TRUE,
                            measure_underrepresentation = FALSE, evcodes = FALSE,
                            user_threshold = 0.05, correction_method = "gSCS",
                            domain_scope = "annotated", custom_bg = NULL,
                            numeric_ns = "", sources = NULL)
    },
    error = function(e) 
    {
      message(e)
      return(gProfiler_Res <- NULL)
    }
  )
  
  if (is.null(gProfiler_Res) == FALSE) 
  {
    enrich_go_kegg <- gProfiler_Res$result[,c(1,3,9,10,11)] %>% filter(source %in% c("GO:BP", "GO:CC", "GO:MF", "KEGG"))
    write.csv(enrich_go_kegg, file = paste("4_DEG_Analysis/Hsapiens/", gPr_iSpecies, "_gProfiler_enrich.csv", sep = ""))
    
    ##### REVIGO Analysis #####
    revigo_session <- html_session(baseurl)
    revigo_form <- html_form(revigo_session)[[1]] 
    
    # Loop through the GO terms
    for (go in c("GO:BP", "GO:CC", "GO:MF")) 
    {
      # Number of Cluster
      nCluster <- 40
      
      # Gene Ontology
      filterGO <- enrich_go_kegg %>% filter (source==go)
      goList <- paste(unique(filterGO$term_id), collapse="\n")
      
      filled_form <- set_values(revigo_form, 'goList'= goList, 
                                'cutoff'= cutoff, 'isPValue'= isPValue, 
                                'whatIsBetter'= whatIsBetter, 'goSizes'= goSizes, 'measure'= measure)
      
      # Get the Revigo Summary Table & Save it as a CSV file
      result_page <- submit_form(revigo_session, filled_form, submit='startRevigo')
      GORevigo <- result_page %>% follow_link("Export results to text table (CSV)")
      httr::content(GORevigo$response, as="text") %>% write_lines("revigoSummaryOnline.csv") 
      
      # Get the REVIGO summary result from the saved file
      revigoSummary <- read.csv2("revigoSummaryOnline.csv", header=TRUE, sep=",")
      revigoSummary$frequency <- gsub('%', '', revigoSummary$frequency)
      revigoSummary$frequency <- as.numeric( as.character(revigoSummary$frequency) );
      
      #Get only the GO terms represented in the "semantic space"
      uniqueGO <- revigoSummary %>% filter(plot_X != "null")
      iNum <- c(4,5,6,7,8) # Change the column into numeric: plot_X - plot_Y - plot_size - uniqueness - dispensability
      uniqueGO [ , iNum] <- apply(uniqueGO [ , iNum], 2, function(x) as.numeric(as.character(x)))
      
      # Add "HeadGO" Column which shows which GO REVIGO term has been put as the main GO for that group of terms
      revigoSummary$HeadGO <- NA
      headGO <- "NA"
      
      for (i in 1:nrow(revigoSummary)) 
      {
        if(revigoSummary$plot_X[i]!="null") {headGO <- revigoSummary$term_ID[i]}
        revigoSummary$HeadGO[i] <- headGO
      }
      write.csv(revigoSummary, file = sprintf("4_DEG_Analysis/Hsapiens/%s_revigo_%s.csv", gPr_iSpecies, go))   #Save the REVIGO result
      
      
      # GENE ONTOLOGY REVIGO CLUSTER ANALYSIS
      #Find the clusters present in the REVIGO results based on their sematic space
      clusterDF <- uniqueGO %>% remove_rownames %>% column_to_rownames(var="term_ID")
      if (nCluster >= nrow(clusterDF)) {nCluster <- nrow(clusterDF)-1}
      
      set.seed(100)
      clusters <- kmeans(clusterDF[,c(3,4)], nCluster)
      clusterDF$Cluster <- clusters$cluster
      
      # Group elements from the same cluster
      revigoCluster <- data.frame(matrix(ncol=10, nrow=nCluster))
      colnames(revigoCluster) <- c("Cluster", "PlotX", "PlotY", "RevigoGOs", "RevigoRep", "RevigoDescription", "AllGOs", "AllWords", "AllDescription", "FinalDescription")
      
      for (i in sort(unique(clusterDF$Cluster)))
      {
        filterGO_cluster <- clusterDF %>% filter(Cluster==i)
        
        # REVIGO GOs
        revigoCluster$Cluster[i] <- i
        revigoCluster$PlotX[i] <- mean(filterGO_cluster$plot_X)
        revigoCluster$PlotY[i] <- mean(filterGO_cluster$plot_Y)
        revigoCluster$RevigoGOs[i] <- paste (rownames(filterGO_cluster), collapse=" ")
        revigoCluster$RevigoRep[i] <- filterGO_cluster[which.max(filterGO_cluster$representative),][1,1] #Term with the max value of "representative of the cluster
        
        # REVIGO Cluster Description -> Concagenate the HEADs' term description of the cluster - Split - Delete generic terms (i.e. "of", "a", "an")
        wordOcc <- paste(filterGO_cluster$description, collapse=" ") %>% str_split(pattern = " ") 
        wordOcc <- wordOcc[[1]]
        indexWords <- !(wordOcc %in% stop_words$word)
        wordOcc <- as.data.frame(table (wordOcc[indexWords])) %>% arrange(desc (Freq))
        wordOcc <- as.vector(wordOcc$Var1[1:5])
        revigoCluster$RevigoDescription[i] <- paste (wordOcc[!is.na(wordOcc)], collapse = " ")
        
        #All GOs Columns
        filterAllGO <- revigoSummary %>% filter(HeadGO %in% rownames(filterGO_cluster))
        
        revigoCluster$AllGOs[i] <- paste (filterAllGO$term_ID, collapse=" ")
        revigoCluster$AllWords[i] <- paste(filterAllGO$description, collapse=" ")
        
        wordOcc <- paste(filterAllGO$description, collapse=" ") %>% str_split(pattern = " ") 
        wordOcc <- wordOcc[[1]]
        indexWords <- !(wordOcc %in% stop_words$word)
        wordOcc <- as.data.frame(table (wordOcc[indexWords])) %>% arrange(desc (Freq))
        wordOcc <- as.vector(wordOcc$Var1[1:5])
        revigoCluster$AllDescription[i] <- paste (wordOcc[!is.na(wordOcc)], collapse = " ")
        
        #Final Description of the Cluster
        revigoCluster$FinalDescription[i] <- paste(revigoCluster$RevigoRep[i], revigoCluster$AllDescription[i], nrow(filterAllGO), sep = " * ")
      }
      write.csv(revigoCluster, file = sprintf("4_DEG_Analysis/Hsapiens/%s_revigoCluster_%s.csv", gPr_iSpecies, go))   #Save the REVIGO_Cluster result
      
      
      #Get the lowest pvalue of the group comparison's GO terms and the cluster combination
      heatmapGO <- data.frame(matrix(ncol = length (unique(filterGO$query)), nrow = nrow(revigoCluster)))
      colnames(heatmapGO) <- unique(filterGO$query)
      rownames(heatmapGO) <- revigoCluster$FinalDescription
      
      for (n in 1:nrow(revigoCluster)) 
      {
        allGO <- str_split(revigoCluster$AllGOs[n], " ") [[1]] #Get all GOs 
        filterGOGroup <- filterGO %>% filter(term_id %in% allGO)  #Get all the rows which have those GOs
        
        uniqueGroup <- unique (filterGOGroup$query)
        for (k in uniqueGroup) 
        {
          pvaluesList <- filterGOGroup %>% filter(query == k)
          heatmapGO[n, grep(k, colnames(heatmapGO))] <- min(pvaluesList$p_value)
        }
      }
      write.csv(heatmapGO, file = sprintf("4_DEG_Analysis/Hsapiens/%s_revigoClusterHeatmap_%s.csv", gPr_iSpecies, go))   #Save the REVIGO_Cluster_Heatmap result
      
      #Create the Heatmap plot
      GroupName <- c()
      for (i in 1:nrow(heatmapGO)) 
      {
        groupPresence <- colnames(heatmapGO)[sapply(heatmapGO[i,], function(x) any(!is.na(x)))]
        groupPresence <- paste(groupPresence, collapse = " * ")
        GroupName <- append(GroupName, groupPresence)
      }
      
      heatmapGO$GroupName <- GroupName
      heatmapGO$NGroups <- length (unique(filterGO$query)) - rowSums(is.na(heatmapGO))
      heatmapGO <- heatmapGO %>% arrange(NGroups, GroupName) %>% dplyr::select(-c(NGroups, GroupName))
      
      heatmapGO <- tibble::rownames_to_column(heatmapGO, "Description")
      heatmapGO_melt <- reshape2::melt(heatmapGO,id.vars=c("Description"))
      
      heatmapGO_melt$Description <- factor(heatmapGO_melt$Description, levels = heatmapGO$Description)
      
      plotName <- stri_replace_all_fixed (paste(iSpecies, go), pattern = c(" "), replacement = c("_"), vectorize_all = FALSE)
      enrichPlot[[plotName]] <- heatmapGO_melt %>% 
        ggplot(aes(x=variable, y=factor(Description), fill=value)) +
        geom_raster(aes(fill = value)) +
        scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
        geom_text(size=3, aes(label=scientific(value, digits = 3))) +
        theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
        theme(axis.text.y = element_text(size=10, face="bold")) +
        coord_fixed(ratio=ncol(heatmapGO)/nrow(heatmapGO))
    }
    
    ##### KEGG Heatmap creation #####
    filterKEGG <- filter(enrich_go_kegg, source=="KEGG")
    heatmapKEGG <- data.frame(matrix(ncol = length (unique(filterKEGG$query)), nrow = length (unique (filterKEGG$term_id))))
    colnames(heatmapKEGG) <- unique(filterKEGG$query)
    rownames(heatmapKEGG) <- unique (filterKEGG$term_name)
    
    for (i in 1:nrow(filterKEGG)) 
    {heatmapKEGG[filterKEGG$term_name[i], filterKEGG$query[i]] <- filterKEGG$p_value[i]}
    write.csv(heatmapKEGG, file = paste("4_DEG_Analysis/Hsapiens/", gPr_iSpecies, "_heatmapKEGG.csv", sep = ""))   #Save the gProfiler_Heatmap result
    
    #Create the Heatmap plot
    GroupName <- c()
    for (i in 1:nrow(heatmapKEGG)) 
    {
      groupPresence <- colnames(heatmapKEGG)[sapply(heatmapKEGG[i,], function(x) any(!is.na(x)))]
      groupPresence <- paste(groupPresence, collapse = " * ")
      GroupName <- append(GroupName, groupPresence)
    }
    
    heatmapKEGG$GroupName <- GroupName
    heatmapKEGG$NGroups <- length (unique(filterGO$query)) - rowSums(is.na(heatmapKEGG))
    heatmapKEGG <- heatmapKEGG %>% arrange(NGroups, GroupName) %>% dplyr::select(-c(NGroups, GroupName))
    
    heatmapKEGG <- tibble::rownames_to_column(heatmapKEGG, "Description")
    heatmapKEGG_melt <- reshape2::melt(heatmapKEGG,id.vars=c("Description"))
    heatmapKEGG_melt$Description <- factor(heatmapKEGG_melt$Description, levels = heatmapKEGG$Description)
    heatValues <- heatmapKEGG_melt$value [!heatmapKEGG_melt$value %in% c(NA)]
    
    plotName <- stri_replace_all_fixed (paste(iSpecies, "KEGG"), pattern = c(" "), replacement = c("_"), vectorize_all = FALSE)
    enrichPlot[[plotName]] <- heatmapKEGG_melt %>% 
      ggplot(aes(x=variable, y=Description, fill=value)) +
      geom_raster(aes(fill = value)) +
      scale_fill_continuous(high = "ghostwhite", low = "firebrick1", na.value="gray90") +
      geom_text(size=3, aes(label=scientific(value, digits = 3))) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 10, face = "bold")) +
      theme(axis.text.y = element_text(size=10, face="bold")) +
      coord_fixed(ratio=ncol(heatmapGO)/nrow(heatmapGO))
    
    ##### Delete the "revigoSummaryOnline.csv" file
    if (file.exists("revigoSummaryOnline.csv")) {file.remove("revigoSummaryOnline.csv")}
    
    ##### Save Heatmap plots #####
    pdf(paste("4_DEG_Analysis/Hsapiens/", iSpecies, "_heatmapPlots.pdf", sep = ""), onefile=TRUE)
    for (i in seq(length(enrichPlot))) {grid.arrange (enrichPlot[[i]])}
    dev.off()
  }
}
