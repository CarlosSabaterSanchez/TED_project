
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to calculate differential expression statistics in pseq objects contaning functional metagenomics data
gene_count <- function(genes_families, genes_metadata, metadata_class_no_formula, selected_levels) {
  
  #Determine all levels for metadata class
  genes_metadata[[metadata_class_no_formula]] <- factor(genes_metadata[[metadata_class_no_formula]], levels=selected_levels)
  class_levels <- levels(genes_metadata[[metadata_class_no_formula]])
  
  #Check only non-redundant genes are present
  #Food and human genes were previously discarded in metagenomic pipeline
  length(which(duplicated(rownames(genes_families)) == TRUE))
  #All genes are non-redundant
  nrow(genes_families) == length(which(duplicated(rownames(genes_families)) == FALSE))
  
  #Create an empty list to store the results
  gene_count_results <- list()
  
  #Prepare gene count
  gene_count <- data.frame(Count=numeric(0), stringsAsFactors = FALSE)
  for (i in 1:nrow(t(genes_families))) {
    sample_genes <- which(t(genes_families)[i,] >= 1, arr.ind = TRUE)
    sample_count <- length(sample_genes)
    gene_count <- rbind(gene_count, sample_count)
  }
  names(gene_count) <- c("Count")
  rownames(gene_count) <- rownames(t(genes_families))
  #There is a R problem with one single column data frames
  gene_count <- cbind(gene_count, gene_count)
  
  #Subset significant class data to get the gene count of each group
  for (i in 1:length(class_levels)) {
    #Select i data
    i_meta <- genes_metadata[grep(class_levels[i], genes_metadata[[metadata_class_no_formula]]), ]
    i_count <- gene_count[rownames(i_meta),]
    
    #Check samples are transposed OK
    rownames(i_meta) == rownames(i_count)
    
    #Bind all i data
    i_data <- cbind(i_meta[[metadata_class_no_formula]], i_count)
                    # i_meta[[metadata_class_participant]], i_count)
    #Format i data
    names(i_data)[1] <- metadata_class_no_formula
    # names(i_data)[2] <- metadata_class_participant
    i_data <- i_data[,-c(ncol(i_data))]
    
    #Calculate i statistics
    i_mean <- mean(i_data$Count)
    i_sd <- sd(i_data$Count)
    i_median <- median(i_data$Count)
    i_quantiles <- quantile(i_data$Count, c(0.25, 0.5, 0.75), type=1)
    
    #Append gene count results to the general list
    gene_count_results[[length(gene_count_results)+1]] <- i_data
    names(gene_count_results)[length(gene_count_results)] <- paste0(class_levels[i], "_gene_count_data")
    
    gene_count_results[[length(gene_count_results)+1]] <- i_mean
    names(gene_count_results)[length(gene_count_results)] <- paste0(class_levels[i], "_gene_count_mean")
    
    gene_count_results[[length(gene_count_results)+1]] <- i_sd
    names(gene_count_results)[length(gene_count_results)] <- paste0(class_levels[i], "_gene_count_sd")
    
    gene_count_results[[length(gene_count_results)+1]] <- i_median
    names(gene_count_results)[length(gene_count_results)] <- paste0(class_levels[i], "_gene_count_median")
    
    gene_count_results[[length(gene_count_results)+1]] <- i_quantiles
    names(gene_count_results)[length(gene_count_results)] <- paste0(class_levels[i], "_gene_count_quantiles")

  }
  
  #Return output lists
  return(gene_count_results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to summarise gene count results
gene_count_summary <- function(gene_count_samples, results_name) {
  
  #Generate a summary table
  gene_count_stats <- cbind(unlist(gene_count_samples[grep("_gene_count_mean", names(gene_count_samples))]),
                            unlist(gene_count_samples[grep("_gene_count_sd", names(gene_count_samples))]),
                            unlist(gene_count_samples[grep("_gene_count_median", names(gene_count_samples))]),
                            t(matrix(unlist(gene_count_samples[grep("_gene_count_quantiles", names(gene_count_samples))]), nrow=3)))
  #Set row and col names
  rownames(gene_count_stats) <- gsub("_.*", "", rownames(gene_count_stats))                         
  colnames(gene_count_stats) <- c("Mean", "SD", "Median", "Q25%", "Q50%", "Q75%")                         
  xlsx::write.xlsx(gene_count_stats, paste0("Results/", results_name, "_gene_count_stats.xlsx"))
  #Return results
  return(gene_count_stats)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
