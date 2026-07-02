
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to process HMMER domains
process_hmmer <- function(selected_file_lists, selected_file_paths, remove_taxa_without_domains,
                         metadata, strain_seq, dendo_class) {
  
  #Check the number of output files is the same for all databases
  if (length(unique(lengths(selected_file_lists))) == 1) {
    #Check the order of genomes is the same for all databases
    genomes_names <- list()
    genomes_raw_files <- list()
    #Iterate the results of each database
    for (i in 1:length(selected_file_lists)) {
      #Get i file list
      i_file_list <- selected_file_lists[[i]]
      #Read individual files from i file list and store them in a raw file list
      i_raw_files <- list()
      for (j in 1:length(i_file_list)) {
        i_raw_files[[j]] <- read.table(paste0(selected_file_paths[i], i_file_list[j]), header=F, skip=3)
        names(i_raw_files)[length(i_raw_files)] <- as.character(i_file_list[j])
      }
      #Append i raw files to a single list
      genomes_raw_files[[length(genomes_raw_files)+1]] <- i_raw_files
      names(genomes_raw_files)[[length(genomes_raw_files)]] <- names(selected_file_lists)[i]
      #Get i genomes names
      i_genome_names <- selected_file_lists[[i]]
      # i_genome_names <- gsub("ref_aaORFs.*", "", i_genome_names) #Weisella CERELA: avoid errors
      i_genome_names <- gsub("_aaORFs.*", "", i_genome_names)
      #Append i genomes names to a single list
      genomes_names[[length(genomes_names)+1]] <- i_genome_names
      names(genomes_names)[[length(genomes_names)]] <- names(selected_file_lists)[i]
    }
    #Check the order of genomes is the same for all databases using a custom function
    all_identical <- function(l) all(mapply(identical, head(l, 1), tail(l, -1)))
    #Continue the analysis if the order of genomes is the same for all databases
    if (all_identical(genomes_names) == T) {
      #Create a new list to store combined results of multiple databases
      combined_files <- list()
      #Iterate all genomes
      for (k in 1:length(genomes_raw_files[[1]])) {
        #Get k genome
        k_genomes <- paste0("genomes_raw_files$", names(genomes_names), "[[", k, "]]", collapse=", ")
        #Combine the results of k genome from multiple databases using "parse" and "eval" options
        k_genome_combine_results <- eval(parse(text=paste0("rbind(", k_genomes, ")")))
        combined_files[[length(combined_files)+1]] <- k_genome_combine_results
        names(combined_files)[[length(combined_files)]] <- sub("(.)", "\\U\\1", gsub("_aaORFs.*", "", names(genomes_raw_files[[1]])[k]), perl=T)
      }
    }
  }
  
  #Gather codes indicating domain families from all samples present in the list containing filtered domain data
  domain_data <- data.frame()
  for (i in 1:length(combined_files)) {
    domain_data <- rbind(domain_data, as.data.frame(combined_files[[i]]$V4))
  }
  
  #Remove duplicated domain codes from the domain data frame
  domain_data <- unique(domain_data)
  names(domain_data) <- c("Domain")
  
  #Search for interesting functional domains using manually-curated pfam codes
  selected_codes <- na.omit(unlist(selected_domains$Domain))
  
  #Prune domain data to select interesting functional domains
  domain_data$Domain <- gsub("[.].*", "", domain_data$Domain)
  rownames(domain_data) <- domain_data$Domain
  domain_data <- domain_data[intersect(selected_codes, rownames(domain_data)),, drop=F]
  
  #Generate a new data frame summarising the number of enzyme families per genome sequence
  for (i in 1:length(combined_files)) {
    
    #Determine enzyme families present in i genome sequence
    i_domains <- intersect(domain_data$Domain, gsub("[.].*", "", combined_files[[i]]$V4))
    i_domains <- data.frame(i_domains)
    
    #Continue the analysis if any enzyme are found in i genome sequence
    if (length(i_domains) > 0) {
      names(i_domains) <- names(combined_files)[i]
      
      #If the number of domains found in i genome is lower than total number of domains found in the complete dataset, 
      #Then, add NAs corresponding to the missing domains
      if (nrow(i_domains) < nrow(domain_data)) {
        i_add_nas <- nrow(domain_data) - nrow(i_domains)
        i_domains[nrow(i_domains)+i_add_nas,] <- NA
      }
      
      #Align domains found in i genome to match the total number of domains found in the complete dataset
      i_match_index <- match(domain_data$Domain, i_domains[,1])
      i_domains_match  <- i_domains[i_match_index,]
      
      #glycosidases found in i genomes to the complete dataset
      domain_data <- cbind(domain_data, i_domains_match)
      #Set the column name corresponding to i genome
      names(domain_data)[ncol(domain_data)] <- names(i_domains)
      
    }
  }
  
  #Format the the data frame summarising the number of domains per genome sequence
  #Set row names (i.e. functional domain) 
  rownames(domain_data) <- gsub("[.].*", "", domain_data$Domain)
  #Remove column containing domain names 
  domain_data <- domain_data[,2:ncol(domain_data)]
  #Replace characters in cells by 1
  domain_data_numeric <- domain_data
  domain_data_numeric <- domain_data_numeric %>% mutate_all(funs(str_replace_all(., "^[:alnum:].*", "1")))
  #Replace NAs by 0
  domain_data_numeric[is.na(domain_data_numeric)] <- 0
  #Convert the data frame to numeric
  domain_data_numeric <- mutate_all(domain_data_numeric, function(x) as.numeric(as.character(x)))
  
  #Filter genomes not showing any selected domain
  if (remove_taxa_without_domains == T) {
    if (length(which(colSums(domain_data_numeric) == 0)) != 0) {
      domain_data_numeric <- domain_data_numeric[,-which(colSums(domain_data_numeric) == 0)]
    }
  }
  
  #Remove non-matching domains without subfamilies
  selected_domain_metadata <- selected_domains[match(rownames(domain_data_numeric), selected_domains$Domain),]
  all(selected_domain_metadata$Family == rownames(domain_data_numeric))
  
  #Order annotation
  selected_domain_data_numeric <- as.data.frame(t(domain_data_numeric))
  
  #Select sequence files showing enzyme domains without duplicates
  sequences_annotated_without_duplicates <- intersect(rownames(selected_domain_data_numeric), strain_seq)
  all_removed_strains_are_dupli <- length(grep(paste0("^", setdiff(rownames(selected_domain_data_numeric), strain_seq), "$", collapse="|"),
                                               duplicated_sequences)) == length(setdiff(rownames(selected_domain_data_numeric), strain_seq))
  print(paste0("All srains removed were duplicates: ", all_removed_strains_are_dupli))
  metadata <- metadata[grep(paste0(sequences_annotated_without_duplicates, collapse="|"), metadata$File),]
  
  #Get annotated and non-duplicated strain sequences again
  strain_seq <- metadata$File
  strain_seq <- gsub(".fna$", "", strain_seq)
  strain_seq <- gsub(".fasta$", "", strain_seq)
  length(intersect(strain_seq, rownames(selected_domain_data_numeric))) == length(strain_seq)
  selected_domain_data_numeric <- selected_domain_data_numeric[strain_seq,]
  #Add strain codes
  rownames(selected_domain_data_numeric) <- gsub(" = .*", "", metadata$Taxonomy_metadata)
  selected_domain_data_numeric$Class <- metadata[[dendo_class]]
  
  #Remove alternative strain names
  rownames(selected_domain_data_numeric) <- gsub(".* / ", "", rownames(selected_domain_data_numeric))
  
  #Add strain abbreviations
  substr(rownames(selected_domain_data_numeric), 4, 4) <- "."
  rownames(selected_domain_data_numeric) <- sub("[.].*? ", ". ", rownames(selected_domain_data_numeric))
  
  #Set factor levels
  selected_domain_metadata$Group <- factor(selected_domain_metadata$Group, levels=unique(selected_domain_metadata$Group))
  rownames(selected_domain_metadata) <- selected_domain_metadata$Domain
  selected_domain_metadata <- selected_domain_metadata[,2:ncol(selected_domain_metadata)]
  
  #Get dendogram species and strains
  #Format dendogram to mix both italics and normal fonts in each label
  #https://stackoverflow.com/questions/63788793/how-to-apply-both-italic-and-normal-fonts-in-the-same-label-in-phylogenetic-tree
  #https://r-graph-gallery.com/package/ggtext.html
  #https://stackoverflow.com/questions/72027474/how-to-format-selected-axis-text-words-with-italic-and-bold-using-ggtext
  selected_df <- selected_domain_data_numeric
  species <- ifelse(is.na(word(rownames(selected_df), 1, 2)), rownames(selected_df), word(rownames(selected_df), 1, 2))
  n_words <- str_count(rownames(selected_df), " ") + 1
  strains <- ifelse(is.na(word(rownames(selected_df), 3, n_words)), rownames(selected_df), word(rownames(selected_df), 3, n_words))
  format_species_strains <- paste0("paste(italic('", species, "'),", " ' ", strains, "')")
  format_species_strains <- paste0("*", species, "* ", "**", strains, "**")
  rownames(selected_domain_data_numeric) <- format_species_strains
  
  #Set strain class as the first column
  selected_domain_data_numeric <- selected_domain_data_numeric[,c(length(selected_domain_data_numeric), 1:(length(selected_domain_data_numeric)-1))]
  
  #Append results
  results <- list(selected_domain_data_numeric, selected_domain_metadata)
  names(results) <- c("selected_domain_data_numeric", "selected_domain_metadata")
  #Return results
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to remove empty rows and columns from the heatmap
remove_empty_cols_rows_from_heatmap <- function(selected_heatmap_data) {
  #Remove empty cols
  empty_cols <- which(colSums(selected_heatmap_data[,2:ncol(selected_heatmap_data)]) == 0)
  if (length(empty_cols) != 0) {
    selected_heatmap_data <- selected_heatmap_data[,-(empty_cols+1)]
  }
  #Remove empty rows
  empty_rows <- which(rowSums(selected_heatmap_data[,2:ncol(selected_heatmap_data)]) == 0)
  if (length(empty_rows) != 0) {
    selected_heatmap_data <- selected_heatmap_data[-empty_rows,]
  }
  #Return clean dataframe
  return(selected_heatmap_data)
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Finish
