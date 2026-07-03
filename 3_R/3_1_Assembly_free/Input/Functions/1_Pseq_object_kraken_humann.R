
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to combine kraken and maxikraken outputs
combine_kraken_taxa <- function(taxa_kraken, taxa_maxikraken, metadata) {
  
  #Set taxonomic clades as rownames
  #Kraken
  rownames(taxa_kraken) <- taxa_kraken$name
  #Maxikraken
  rownames(taxa_maxikraken) <- taxa_maxikraken$name
  
  #Remove unnecessary columns
  #Kraken
  taxa_kraken <- taxa_kraken[,4:ncol(taxa_kraken)]
  #Maxikraken
  taxa_maxikraken <- taxa_maxikraken[,4:ncol(taxa_maxikraken)]
  
  #Clean sample names
  #Kraken
  names(taxa_kraken) <- gsub("_.*", "", names(taxa_kraken))
  #Maxikraken
  names(taxa_maxikraken) <- gsub("_.*", "", names(taxa_maxikraken))
  
  #Prune microbial counts and abundances
  #Kraken
  counts_taxa_kraken <- taxa_kraken[,seq(1,ncol(taxa_kraken), by=2)]
  abundance_taxa_kraken <- taxa_kraken[,seq(2,ncol(taxa_kraken), by=2)]
  #Maxikraken
  #Family
  counts_taxa_maxikraken <- taxa_maxikraken[,seq(1,ncol(taxa_maxikraken), by=2)]
  abundance_taxa_maxikraken <- taxa_maxikraken[,seq(2,ncol(taxa_maxikraken), by=2)]
  
  #Find taxa determined by each database
  common_taxa <- intersect(rownames(abundance_taxa_maxikraken), rownames(abundance_taxa_kraken))
  unique_maxikraken <- setdiff(rownames(abundance_taxa_maxikraken), rownames(abundance_taxa_kraken))
  unique_kraken <- setdiff(rownames(abundance_taxa_kraken), rownames(abundance_taxa_maxikraken))
  #Combine taxonomic assignations from both databases
  for (i in 1:length(common_taxa)) {
    #Find i common taxa
    i_taxa <- common_taxa[i]
    i_taxa_kraken <- rowMeans(counts_taxa_kraken[i_taxa,])
    i_taxa_maxikraken <- rowMeans(counts_taxa_maxikraken[i_taxa,])
    #If bacterial counts are higher in kraken, remove this taxa from maxikraken
    if (i_taxa_kraken >= i_taxa_maxikraken) {
      counts_taxa_maxikraken <- counts_taxa_maxikraken[!(row.names(counts_taxa_maxikraken) %in% i_taxa),]
    }
    #If bacterial counts are higher in maxikraken, remove this taxa from kraken
    if (i_taxa_kraken < i_taxa_maxikraken) {
      counts_taxa_kraken <- counts_taxa_kraken[!(row.names(counts_taxa_kraken) %in% i_taxa),]
    }
  }
  #Check all duplicated taxa are removed
  duplicate_check <- (nrow(counts_taxa_maxikraken)+nrow(counts_taxa_kraken)) == (length(common_taxa)+length(unique_maxikraken)+length(unique_kraken))
  print(paste0("Taxonomic data combined: ", duplicate_check))
  #Combine both datasets
  combined_taxonomy <- rbind.data.frame(counts_taxa_kraken, counts_taxa_maxikraken)
  #Transpose taxonomy table to match metadata format
  taxa_counts <- t(combined_taxonomy)
  
  #Select samples corresponding to this study
  select_samples <- intersect(rownames(metadata), rownames(taxa_counts))
  selected_metadata <- metadata[select_samples, ]
  #Order taxa
  taxa_counts <- taxa_counts[select_samples, ]
  metadata_check <- all(rownames(taxa_counts) == rownames(selected_metadata))
  print(paste0("Same samples in metadata and combined taxonomy file: ", metadata_check))
  #Convert microbial counts to percentages
  taxa_percentage <- taxa_counts*100/rowSums(taxa_counts)
  # taxa_percentage <- taxa/rowSums(taxa)
  print(paste0("Sum of taxa abundances in each sample: ", all(rowSums(taxa_percentage))))
  
  #Return results
  results <- list(taxa_counts, taxa_percentage, selected_metadata)
  names(results) <- c("taxa_counts", "taxa_percentage", "selected_metadata")
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to generate pseq objects for gene families and metabolic pathways data
generate_pseq_kraken_humann <- function(pseq_input, selected_metadata, abundance_filter, frequency_filter) {
  
  #Create a list to store pseq objects
  pseq_list <- list()
  
  #Iterate pseq inputs
  for (i in 1:length(pseq_input)) {
    #Get i pseq input
    i_pseq_input <- pseq_input[[i]]
    
    #Process kraken taxonomy
    if (length(grep("kraken", names(pseq_input)[i], ignore.case=T)) != 0) {
      
      #PREPARE TAXA
      #Remove taxa showing near zero variances to avoid errors (specially at strain level)
      i_pseq_input <- i_pseq_input[,-caret:::nearZeroVar(i_pseq_input)] 
      taxa <- t(i_pseq_input)
      #PSEQ OBJECT OF TAXA 
      #Create taxonomy file
      taxonomy <- as.matrix(rownames(taxa))
      colnames(taxonomy) <- c("Species")
      rownames(taxonomy) <- taxonomy[,1]
      #Prepare pseq metadata format
      pseq_metadata <- sample_data(selected_metadata)
      rownames(pseq_metadata) <- pseq_metadata$Shotgun_code
      #Format otu data
      abundances <- as.data.frame(taxa)
      #Create phyloseq object
      pseq_taxa <- phyloseq(tax_table(taxonomy), 
                            otu_table(abundances, 
                                      taxa_are_rows=TRUE), 
                            pseq_metadata)
      #Filter abundance and prevalence
      if (all(unique(colSums(pseq_taxa@otu_table)) == 100)) {
        #Filter abundance data
        pseq_taxa_filtered <- core(pseq_taxa, detection=abundance_filter, prevalence=frequency_filter)
      } else {
        #Convert pseq counts to abundances
        pseq_taxa_abundance <- pseq_taxa %>%  
          microbiome::transform(transform = "compositional")
        #Get abundances (%)
        pseq_taxa_abundance@otu_table <- pseq_taxa_abundance@otu_table*100
        #Filter abundance data
        pseq_taxa_abundance_filtered <- core(pseq_taxa_abundance, detection=abundance_filter, prevalence=frequency_filter)
        pseq_taxa_filtered <- prune_taxa(rownames(pseq_taxa_abundance_filtered@otu_table), pseq_taxa)
        }

      #Append i pseq object to pseq list
      pseq_list[[length(pseq_list)+1]] <- pseq_taxa
      names(pseq_list)[[length(pseq_list)]] <- names(pseq_input)[i]
      pseq_list[[length(pseq_list)+1]] <- pseq_taxa_filtered
      names(pseq_list)[[length(pseq_list)]] <- paste0(names(pseq_input)[i], "_filtered")
      
    }
    
    #Process humann taxonomy
    if (length(grep("humann_taxa", names(pseq_input)[i], ignore.case=T)) != 0) {
      
      #PREPARE TAXA
      taxa <- t(i_pseq_input)
      #PSEQ OBJECT OF TAXA 
      #Continue with filtered data
      data_nzv <- as.data.frame(t(taxa))
      #Generate data frame with separated names of taxonomic classes
      taxonomy <- data.frame(Kingdom= numeric(1), Phylum= numeric(1), Class= numeric(1),
                             Order= numeric(1), Family= numeric(1), Genus= numeric(1),
                             Species= numeric(1),
                             stringsAsFactors = FALSE)
      #Iterate dataset
      for (j in 1:length(data_nzv)){
        names_taxonomy <- names(data_nzv)[j]
        
        if (isTRUE(grep("p__", names_taxonomy, value=FALSE)!=0)) {
          names_taxonomy <- unlist(strsplit(names_taxonomy, "p__"))
        }
        
        if (isTRUE(grep("c__", names_taxonomy, value=FALSE)!=0)) {
          names_taxonomy <- unlist(strsplit(names_taxonomy, "c__"))
        }
        
        if (isTRUE(grep("o__", names_taxonomy, value=FALSE)!=0)) {
          names_taxonomy <- unlist(strsplit(names_taxonomy, "o__"))
        }
        
        if (isTRUE(grep("f__", names_taxonomy, value=FALSE)!=0)) {
          names_taxonomy <- unlist(strsplit(names_taxonomy, "f__"))
        }
        
        if (isTRUE(grep("g__", names_taxonomy, value=FALSE)!=0)) {
          names_taxonomy <- unlist(strsplit(names_taxonomy, "g__"))
        }
        
        if (isTRUE(grep("s__", names_taxonomy, value=FALSE)!=0)) {
          names_taxonomy <- unlist(strsplit(names_taxonomy, "s__"))
        }
        
        if (isTRUE(grep("[0-9]", names_taxonomy, value=FALSE)!=0)) {
          names_taxonomy[length(names_taxonomy)] <- paste0(names_taxonomy[length(names_taxonomy)], "_strain")
        }
        
        na <- as.character(rep(NA, length(taxonomy)-length(names_taxonomy)))
        filled_na <- c(names_taxonomy, na)
        taxonomy <- rbind(taxonomy, filled_na)
        
      }
      
      #Clean names for every taxonomy clade
      taxonomy <- taxonomy[-c(1),]
      taxonomy[,1] <- gsub("k__", "", taxonomy[,1])
      for (j in 1:ncol(taxonomy)) {
        taxonomy[,j] <- gsub("[.]", "", taxonomy[,j])
      }
      
      #Reaname taxonomy rows with the last clade
      taxonomy_rows <-c()
      
      for(j in 1:nrow(taxonomy)) {
        for (k in 1:ncol(taxonomy)) {
          
          #Add last identified clade
          if (is.na(taxonomy[j,k]) == TRUE) {
            if (is.na(taxonomy[j,c(k-1)]) == FALSE) {
              taxonomy_rows <- c(taxonomy_rows, taxonomy[j,c(k-1)])
              taxonomy_rows[j] <- paste0("_", taxonomy_rows[j])
              taxonomy_rows[j] <- paste0(j, taxonomy_rows[j])
            }
          }
          
          #Add species if available
          if (k == ncol(taxonomy)) {
            if (is.na(taxonomy[j,k]) == FALSE) {
              taxonomy_rows <- c(taxonomy_rows, taxonomy[j,k])
              taxonomy_rows[j] <- paste0("_", taxonomy_rows[j])
              taxonomy_rows[j] <- paste0(i, taxonomy_rows[j])
            }
          }
          
        }
      }
      #Set taxonomy rownames
      rownames(taxonomy) <- taxonomy_rows
      
      #Substet metaphlan2 output data
      rownames(taxonomy)
      names(data_nzv) <- rownames(taxonomy)
      
      #Select kingdom data
      kingdom <- taxonomy[which(is.na(taxonomy$Phylum)), ]
      kingdom_abundances <- data_nzv[,rownames(kingdom)]
      
      #Select phylum data
      phylum <- taxonomy[which(is.na(taxonomy$Class)), ]
      rownames_non_phylum <- rownames(kingdom)
      phylum <- phylum[! rownames(phylum) %in% rownames_non_phylum, drop = F, ]
      phylum_abundances <- data_nzv[,rownames(phylum)]
      
      #Select class data
      class <- taxonomy[which(is.na(taxonomy$Order)), ]
      rownames_non_class <- c(rownames(kingdom), rownames(phylum))
      class <- class[! rownames(class) %in% rownames_non_class, drop = F, ]
      class_abundances <- data_nzv[,rownames(class)]
      
      #Select order data
      order <- taxonomy[which(is.na(taxonomy$Family)), ]
      rownames_non_order <- c(rownames(kingdom), rownames(phylum),
                              rownames(class))
      order <- order[! rownames(order) %in% rownames_non_order, drop = F, ]
      order_abundances <- data_nzv[,rownames(order)]
      
      #Select family data
      family <- taxonomy[which(is.na(taxonomy$Genus)), ]
      rownames_non_family <- c(rownames(kingdom), rownames(phylum),
                               rownames(class), rownames(order))
      family <- family[! rownames(family) %in% rownames_non_family, drop = F, ]
      family_abundances <- data_nzv[,rownames(family)]
      
      #Select genus data
      genus <- taxonomy[which(is.na(taxonomy$Species)), ]
      rownames_non_genus <- c(rownames(kingdom), rownames(phylum),
                              rownames(class), rownames(order),
                              rownames(family))
      genus <- genus[! rownames(genus) %in% rownames_non_genus, drop = F, ]
      genus_abundances <- data_nzv[,rownames(genus)]
      
      #Select species data
      species <- taxonomy
      rownames_non_species <- c(rownames(kingdom), rownames(phylum),
                                rownames(class), rownames(order),
                                rownames(family), rownames(genus))
      species <- species[! rownames(species) %in% rownames_non_species, drop = F, ]
      species_abundances <- data_nzv[,rownames(species)]
      
      #Select strain data
      rownames_strain <- grep("_strain", rownames(species), value=TRUE)
      strain <- species[rownames_strain,]
      strain_abundances <- data_nzv[,rownames(strain)]
      
      #Create pseq objects to study the selected taxonomic class:
      #Complete data: taxonomy and data_nzv (OTU)
      #Individual groups:
      #kingdom, phylum, class, order, family, genus, species, strain
      pseq_taxonomy <- taxonomy
      pseq_abundances <- data_nzv
      
      #"OTU" Table
      otu <- t(pseq_abundances)
      taxonomy <- as.matrix(pseq_taxonomy)
      rownames(otu) <- rownames(taxonomy)
      
      #Multiply percentages * reads yo get counts and round the result
      otu_counts <- otu
      for (j in 1:ncol(otu)) {
        otu_counts[,j] <- otu_counts[,j]*selected_metadata$Reads_Knead[j]
      }
      
      #Scale OTU data: multiply relative proportion per 1e+06
      otu_counts <- otu_counts*1e+06/max(otu_counts)
      otu_counts <- round(otu_counts)
      
      #Prepare pseq metadata format
      pseq_metadata <- sample_data(selected_metadata)
      rownames(pseq_metadata) <- gsub("_", "", pseq_metadata$Shotgun_code)
      
      #Create phyloseq object
      #Abundance
      pseq_taxa_abundance <- phyloseq(tax_table(taxonomy), otu_table(otu, taxa_are_rows=TRUE), 
                                   pseq_metadata)
      #Counts
      pseq_taxa_counts <- phyloseq(tax_table(taxonomy), otu_table(otu_counts, taxa_are_rows=TRUE), 
                            pseq_metadata)
      #Filter abundance and prevalence
      pseq_taxa_abundance_filtered <- core(pseq_taxa_abundance, detection=abundance_filter, prevalence=frequency_filter)
      pseq_taxa_counts_filtered <- prune_taxa(rownames(pseq_taxa_abundance_filtered@otu_table), pseq_taxa_counts)
      
      #Append i pseq object to pseq list
      pseq_list[[length(pseq_list)+1]] <- pseq_taxa_abundance
      names(pseq_list)[[length(pseq_list)]] <- c("humann_taxa_abundance")
      pseq_list[[length(pseq_list)+1]] <- pseq_taxa_counts
      names(pseq_list)[[length(pseq_list)]] <- c("humann_taxa_counts")
      pseq_list[[length(pseq_list)+1]] <- pseq_taxa_abundance_filtered
      names(pseq_list)[[length(pseq_list)]] <- c("humann_taxa_abundance_filtered")
      pseq_list[[length(pseq_list)+1]] <- pseq_taxa_counts_filtered
      names(pseq_list)[[length(pseq_list)]] <- c("humann_taxa_counts_filtered")
      
    }
      
    #Process humann functional data
    if (length(grep("humann_func_", names(pseq_input)[i], ignore.case=T)) != 0) {
      
      #PREPARE GENES FAMILIES 
      genes_families <- i_pseq_input
      #PSEQ OBJECT OF GENES FAMILIES 
      #Create "OTU" Table with genes associated to microbial species
      taxonomy_genes_families <- grep("s__", rownames(genes_families))
      final_genes_families <- genes_families[taxonomy_genes_families,]
      rownames(final_genes_families) <- gsub("UniRef*.0_", "", rownames(final_genes_families))
      rownames(final_genes_families) <- gsub("[.]", "|", rownames(final_genes_families))
      #Create taxonomy file
      taxonomy_families <- rownames(final_genes_families)
      taxonomy_families <- as.matrix(gsub(".*s__","",taxonomy_families))
      colnames(taxonomy_families) <- c("Species")
      rownames(taxonomy_families) <- rownames(final_genes_families)
      #Prepare pseq metadata format
      pseq_metadata <- sample_data(selected_metadata)
      rownames(pseq_metadata) <- gsub("_", "", pseq_metadata$Shotgun_code)
      #Create phyloseq object
      pseq_families <- phyloseq(tax_table(taxonomy_families), 
                                otu_table(final_genes_families, taxa_are_rows=TRUE), 
                                pseq_metadata)
      #Filter abundance and prevalence
      pseq_families_filtered <- core(pseq_families, detection=1, prevalence=frequency_filter)
      
      #Append i pseq object to pseq list
      pseq_list[[length(pseq_list)+1]] <- pseq_families
      names(pseq_list)[[length(pseq_list)]] <- names(pseq_input)[i]
      pseq_list[[length(pseq_list)+1]] <- pseq_families_filtered
      names(pseq_list)[[length(pseq_list)]] <- paste0(names(pseq_input)[i], "_filtered")
      
    }
      
  }

  #Return results
  return(pseq_list)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
