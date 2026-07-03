
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to perform microbiota statistical analyses (differential abundance)
differential_abundance_tests <- function(j_i_data, categorical_var, p_threshold, padj_threshold,
                                         run_deseq, run_microbiomemarker, run_maaslin, analyze_genes, j_i_name, i_name) {
  
  #Advanced models of differential abundance - for abundant data
  print("Running differential abundance tests...")
  #Deseq2
  if (run_deseq == T) {
    #In order to run Deseq2 method:
    print("Running Deseq2...")
    #Add 1 pseudocount to all gene families and metabolic pathways abundance in order to avoid zeros for DESEQ2 calculation
    j_i_data@otu_table <- j_i_data@otu_table+1
    #DESeq2 abundance testing
    #Start by converting phyloseq object to deseq2 format
    ds2 <- phyloseq_to_deseq2(j_i_data, formula(paste("~ ", categorical_var)))
    #Run DESeq2 analysis (all taxa at once!)
    dds <- try(DESeq(ds2))
    if (class(dds) == "try-error") {
      df_filt <- NULL
    } else {
      #Alternative method to avoid errors
      # dds <- estimateSizeFactors(ds2)
      # dds <- estimateDispersionsGeneEst(dds)
      # dispersions(dds) <- mcols(dds)$dispGeneEst
      #Investigate results
      res <- results(dds)
      deseq.results <- as.data.frame(res)
      df <- deseq.results
      df$taxon <- rownames(df)
      df <- df %>% arrange(log2FoldChange, padj)
      #Select genes with significant p-value adj
      df_filt <- df[which(df$pvalue > 0.0 & df$pvalue < p_threshold),]
      df_filt <- df[which(df$padj > 0.0 & df$padj < padj_threshold),]
    }
    #Remove pseudocounts
    j_i_data@otu_table <- j_i_data@otu_table-1
  }
  
  #Microbiomemarker
  if (run_microbiomemarker == T) {
    #Format input pseq object
    print("Running Microbiomemarker...")
    stats_input <- j_i_data
    # stats_input <- pseq_div
    #Format input pseq object
    # colnames(stats_input@tax_table) <- c("Species")
    #Normalize data prior statistical analysis
    # stats_input <- normalize(stats_input, method="TSS")
    #Run multiple statistical tests
    # aldex <- run_aldex(stats_input, group=categorical_var)
    #ANCOM
    print("Microbiomemarker ANCOM")
    ancom <- run_ancom(stats_input, group=categorical_var)
    #ANCOMBC
    print("Microbiomemarker ANCOMBC")
    ancombc <- run_ancombc(stats_input, group=categorical_var)
    #LEFSE
    print("Microbiomemarker LEFSE")
    lefse <- try(run_lefse(stats_input, group=categorical_var), silent=T)
    if (class(lefse) == "try-error") {
      setClass("results", representation(name="character", marker_table="data.frame"))
      lefse <- new("results", name="", marker_table=data.frame())
    }
    #METAGENOMESEQ
    print("Microbiomemarker METAGENOMESEQ")
    metagenomeseq <- try(run_metagenomeseq(stats_input, group=categorical_var, method="ZIG"))
    if (class(metagenomeseq) == "try-error") {
      setClass("results", representation(name="character", marker_table="data.frame"))
      metagenomeseq <- new("results", name="", marker_table=data.frame())
    }
    
    #Select significant taxa
    marker_taxa <- unique(sort(c(
      # aldex@marker_table$feature,
      ancom@marker_table$feature,
      ancombc@marker_table$feature,
      lefse@marker_table$feature,
      metagenomeseq@marker_table$feature
    )))
    #Format microbiomeMarker outputs
    marker_taxa_format <- gsub("\\|", "--", marker_taxa)
    marker_taxa_format <- gsub("--$", "", marker_taxa_format)
    marker_taxa_format <- gsub("--.?__.?__.*", "", marker_taxa_format)
    marker_taxa_format <- gsub(".*--", "", marker_taxa_format)
    marker_taxa_format <- gsub("^.?__", "", marker_taxa_format)
    marker_taxa_format <- gsub("_.?__.*", "", marker_taxa_format)
  }
  
  #MaAsLiN2
  if (run_maaslin == T) {
    #Format input pseq object
    print("Running MaAsLiN2...")
    maaslin_input <- j_i_data
    #Create a new vector to store significant taxa
    all_maaslin_taxa <- c()
    #Iterate factor levels to set multiple references
    for (k in 1:length(levels(maaslin_input@sam_data[[categorical_var]]))) {
      #Get k level
      k_level <- levels(maaslin_input@sam_data[[categorical_var]])[k]
      #Run MaAsLiN2
      #Normalization method: CLR more restrictive than CSS, TMM and TSS and NONE
      k_maaslin=Maaslin2(
        input_data=maaslin_input@otu_table, plot_heatmap=F, plot_scatter=F,
        input_metadata=data.frame(maaslin_input@sam_data),
        output=paste0("Maaslin/Maaslin_", i_name, "_", j_i_name, "_ref_", k_level),
        min_prevalence=0.001, max_significance=0.25,
        fixed_effects=categorical_var, reference=paste0(categorical_var, ",", k_level))
      #Read results
      k_maaslin_results <- read.table(file=paste0("Maaslin/Maaslin_", i_name, "_", j_i_name, "_ref_", k_level, "/significant_results.tsv"), 
                                      sep='\t', header=TRUE)
      #Select significant taxa
      k_maaslin_results <- subset(k_maaslin_results, k_maaslin_results$qval < 0.25)
      k_maaslin_taxa <- k_maaslin_results$feature
      #Append the results to a single vector
      all_maaslin_taxa <- sort(unique(c(all_maaslin_taxa, k_maaslin_taxa)))
    }
  }
  
  #Combine results from different statistical tests
  marker_intersect <- c()
  #DESeq2 results
  if (exists("df_filt")) {
    if (analyze_genes == T) {
      marker_intersect <- c(marker_intersect, df_filt$taxon)
    } else {
      marker_intersect <- c(marker_intersect, gsub(".*remove", "", sub("_", "remove", df_filt$taxon)))
    }
  }
  #Microbiomemarker results
  if (exists("marker_taxa_format")) {
    if (analyze_genes == T) {
      print("Microbiomemarker not suitable for functional data")
    } else {
      marker_intersect <- c(marker_intersect, marker_taxa_format)
    }
  }
  #MaAsLiN2 results
  if (exists("all_maaslin_taxa")) {
    if (analyze_genes == T) {
      marker_intersect <- c(marker_intersect, all_maaslin_taxa)
    } else {
      marker_intersect <- c(marker_intersect, all_maaslin_taxa)
    }
  }
  #Remove duplicates
  marker_intersect <- sort(unique(marker_intersect))
  
  #Return results
  return(marker_intersect)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to calculate taxa increments and filter genes with taxa
increments_and_genes_with_taxa <- function(j_i_data, significant_rownames, selected_clades, calculate_increments, increment_threshold, filter_genes, i_name) {
  
  #Create a list to store results
  results <- list()
  #Calculate increments based on taxonomic data
  if (calculate_increments == T) {
    #Select only those gene families and metabolic pathways associated to relevant taxa
    taxa_per_level <- list()
    rowmeans_taxa_per_level <- list()
    #Calculate percentages
    if (any(colSums(j_i_data@otu_table) > 100)) {
      #Convert pseq counts to abundances
      j_i_data_abundance <- j_i_data %>%  
        microbiome::transform(transform = "compositional")
      #Get abundances (%)
      j_i_data_abundance@otu_table <- j_i_data_abundance@otu_table*100
    } else {
      j_i_data_abundance <- j_i_data
    }
    
    #Iterate factor levels
    for (l in 1:length(levels(j_i_data_abundance@sam_data[[i_name]]))) {
      #Prune l level
      l_level <- levels(j_i_data_abundance@sam_data[[i_name]])[l]
      l_level_taxa <- prune_samples(j_i_data_abundance@sam_data[[i_name]]==l_level, j_i_data_abundance)
      #Prune significant taxa determined by the statistical analysis
      l_level_taxa <- prune_taxa(unique(significant_rownames), l_level_taxa)
      #Get l rowmeans
      l_rowmeans <- rowMeans(l_level_taxa@otu_table)
      #Append i results to a single list
      #l pseq
      taxa_per_level[[length(taxa_per_level)+1]] <- l_level_taxa
      names(taxa_per_level)[length(taxa_per_level)] <- l_level
      #l rowmeans
      rowmeans_taxa_per_level[[length(rowmeans_taxa_per_level)+1]] <- l_rowmeans
      names(rowmeans_taxa_per_level)[length(rowmeans_taxa_per_level)] <- l_level
    }
    
    #Calculate increments
    increment <- eval(parse(text=paste0("rowmeans_taxa_per_level[[", seq(1:length(rowmeans_taxa_per_level)), "]]", collapse=" - ")))
    #Prune increments showing high values (i.e. higher than increment_threshold%)
    filter_increment <- which(abs(increment) > increment_threshold)
    #Select taxa showing these increments
    increment_taxa <- names(filter_increment)
    #Remove numbers to search for these taxa in functional datasets
    relevant_taxa <- unique(sub(".*[0-9]+_", "", increment_taxa))
    #Remove strains from relevant taxa to avoid errors (i.e. do not include non-relevant taxa)
    if (length(grep("strain", relevant_taxa)) > 0) {
      relevant_taxa <- relevant_taxa[-grep("strain", relevant_taxa)]
    }
    
    #Append the results to a single list
    results[[length(results)+1]] <- relevant_taxa
    names(results)[[length(results)]] <- c("taxa_showing_relevant_increments")
  }
  
  #Get genes and pathways associated to relevant taxa that increment significantly
  if (filter_genes == T) {
    #Filter genes
    genes_with_taxa <- unique(grep(paste(gsub("sp._", "sp_", gsub(" ", "_", selected_clades)),collapse="|"), 
                                   significant_rownames, value=T, ignore.case=T))
    #Subset significant genes and pathways dataframes
    j_i_data_with_taxa <- prune_taxa(genes_with_taxa, j_i_data) 
    
    #Append the results to a single list
    results[[length(results)+1]] <- j_i_data_with_taxa
    names(results)[[length(results)]] <- c("functional_data_with_taxa")
  }
  
  #Return results
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to calculate differential expression statistics in pseq objects contaning functional metagenomics data
differential_tables <- function(pseq_taxa, significant_rownames, metadata_class_no_formula) {
  
  #Filter significant rownames
  if (length(significant_rownames) != 0) {
    #Prepare pseq object with significant taxa
    significant_pseq_taxa <- prune_taxa(significant_rownames, pseq_taxa)
    
    #Determine all levels for metadata class
    class_levels <- levels(as.factor(significant_pseq_taxa@sam_data[[metadata_class_no_formula]]))
    #Generate and empty dataframe to store mean values per class level
    means_levels <- data.frame(matrix(, nrow=nrow(significant_pseq_taxa@otu_table), ncol=0))
    #Subset significant class data to get the means of each group
    for (i in 1:length(class_levels)) {
      #Generate a pseq object for each level
      significant_i <- prune_samples(significant_pseq_taxa@sam_data[[metadata_class_no_formula]] == class_levels[i], significant_pseq_taxa)
      #Calculate means for all functions
      i_means <- data.frame(rowMeans(significant_i@otu_table))
      names(i_means) <- paste0(class_levels[i], "_Mean") 
      #Calculate deviations for all functions
      i_sd <- as.data.frame(apply(significant_i@otu_table, 1, sd)) 
      names(i_sd) <- paste0(class_levels[i], "_SD") 
      #Append i means to common file
      means_levels <- cbind(means_levels, i_means, i_sd) 
    }
    
    #Create different datasets to store functions higher in each metadata class
    #Create an empty list
    max_list <- list()
    #Iterate dataframe containing mean values
    for (i in 1:ncol(means_levels[,grep("_Mean", names(means_levels))])) { 
      #Find functions showing the highest value at i level
      i_max_functions <- which(max.col(means_levels[,grep("_Mean", names(means_levels))]) == i) 
      i_max_functions_df <- means_levels[i_max_functions,]
      #Append functions showing the highest value at i level to the common list
      max_list[[length(max_list)+1]] <- i_max_functions_df
      names(max_list)[i] <- paste("Higher in: ", names(means_levels[,grep("_Mean", names(means_levels))])[i]) 
    }
    
    #Create a new list for taxa data
    taxa_max_list <- max_list
    significant_pseq_taxa_data <- significant_pseq_taxa@otu_table
    significant_pseq_taxa_metadata <- significant_pseq_taxa@sam_data
    
  } else {
    significant_pseq_taxa <- NULL
    significant_pseq_taxa_data <- NULL
    significant_pseq_taxa_metadata <- NULL
    taxa_max_list <- NULL
  }
  
  #Return output list
  results <- list(taxa_max_list, significant_pseq_taxa,
                  significant_pseq_taxa_data, significant_pseq_taxa_metadata)
  names(results) <- c("tables_per_level", "pseq_object_sig_filtered", 
                      "otu_table_sig_filtered", "metadata_sig_filtered")
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to summarise differential expression statistics in pseq objects contaning functional metagenomics data
differential_tables_summary <- function(selected_metadata, metadata_class_no_formula, dif_tables, pseq_taxa, results_name,
                                        analyze_taxonomy, analyze_families, analyze_pathways) {
  
  #Create new tables to summarize statistical differences
  taxa_summary <- data.frame()
  families_species_summary <- data.frame()
  pathways_species_summary <- data.frame()
  pathways_function_summary <- data.frame()
  pathways_function_summary_simple <- data.frame()
  
  #Iterate factor levels to export the results
  for (i in 1:length(levels(selected_metadata[[metadata_class_no_formula]]))) {
    #Get i level
    i_level <- levels(selected_metadata[[metadata_class_no_formula]])[i]
    #Get i differential tables
    i_taxa <- dif_tables$tables_per_level[[grep(paste0(" ", i_level, "_Mean"), names(dif_tables$tables_per_level))]]
    #Prevent errors if taxa/families/pathways are missing
    if (length(i_taxa) == 0) {
      i_taxa <- data.frame(matrix(nrow=0, ncol=0))
    }
    #Export differential tables
    #Taxa
    if (analyze_taxonomy == T) {
      if (nrow(i_taxa) > 0) {
        i_taxa <- i_taxa[order(i_taxa[[paste0(i_level, "_Mean")]], decreasing=T),]
        #Add taxonomic levels to i taxa
        i_taxa_table <- pseq_taxa@tax_table[grep(paste0(rownames(i_taxa), collapse="|"), rownames(pseq_taxa@tax_table)),]
        i_taxa_table <- i_taxa_table[rownames(i_taxa),]
        i_taxa_level <- c()
        for (j in 1:nrow(i_taxa_table)) {
          j_i_taxa_level <- colnames(i_taxa_table[j, which(i_taxa_table[j,] != "NA")[length(which(i_taxa_table[j,] != "NA"))]])
          i_taxa_level <- c(i_taxa_level, j_i_taxa_level)
        }
        i_taxa$'Taxonomic level:' <- i_taxa_level 
        i_taxa$Taxa <- gsub("_", " ", gsub(".*replace", "", sub("_", "replace", rownames(i_taxa))))
        i_taxa$'Higher in:' <- i_level
        i_taxa <- i_taxa[,c(ncol(i_taxa), ncol(i_taxa)-1, ncol(i_taxa)-2, 1:(ncol(i_taxa)-3))]
        colnames(i_taxa) <- gsub("_Mean", " (Mean)", colnames(i_taxa))
        colnames(i_taxa) <- gsub("_SD", " (SD)", colnames(i_taxa))
        #Sort by taxonomic levels
        i_taxa$`Taxonomic level:` <- factor(i_taxa$`Taxonomic level:`,
                                            levels=colnames(pseq_taxa@tax_table))
        i_taxa <- i_taxa[order(i_taxa$`Taxonomic level:`),]
        #Append i taxa to a single data.frame
        taxa_summary <- rbind(taxa_summary, i_taxa)
        #Export tables
        xlsx::write.xlsx(taxa_summary, paste0("Results/", results_name, "_taxa_summary.xlsx"))     
      } else {
        i_taxa <- NULL
      }
    }
    
    #Genes
    if (analyze_families == T && is.null(i_taxa) == F) {
      if (nrow(i_taxa) > 0) {
        #Genes - contributing species
        i_families_species <- table(gsub("*.__", "", gsub(".*[|]", "", rownames(i_taxa)))) %>% 
          as.data.frame() %>% 
          arrange(desc(Freq)) %>% `colnames<-`(c("Taxa", "Frequency"))
        i_families_species$'Higher in:' <- i_level
        i_families_species <- i_families_species[,c(ncol(i_families_species), 1:(ncol(i_families_species)-1))]
        i_families_species$Taxa <- gsub("_", " ", i_families_species$Taxa)
        #Append i taxa to a single data.frame
        families_species_summary <- rbind(families_species_summary, i_families_species)
        #Export tables
        xlsx::write.xlsx(families_species_summary, paste0("Results/", results_name, "_families_species_summary.xlsx"))
      } else {
        families_species_summary <- i_families_species
      }
    }
    
    #Pathways
    if (analyze_pathways == T && is.null(i_taxa) == F) {
      if (nrow(i_taxa) > 0) {
        #Pathways - contributing species
        i_pathways_species <- table(gsub("*.__", "", gsub(".*[|]", "", rownames(i_taxa)))) %>%
          as.data.frame() %>%
          arrange(desc(Freq)) %>% `colnames<-`(c("Taxa", "Frequency"))
        i_pathways_species$'Higher in:' <- i_level
        i_pathways_species <- i_pathways_species[,c(ncol(i_pathways_species), 1:(ncol(i_pathways_species)-1))]
        i_pathways_species$Taxa <- gsub("_", " ", i_pathways_species$Taxa)
        #Append i taxa to a single data.frame
        pathways_species_summary <- rbind(pathways_species_summary, i_pathways_species)
        
        #Pathways - contributing functions
        i_pathways_function <- table(gsub(".*:", "", gsub("[|].*", "", rownames(i_taxa)))) %>%
          as.data.frame() %>%
          arrange(desc(Freq)) %>% `colnames<-`(c("Pathway", "Frequency")) 
        i_pathways_function$'Higher in:' <- i_level
        i_pathways_function$Pathway <- gsub("&", "", i_pathways_function$Pathway)
        i_pathways_function$Pathway <- gsub("\\(.*\\)", "", i_pathways_function$Pathway)
        i_pathways_function$Pathway <- gsub("^ ", "", i_pathways_function$Pathway)
        i_pathways_function$Pathway <- gsub("UNINTEGRATED", "unidentified", i_pathways_function$Pathway)
        
        #Pathways - contributing functions (simple)
        #Summarise metabolic function using specific terms
        #Amino acid metabolism
        aminoacid_terms <- c(
          "alanine", "arginine", "asparagine", "aspartic acid", "cysteine", 
          "glutamic acid", "glutamine", "glycine", "histidine", "isoleucine",
          "leucine", "lysine", "methionine", "phenylalanine", "proline",
          "serine", "threonine", "tryptophan", "tyrosine", "valine",
          "selenocysteine", "pyrrolysine",
          
          "phenylalanine", "valine", "tryptophan", "threonine", "isoleucine",
          "methionine", "histidine", "leucine", "lysine",
          
          "amino acid", "chorismate", "ornithine"
        )
        #Nucleotide metabolism
        nucleotide_terms <- c(
          "adenine", "adenosine", "guanosine", "inosine", "queuosine",
          "nucleotide", "ribonucleotide", "ribonucleoside", "pyrimidine",
          "UMP", "UDP"
        )
        #Carbohydrate metabolism
        carbohydrate_terms <- c(
          "acetylneuraminate", "acetylglucosamine", "arabinose", "fucose", "fructose", "galactose", "galacturonate", 
          "glucose", "glucuronic", "glucuronoside", "glycolysis", "glycogen", "lactose", "lactate", "mannose",
          "pectic", "pectin", "pentose phosphate", "rhamnose", "sucrose", "stachyose", "starch", "xylose"
          
        )
        #Vitamin metabolism
        vitamin_terms <- c("folate", "thiamine", "flavin", "coenzyme A")
        #Lipid metabolism
        lipid_terms <- c("methylerythritol")
        #Cell wall metabolism
        cell_wall_terms <- c("peptidoglycan")
        
        #Find these terms
        i_pathways_function$Metabolism <- c(rep("Other", nrow(i_pathways_function)))
        i_pathways_function$Metabolism[grep(paste0(aminoacid_terms, collapse="|"), i_pathways_function$Pathway, ignore.case=T)] <- c("Amino acids")
        i_pathways_function$Metabolism[grep(paste0(nucleotide_terms, collapse="|"), i_pathways_function$Pathway, ignore.case=T)] <- c("Nucleotides")
        i_pathways_function$Metabolism[grep(paste0(carbohydrate_terms, collapse="|"), i_pathways_function$Pathway, ignore.case=T)] <- c("Carbohydrates")
        i_pathways_function$Metabolism[grep(paste0(vitamin_terms, collapse="|"), i_pathways_function$Pathway, ignore.case=T)] <- c("Vitamins")
        i_pathways_function$Metabolism[grep(paste0(lipid_terms, collapse="|"), i_pathways_function$Pathway, ignore.case=T)] <- c("Lipids")
        i_pathways_function$Metabolism[grep(paste0(cell_wall_terms, collapse="|"), i_pathways_function$Pathway, ignore.case=T)] <- c("Cell wall")
        #Order columns
        i_pathways_function <- i_pathways_function[,c(ncol(i_pathways_function)-1, ncol(i_pathways_function), 1:(ncol(i_pathways_function)-2))]
        
        #Group metabolic functions by simple categories
        i_pathways_function_simple <- i_pathways_function %>% group_by(Metabolism, 
                                                                       `Higher in:`) %>%
          summarise_each(funs(paste(., collapse = " ; ")))
        #Sum the number of occurrences of each metabolic category
        for (k in 1:nrow(i_pathways_function_simple)) {
          k_occurrences <- as.numeric(eval(parse(text=gsub(" ; ", " + ", i_pathways_function_simple$Frequency[k]))))
          i_pathways_function_simple$Frequency[k] <- k_occurrences
        }
        #Order table by decreasing frequencies
        i_pathways_function_simple$Frequency <- as.numeric(i_pathways_function_simple$Frequency)
        i_pathways_function_simple <- i_pathways_function_simple[order(i_pathways_function_simple$Frequency, decreasing=T),]
        i_pathways_function_simple <- i_pathways_function_simple[,c(2,1,3:ncol(i_pathways_function_simple))]
        #Append i taxa to a single data.frame
        pathways_function_summary <- rbind(pathways_function_summary, i_pathways_function)
        pathways_function_summary_simple <- rbind(pathways_function_summary_simple, i_pathways_function_simple)
        
        #Export tables
        xlsx::write.xlsx(pathways_species_summary, paste0("Results/", results_name, "_pathways_species_summary.xlsx"))               
        xlsx::write.xlsx(pathways_function_summary, paste0("Results/", results_name, "_pathways_function_summary.xlsx"))               
        xlsx::write.xlsx(as.data.frame(pathways_function_summary_simple), paste0("Results/", results_name, "_pathways_function_summary_simple.xlsx"))  
      } else {
        i_pathways_species <- NULL
        i_pathways_function_simple <- NULL
        i_pathways_function_simple <- NULL
      }
    }
  } 
  
  #Return output list
  results_list <- list(taxa_summary, families_species_summary, pathways_species_summary,
                       pathways_function_summary, as.data.frame(pathways_function_summary_simple))
  names(results_list) <- c("taxa_summary", "families_species_summary", "pathways_species_summary",
                           "pathways_function_summary", "pathways_function_summary_simple")
  return(results_list)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Perform multiple microbiome statistical tests using a custom function
microbiota_multiple_tests <- function(subset_list, categorical_var, increment_threshold, p_threshold, padj_threshold, all_clades) {
  
  #Create a list to summarise the results
  results <- list()
  
  #Iterate subset list
  for (i in 1:length(subset_list)) {
    #Get i subset
    i_subset <- subset_list[[i]]
    i_name <- names(subset_list)[i]
    #Iterate each dataset within each subset
    for (j in 1:length(i_subset)) {
      #Get j i dataset
      j_i_data <- i_subset[[j]]
      j_i_name <- names(i_subset)[j]
      print(paste0("Analysing subset: ", names(subset_list)[i], " - ", "Dataset: ", names(i_subset)[j]))
      #Set categorical variable levels
      j_i_data@sam_data[[categorical_var[i]]] <- factor(j_i_data@sam_data[[categorical_var[i]]],
                                                        levels=sort(unique(j_i_data@sam_data[[categorical_var[i]]])))
      
      #Process kraken taxonomy
      if (length(grep("kraken", names(i_subset)[j], ignore.case=T)) != 0) {
        #Kraken stats
        kraken_significant <- differential_abundance_tests(j_i_data, categorical_var[i], p_threshold, padj_threshold,
                                                           run_deseq=T, run_microbiomemarker=T, run_maaslin=T, analyze_genes=F, j_i_name, i_name)
        significant <- kraken_significant
        #Get significant taxa (search for exacts taxa names at the beginning and the end of string)
        kraken_significant_no_special <- unique(gsub(" ", "", gsub("=", "", stringr:::str_replace_all(kraken_significant, "[[:punct:]]", ""))))
        rownames_no_special <- gsub(" ", "", gsub("=", "", stringr:::str_replace_all(rownames(j_i_data@otu_table), "[[:punct:]]", "")))
        matches <- grep(paste0("^", kraken_significant_no_special, "$", collapse="|"), rownames_no_special, ignore.case=F, value=F)
        print(paste0(names(i_subset)[j], " - Significant data: ", length(matches) == length(kraken_significant_no_special)))
        #Select matching rownames
        kraken_significant_rownames <- rownames(j_i_data@otu_table)[matches]
        #Post-processing of statistical results
        increments_and_taxa_filtering <- increments_and_genes_with_taxa(j_i_data, kraken_significant_rownames, all_clades, calculate_increments=T, increment_threshold, filter_genes=F, i_name)
        #Generate differential tables
        dif_tables <- differential_tables(j_i_data, kraken_significant_rownames, categorical_var[i])
        #Summarise results
        summary <- differential_tables_summary(j_i_data@sam_data, categorical_var[i], dif_tables, j_i_data, 
                                               paste0(names(subset_list)[i], "_", names(i_subset)[j]), 
                                               analyze_taxonomy=T, analyze_families=F, analyze_pathways=F)
      }
      
      #Process humann taxonomy
      if (length(grep("humann_taxa", names(i_subset)[j], ignore.case=T)) != 0) {
        #Humann stats
        humann_significant <- differential_abundance_tests(j_i_data, categorical_var[i], p_threshold, padj_threshold,
                                                           run_deseq=T, run_microbiomemarker=T, run_maaslin=T, analyze_genes=F, j_i_name, i_name)
        significant <- humann_significant
        #Get significant taxa (search for exacts taxa names at the beginning and the end of string)
        humann_significant_no_special <- unique(stringr:::str_replace_all(gsub(".*remove", "", sub("^X.*?_", "remove", humann_significant)), "[[:punct:]]", ""))
        rownames_no_special <- stringr:::str_replace_all(gsub(".*remove", "", sub("^.*?_", "remove", rownames(j_i_data@otu_table))), "[[:punct:]]", "")
        matches <- grep(paste0("^", humann_significant_no_special, "$", collapse="|"), rownames_no_special, ignore.case=F, value=F)
        print(paste0(names(i_subset)[j], " - Significant data: ", length(matches) == length(humann_significant_no_special)))
        #Select matching rownames
        humann_significant_rownames <- rownames(j_i_data@otu_table)[matches]
        #Post-processing of statistical results
        increments_and_taxa_filtering <- increments_and_genes_with_taxa(j_i_data, humann_significant_rownames, all_clades, calculate_increments=T, increment_threshold, filter_genes=F, i_name)
        #Generate differential tables
        dif_tables <- differential_tables(j_i_data, humann_significant_rownames, categorical_var[i])
        #Summarise results
        summary <- differential_tables_summary(j_i_data@sam_data, categorical_var[i], dif_tables, j_i_data, 
                                               paste0(names(subset_list)[i], "_", names(i_subset)[j]), 
                                               analyze_taxonomy=T, analyze_families=F, analyze_pathways=F)
      }
      
      #Process humann functional data
      if (length(grep("humann_func_", names(i_subset)[j], ignore.case=T)) != 0) {
        #Humann stats
        humann_significant <- differential_abundance_tests(j_i_data, categorical_var[i], p_threshold, padj_threshold,
                                                           run_deseq=T, run_microbiomemarker=F, run_maaslin=T, analyze_genes=T, j_i_name, i_name)
        significant <- humann_significant
        #Get significant taxa (search for exacts taxa names at the beginning and the end of string)
        #Remove special characters
        humann_significant_no_special <- unique(sort(gsub(" ", "", gsub("\\|", "", stringr:::str_replace_all(humann_significant, "[[:punct:]]", "")))))
        rownames_no_special <- gsub(" ", "", gsub("\\|", "", stringr:::str_replace_all(rownames(j_i_data@otu_table), "[[:punct:]]", "")))
        matches <- match(humann_significant_no_special, rownames_no_special)
        # matches <- grep(paste0("^", humann_significant_no_special, "$", collapse="|"), rownames_no_special, ignore.case=F, value=F)
        print(paste0(names(i_subset)[j], " - Significant data: ", length(matches) == length(humann_significant_no_special)))
        #Select matching rownames
        humann_significant_rownames <- rownames(j_i_data@otu_table)[matches]
        #Post-processing of statistical results
        increments_and_taxa_filtering <- increments_and_genes_with_taxa(j_i_data, humann_significant_rownames, all_clades, calculate_increments=T, increment_threshold, filter_genes=T, i_name)
        #Generate differential tables
        dif_tables <- differential_tables(increments_and_taxa_filtering$functional_data_with_taxa, humann_significant_rownames, categorical_var[i])
        #Summarise results
        if (length(grep("pathways", names(i_subset)[j])) != 0) {
          summary <- differential_tables_summary(increments_and_taxa_filtering$functional_data_with_taxa@sam_data, categorical_var[i], dif_tables, increments_and_taxa_filtering$functional_data_with_taxa, 
                                                 paste0(names(subset_list)[i], "_", names(i_subset)[j]), 
                                                 analyze_taxonomy=T, analyze_families=F, analyze_pathways=T)
        } else {
          summary <- differential_tables_summary(increments_and_taxa_filtering$functional_data_with_taxa@sam_data, categorical_var[i], dif_tables, increments_and_taxa_filtering$functional_data_with_taxa, 
                                                 paste0(names(subset_list)[i], "_", names(i_subset)[j]), 
                                                 analyze_taxonomy=F, analyze_families=T, analyze_pathways=F)
        }
      }
      
      #Generate a list of results
      i_results <- list(significant, increments_and_taxa_filtering, summary)
      names(i_results) <- c("microbiota_statistics", "increments_and_taxa_filtering", "dif_tables_summary")
      #Append i results to a single list
      results[[length(results)+1]] <- i_results
      names(results)[[length(results)]] <- paste0(names(subset_list)[i], "_", names(i_subset)[j])
      
    }
  }
  
  #Return output list
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
