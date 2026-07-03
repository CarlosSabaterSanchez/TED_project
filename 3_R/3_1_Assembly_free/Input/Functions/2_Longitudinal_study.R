
####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3

#Custom function to process longitudinal 16S taxonomic data - Less restrictive statistical criteria - trends
longitudinal_study_trends <- function(selected_pseq, substrate_label, health_label, time_label, p_value_threshold) {
  
  #Load libraries
  library(DESeq2)
  library(tidyverse)
  library(vegan)
  library(ALDEx2)
  library(ANCOMBC)
  library(microbiomeMarker)
  
  #Remove the negative control
  selected_pseq <- prune_samples(selected_pseq@sam_data[[substrate_label]] != "-", selected_pseq)
  #Get health status levels
  health <- unique(selected_pseq@sam_data[[health_label]])
  #Generate an empty list to store pseq objects for different health status
  health_pseqs <- list() 
  #Generate pseq objects for different health status
  for (i in 1:length(health)) {
    pseq_i <- prune_samples(selected_pseq@sam_data[[health_label]] == health[i], selected_pseq)
    health_pseqs <- c(health_pseqs, pseq_i)
    names(health_pseqs)[i] <- health[i]
  }
  
  #Get initial time points for different health status
  #Generate an empty list to store pseq objects for initial timepoints of different health status
  initial_health_pseqs <- list() 
  #Generate pseq objects with initial timepoints for different health status
  for (i in 1:length(health_pseqs)) {
    pseq_initial_i <- prune_samples(health_pseqs[[i]]@sam_data[[time_label]] == min(as.numeric(health_pseqs[[i]]@sam_data[[time_label]])),
                                    health_pseqs[[i]])
    #Append the results to the list
    initial_health_pseqs <- c(initial_health_pseqs, pseq_initial_i)
    names(initial_health_pseqs)[i] <- paste0(names(health_pseqs)[i], "_initial")
  }
  
  #Generate pseq objects for each substrate
  #Get substrate levels
  substrate <- unique(selected_pseq@sam_data[[substrate_label]])
  #Generate an empty list to store pseq objects for each substrate
  substrate_with_initial_health_pseqs <- list() 
  #Generate pseq objects for each substrate
  for (i in 1:length(health_pseqs)) {
    for (j in 1:length(substrate)) {
      pseq_health_substrate_i <- prune_samples(health_pseqs[[i]]@sam_data[[substrate_label]] == substrate[j],
                                               health_pseqs[[i]])
      #Append initial time point if it is nor present in the pseq subset
      if (length(grep(min(as.numeric(health_pseqs[[i]]@sam_data[[time_label]])), 
                      pseq_health_substrate_i@sam_data[[time_label]])) == 0) {
        pseq_health_substrate_i <- merge_phyloseq(initial_health_pseqs[[i]], pseq_health_substrate_i)
      }
      #Append the results to the list
      substrate_with_initial_health_pseqs <- c(substrate_with_initial_health_pseqs, pseq_health_substrate_i)
      names(substrate_with_initial_health_pseqs)[length(substrate_with_initial_health_pseqs)] <- paste0(names(health_pseqs[i]), "_",
                                                                                                        substrate[j])
    }
  }

  #Create empty list to store significant taxa
  significant_list <- list()
  #Calculate statistically significant differences during fermentation for each substrate
  for (i in 1:length(substrate_with_initial_health_pseqs)) {
    #Remove non-numeric characters from time values
    substrate_with_initial_health_pseqs[[i]]@sam_data[[time_label]] <- as.numeric(gsub("[^0-9.-]", "", 
                                                                                       substrate_with_initial_health_pseqs[[i]]@sam_data[[time_label]]))
    #Get time levels for each substrate
    times_i <- sort(unique(substrate_with_initial_health_pseqs[[i]]@sam_data[[time_label]]))
    
    #Iterate time values to calculate statistically significant differences
    for (j in 1:length(times_i)) {
      #Use minimum fermentation time as reference to calculate significant increments during fermentation
      if (times_i[j] != min(times_i)) {
        compare_group_i_j <- merge_phyloseq(prune_samples(substrate_with_initial_health_pseqs[[i]]@sam_data[[time_label]] == min(times_i),
                                                          substrate_with_initial_health_pseqs[[i]]),
                                            prune_samples(substrate_with_initial_health_pseqs[[i]]@sam_data[[time_label]] == times_i[j],
                                                          substrate_with_initial_health_pseqs[[i]]))
        
        #Run ancom
        ancom_i_j <- ancombc(phyloseq = compare_group_i_j, formula = time_label,
                             p_adj_method = "holm", lib_cut = 0, #zero_cut = 0.90, 
                             group = time_label, struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                             max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
        
        #Get significant_taxa
        #Select taxa with significant p-value 
        df_filt <- ancom_i_j$res$p_val[which(ancom_i_j$res$p_val[,3] < p_value_threshold),]
        #######################################################################################
        #Microbiomemarker
        #Format input pseq object
        stats_input <- compare_group_i_j
        # stats_input <- pseq_div
        #Format input pseq object (optional?)
        # colnames(stats_input@tax_table) <- c("Species")
        #Normalize data prior statistical analysis
        # stats_input <- normalize(stats_input, method="TSS")
        #Run multiple statistical tests
        # aldex <- run_aldex(stats_input, group=time_label)
        ancom <- run_ancom(stats_input, group=time_label)
        # ancombc <- run_ancombc(stats_input, group=time_label)
        lefse <- try(run_lefse(stats_input, group=time_label), silent=T)
        if (class(lefse) == "try-error") {
          setClass("results", representation(name="character", marker_table="data.frame"))
          lefse <- new("results", name="", marker_table=data.frame())
        }
        metagenomeseq <- run_metagenomeseq(stats_input, group=time_label)
        #Select significant taxa
        marker_taxa <- unique(sort(c(
          # aldex@marker_table$feature,
          ancom@marker_table$feature,
          # ancombc@marker_table$feature,
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
        #Combine the results of all statistical methods
        marker_intersect <- sort(unique(c(marker_taxa_format, df_filt$taxon)))
        # marker_intersect <- sort(unique(c(marker_taxa_format, 
        #                                   gsub(".*remove", "", sub("_", "remove", df_filt$taxon)))))
        
        #######################################################################################
        
        #Select taxa higher after fermentation
        #Continue the analysis if there are any significant values
        if (length(marker_intersect) > 0) {
          for(k in 1:length(marker_intersect)) {
            #Get k taxa
            k_taxa <- marker_intersect[k]
            
            #Add humann numbers to taxa
            k_taxa <- rownames(compare_group_i_j@otu_table)[grep(paste0(k_taxa, "$"), rownames(compare_group_i_j@otu_table))]
            
            #Calculate taxa counts at the beginning of the fermentation
            count_zero <- prune_samples(compare_group_i_j@sam_data[[time_label]] != times_i[j], 
                                        compare_group_i_j)@otu_table[k_taxa,]
            
            #Calculate taxa counts at the end of the fermentation
            count_ferm <- prune_samples(compare_group_i_j@sam_data[[time_label]] == times_i[j],
                                        compare_group_i_j)@otu_table[k_taxa,]
            
            #Calculate the increment in taxa counts per individual
            count_dif <- count_ferm - count_zero
            
            #Print taxa higher after fermentation in all individuals (no zero nor negative increments) for those taxa higher than 0.1% abundance
            if(length(which(count_dif < 0)) == 0 & round(mean(count_ferm),1) >= min(colSums(compare_group_i_j@otu_table))/1000) {
              print(c(k_taxa, "Initial:", round(mean(count_zero)), "Final:", round(mean((count_ferm)))))
              some_indiv <- which(count_ferm > 0)
              increment <- round(mean(count_dif[count_dif != 0]),2)
              sd_increment <- round(sd(count_dif[count_dif != 0]),2)
              significant_list <- c(significant_list, paste0(k_taxa, " + Time: ", times_i[j],
                                                             " + Individuals: ", 
                                                             paste(some_indiv, collapse = " and "), 
                                                             paste0(" + Increment: ", increment),
                                                             paste0(" + SD: ", sd_increment)))
              names(significant_list)[length(significant_list)] <- paste0(names(substrate_with_initial_health_pseqs[i]), 
                                                                          "_Time_", times_i[j])
            }
          }
        }
      }
    }
  }
  
  #Combine list elements from the same sample (i.e. same health status, substrate and time)
  significant_list <- with(stack(significant_list), split(values, ind))
  #Replace NaN in bacterial increment (significant taxa showing 0 rounded increments)
  significant_list <- lapply(X = significant_list, FUN = function(t) gsub(pattern = "NaN", replacement = "0", x = t, fixed = TRUE))
  #Get all significant taxa
  significant_taxa <- sort(unique(gsub(" \\+ .*", "", unlist(significant_list))))
  #Find the total number of significant taxa
  n_taxa <- length(significant_taxa)
  #Iterate list elements to add NAs until maximum length is reached
  for (i in 1:length(significant_list)) {
    #Determine the number of NAs to add until maximum length is reached
    i_significant_list <- significant_list[[i]]
    if (length(which(duplicated(i_significant_list))) != 0) {
      i_significant_list <- i_significant_list[-which(duplicated(i_significant_list))]
    }
    nas_to_add <- n_taxa - length(i_significant_list)
    #Add NAs until maximum length is reached
    significant_list[[i]] <- c(significant_list[[i]], rep(NA, nas_to_add))
    #Reorder list elements
    genera_no_indiv <- gsub(" \\+ .*", "", significant_list[[i]])
    reorder_idx <- match(significant_taxa, genera_no_indiv)
    significant_list[[i]] <- significant_list[[i]][reorder_idx]
  }
  
  #Combine significant taxa found in each sample into a single data frame
  significant_all <- data.frame(do.call(rbind, significant_list))
  colnames(significant_all) <- significant_taxa
  #Add substrate and time metadata to the newly generated data frame
  #Time
  significant_all <- tibble::add_column(significant_all, Time = gsub(".*Time_", "", rownames(significant_all)), 
                                        .before = 1)
  #Substrate
  significant_all <- tibble::add_column(significant_all, Substrate = sub(".*?_", "", gsub("_Time_.*", "", rownames(significant_all))), 
                                        .before = 1)
  #Health status
  significant_all <- tibble::add_column(significant_all, Status = gsub("_.*", "", rownames(significant_all)), 
                                        .before = 1)
  
  #Stack taxonomy results per substrate (combine data from all three points)
  significant_all_per_substrate <- significant_all %>% group_by(Substrate, Status) %>%
    summarise_each(funs(paste(., collapse = ";")))
  
  #Select the  time point where the maximum significant increase of these taxa is observed
  for (i in 1:nrow(significant_all_per_substrate)) {
    for (j in 1:ncol(significant_all_per_substrate)) {
      if (str_count(significant_all_per_substrate[i, j], "Increment") >= 1) {
        #Determine samples showing significant increments
        i_showing_increments <- significant_all_per_substrate[i, j]
        #Split samples showing significant increments
        i_splits <- unlist(str_split(i_showing_increments, "(?<=[:alnum:];)"))
        #Select non-NA samples
        i_splits <- i_splits[grep("Increment", i_splits)]
        #Convert increments to numeric
        i_increments <- as.numeric(gsub(" +.*", "", 
                                        gsub("-", "", gsub(".*Increment: |;.*", "", i_splits))))
        #Select maximum increments
        i_max_increment <- which(i_increments == max(i_increments))
        #If there are two or more samples showing the maximum increment values, select the first time
        if (length(i_max_increment) > 1) {
          i_max_increment <- i_max_increment[1]
        }
        #Get maximum value in the cell
        i_selected <- i_splits[i_max_increment]
        #Clean cell format
        i_selected <- sub(";$", "", i_selected)
        #Select cell values
        significant_all_per_substrate[i, j] <- i_selected  
      }
    }
  }
  
  #Replace "NA-NA-NA*" by "-"
  significant_all_per_substrate <- map_df(significant_all_per_substrate, 
                                          ~ gsub("NA.*", "-", .x))
  #Format health status table
  significant_all_per_substrate$Status <- gsub(";.*", "", significant_all_per_substrate$Status)
  #Remove time column (this information is already included in the cells)
  significant_all_per_substrate <- subset(significant_all_per_substrate, select=-c(Time))
  
  #Return results table and pseq objects (subsetted by health status and substrate) as list
  results <- list(significant_all_per_substrate, substrate_with_initial_health_pseqs)
  names(results) <- c("significant_all_per_substrate", "substrate_with_initial_health_pseqs")
  return(results)
  
}

####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3

#Custom function to process longitudinal 16S taxonomic data
longitudinal_heatmap <- function(selected_pseq, selected_longitudinal, time_label, substrate_label, health_label) {
  
  #Load libraries
  library(reshape2)
  library(tidyverse)
  
  #Remove the negative control
  selected_pseq <- prune_samples(selected_pseq@sam_data[[substrate_label]] != "-", selected_pseq)
  #Transform OTU abundances to percentages
  pseq_percentage <- microbiome::transform(selected_pseq, "compositional") 
  #Multiply per 100
  pseq_percentage@otu_table <- pseq_percentage@otu_table*100
  #Remove non-numeric characters from time values
  pseq_percentage@sam_data[[time_label]] <- as.numeric(gsub("[^0-9.-]", "", 
                                                            pseq_percentage@sam_data[[time_label]]))
  
  #Calculate mean taxa abundance at the end of the fermentation
  for (i in 1:nrow(selected_longitudinal)) {
    for (j in 3:ncol(selected_longitudinal)) {
      if (selected_longitudinal[i,j] != "-") {
        #Get taxa
        i_taxa <- colnames(selected_longitudinal)[j] 
        #Get fermentation time
        i_time <- as.numeric(gsub(" .*", "", gsub(".*Time: ", "", selected_longitudinal[i,j])))
        #Get individuals where changes are observed
        i_individuals <- gsub(" \\+.*", "", gsub(".*Individuals: ", "", selected_longitudinal[i,j]))
        i_individuals <- as.numeric(strsplit(i_individuals, " and ")[[1]])
        #Get status
        i_status <- selected_longitudinal$Status[i]
        #Get substrate
        i_substrate <- selected_longitudinal$Substrate[i]
        #Calculate taxa percentages at the end of fermentation
        i_pseq <- prune_samples(pseq_percentage@sam_data[[health_label]] == i_status &
                                  pseq_percentage@sam_data[[substrate_label]] == i_substrate &
                                  pseq_percentage@sam_data[[time_label]] == i_time,
                                pseq_percentage)
        #Get OTU table for the selected taxa
        i_otu <- as.data.frame(i_pseq@otu_table[i_taxa,])
        #Calculate percentage (%) of each taxa at the end of fermentation
        i_percentage_end <- round(mean(i_otu[i_otu != 0]), 2)
        i_sd_percentage_end <- round(sd(i_otu[i_otu != 0]), 2)
        #Add percentage data to each individual cell
        selected_longitudinal[i,j] <- paste0(selected_longitudinal[i,j], " + Percentage: ", i_percentage_end, " + SD: ", i_sd_percentage_end)
      }
    }
  }
  
  #Increment data frame
  increment_counts <- selected_longitudinal %>% mutate_all(funs(str_replace(., ".*Increment: ", "")))
  increment_counts <- increment_counts %>% mutate_all(funs(str_replace(., "\\+.*", "")))
  increment_counts[,3:ncol(increment_counts)] <- increment_counts[,3:ncol(increment_counts)] %>% mutate_all(funs(str_replace(., "-", "0")))
  #Set dataframe names
  rownames(increment_counts) <- paste0(#increment_counts$Status, " + ", 
    gsub("_", " ", increment_counts$Substrate))
  names(increment_counts) <- gsub("_", " ", names(increment_counts))
  
  #Get heatmap labels
  counts_labels <- selected_longitudinal %>% mutate_all(funs(str_replace_all(., ".*Time: ", "T:")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "Individuals: ", "I:")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "Increment: ", "Incr:")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "Percentage: ", "P:")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., " and ", ",")))
  #Remove standard deviations
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., " \\+ SD\\:[:space:].*[:space:]\\+[:space:]P", " \\+ P")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., " \\+ SD\\:[:space:].*", "")))
  #Remove separators (replace + by ;)
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "[:space:]\\+[:space:]", "; ")))
  #Remove - from empty cells
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "-", "")))
  #Set carriage return
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "P:", "\nP:")))
  
  #Set dataframe names
  names(counts_labels) <- sub(".*? ", "", names(increment_counts)) #Remove humann numbers
  rownames(counts_labels) <- rownames(increment_counts)
  
  #COUNTS HEATMAP
  #Prepare heatmap data
  heatmap_data <- t(increment_counts)[3:ncol(increment_counts),]
  heatmap_data <- mutate_all(as.data.frame(heatmap_data), function(x) as.numeric(as.character(x)))
  #Prepare heatmap labels
  heatmap_labels <- t(counts_labels)[3:ncol(counts_labels),]
  
  #Remove special characters from microbial taxa in order to generate clean plots
  rownames(heatmap_data) <- str_replace_all(rownames(heatmap_data), "[^[:alnum:]- ]", "")
  rownames(heatmap_labels) <- str_replace_all(rownames(heatmap_labels), "[^[:alnum:]- ]", "")
  
  #Re-order alphabetically bacterial clade names
  heatmap_data <- heatmap_data[sort(rownames(heatmap_data)),]
  heatmap_labels <- heatmap_labels[sort(rownames(heatmap_labels)),]
  
  #Create a new vector to store maximum values per taxa
  max_increments_per_taxa <- c()
  #Determine different ranges of increments in bacterial values (for high and low-abundant genera)
  for (i in 1:nrow(heatmap_data)) {
    max_increments_per_taxa <- c(max_increments_per_taxa, max(heatmap_data[i,]))
  }
  
  #Divide bacterial value ranges into three different scales using Q50%, Q75% and Q100%
  scale_low <- c(0 , quantile(max_increments_per_taxa)[3])
  scale_medium <- c(quantile(max_increments_per_taxa)[3] , quantile(max_increments_per_taxa)[4])
  scale_high <- c(quantile(max_increments_per_taxa)[4] , quantile(max_increments_per_taxa)[5])
  
  #Generate new dataframes containing low, medium and high increments
  #Low increments
  data_low <- heatmap_data[which(max_increments_per_taxa <= max(scale_low)),]
  rownames(data_low) <- sub(".*? ", "", rownames(data_low)) #Remove humann numbers
  labels_low <- heatmap_labels[rownames(data_low),]
  #Medium increments
  data_medium <- heatmap_data[which(max_increments_per_taxa > max(scale_low) & max_increments_per_taxa <= max(scale_medium)),]
  rownames(data_medium) <- sub(".*? ", "", rownames(data_medium)) #Remove humann numbers
  labels_medium <- heatmap_labels[rownames(data_medium),]
  #High increments
  data_high <- heatmap_data[which(max_increments_per_taxa > max(scale_medium) & max_increments_per_taxa <= max(scale_high)),]
  rownames(data_high) <- sub(".*? ", "", rownames(data_high)) #Remove humann numbers
  labels_high <- heatmap_labels[rownames(data_high),]
  
  #Format data for heatmap plots
  #Low increments
  low_format <- melt(t(data_low), id = colnames(t(data_low)))
  low_labels_format <- melt(t(labels_low), id = colnames(t(labels_low)))
  #Remove increments from label (will be added to the plot separately)
  low_labels_format <- low_labels_format %>% mutate_all(funs(str_replace_all(., " Incr:.*;", "")))
  
  #Medium increments
  medium_format <- melt(t(data_medium), id = colnames(t(data_medium)))
  medium_labels_format <- melt(t(labels_medium), id = colnames(t(labels_medium)))
  #Remove increments from label (will be added to the plot separately)
  medium_labels_format <- medium_labels_format %>% mutate_all(funs(str_replace_all(., " Incr:.*;", "")))
  
  #High increments
  high_format <- melt(t(data_high), id = colnames(t(data_high)))
  high_labels_format <- melt(t(labels_high), id = colnames(t(labels_high)))
  #Remove increments from label (will be added to the plot separately)
  high_labels_format <- high_labels_format %>% mutate_all(funs(str_replace_all(., " Incr:.*;", "")))
  
  
  
  
  #Plot heatmap
  #Low increments
  low_plot <- ggplot(data = low_format, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(aes(fill = value)) + scale_y_discrete(limits=rev) +
    scale_fill_gradient2(name = "Low significant increments in read counts", 
                         low = "white", high = "green", breaks = c(floor(min(low_format$value[low_format$value > 0])),
                                                                   floor(max(low_format$value)/2), 
                                                                   floor(max(low_format$value)))) +  
    labs(x = "") +
    geom_tile(colour = "gray50") +
    theme(text=element_text(size=20), 
          axis.title.y = element_blank(),
          legend.position = "top", legend.justification="center",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"),
          legend.title = element_text(size=13), legend.text = element_text(size=10)) +
    geom_text(aes(label = gsub("\\(\\+0\\)", "", paste0(low_labels_format$value, " \n(+", round(low_format$value), ")"))
    ), size=4.0)
  
  
  
  
  #Medium increments
  medium_plot <- ggplot(data = medium_format, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(aes(fill = value)) + scale_y_discrete(limits=rev) +
    scale_fill_gradient2(name = "Moderate significant increments in read counts", 
                         low = "white", high = "yellow", breaks = c(floor(min(medium_format$value[medium_format$value > 0])),
                                                                    floor(max(medium_format$value)/2), 
                                                                    floor(max(medium_format$value)))) + 
    labs(x = "") +
    geom_tile(colour = "gray50") +
    theme(text=element_text(size=20), 
          axis.title.y = element_blank(),
          legend.position = "top", legend.justification="center",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"),
          legend.title = element_text(size=13), legend.text = element_text(size=10)) +
    geom_text(aes(label = gsub("\\(\\+0\\)", "", paste0(medium_labels_format$value, " \n(+", round(medium_format$value), ")"))), size=4.0)
  
  
  
  
  #High increments
  high_plot <- ggplot(data = high_format, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(aes(fill = value)) + scale_y_discrete(limits=rev) +
    scale_fill_gradient2(name = "High significant increments in read counts", 
                         low = "white", high = "red", breaks = c(floor(min(high_format$value[high_format$value > 0])),
                                                                 floor(max(high_format$value)/2), 
                                                                 floor(max(high_format$value)))) + 
    labs(x = "") +
    geom_tile(colour = "gray50") +
    theme(text=element_text(size=20), 
          axis.title.y = element_blank(),
          legend.position = "top", legend.justification="center",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"), 
          legend.title = element_text(size=13), legend.text = element_text(size=10)) +
    geom_text(aes(label = gsub("\\(\\+0\\)", "", paste0(high_labels_format$value, " \n(+", round(high_format$value), ")"))), size=4.0)
  
  
  
  
  #Plot multiple heatmaps simultaneously
  dev.new(width=5, height=4)
  gridExtra::grid.arrange(low_plot)
  low_plot_record <- recordPlot()
  dev.new(width=5, height=4)
  gridExtra::grid.arrange(medium_plot, high_plot, nrow=1, widths=c(5, 5))
  high_medium_plot_record <- recordPlot()
  
  #Save heatmap plots
  #g <- gridExtra::arrangeGrob(low_plot) #Generates g
  #ggsave(file="heatmaps.png", g, device='png', dpi=700) #Saves g
  
  #g <- gridExtra::arrangeGrob(medium_plot, high_plot, nrow=1, widths=c(5, 5)) #Generates g
  #ggsave(file="heatmaps.png", g, device='png', dpi=700) #Saves g
  
  #Return longitudinal heatmaps as list
  heatmap_plots <- list(low_plot_record, high_medium_plot_record)
  names(heatmap_plots) <- c("low_plot_record", "high_medium_plot_record")
  return(heatmap_plots)
  
}

####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3

#Custom function to process longitudinal 16S taxonomic data
longitudinal_single_heatmap <- function(selected_pseq, selected_longitudinal, time_label, substrate_label, health_label) {
  
  #Load libraries
  library(reshape2)
  library(tidyverse)
  
  #Remove the negative control
  selected_pseq <- prune_samples(selected_pseq@sam_data[[substrate_label]] != "-", selected_pseq)
  #Transform OTU abundances to percentages
  pseq_percentage <- microbiome::transform(selected_pseq, "compositional") 
  #Multiply per 100
  pseq_percentage@otu_table <- pseq_percentage@otu_table*100
  colSums(pseq_percentage@otu_table)
  #Remove non-numeric characters from time values
  pseq_percentage@sam_data[[time_label]] <- as.numeric(gsub("[^0-9.-]", "", 
                                                            pseq_percentage@sam_data[[time_label]]))
  
  #Calculate mean taxa abundance at the end of the fermentation
  for (i in 1:nrow(selected_longitudinal)) {
    for (j in 3:ncol(selected_longitudinal)) {
      if (selected_longitudinal[i,j] != "-") {
        #Get taxa
        i_taxa <- colnames(selected_longitudinal)[j] 
        #Get fermentation time
        i_time <- as.numeric(gsub(" .*", "", gsub(".*Time: ", "", selected_longitudinal[i,j])))
        #Get individuals where changes are observed
        i_individuals <- gsub(" \\+.*", "", gsub(".*Individuals: ", "", selected_longitudinal[i,j]))
        i_individuals <- as.numeric(strsplit(i_individuals, " and ")[[1]])
        #Get status
        i_status <- selected_longitudinal$Status[i]
        #Get substrate
        i_substrate <- selected_longitudinal$Substrate[i]
        #Calculate taxa percentages at the end of fermentation
        i_pseq <- prune_samples(pseq_percentage@sam_data[[health_label]] == i_status &
                                  pseq_percentage@sam_data[[substrate_label]] == i_substrate &
                                  pseq_percentage@sam_data[[time_label]] == i_time,
                                pseq_percentage)
        #Get OTU table for the selected taxa
        i_otu <- as.data.frame(i_pseq@otu_table[i_taxa,])
        #Calculate percentage (%) of each taxa at the end of fermentation
        i_percentage_end <- round(mean(i_otu[i_otu != 0]), 2)
        i_sd_percentage_end <- round(sd(i_otu[i_otu != 0]), 2)
        #Add percentage data to each individual cell
        selected_longitudinal[i,j] <- paste0(selected_longitudinal[i,j], " + Percentage: ", i_percentage_end, " + SD: ", i_sd_percentage_end)
      }
    }
  }
  
  #Increment data frame
  increment_counts <- selected_longitudinal %>% mutate_all(funs(str_replace(., ".*Increment: ", "")))
  increment_counts <- increment_counts %>% mutate_all(funs(str_replace(., "\\+.*", "")))
  increment_counts[,3:ncol(increment_counts)] <- increment_counts[,3:ncol(increment_counts)] %>% mutate_all(funs(str_replace(., "-", "0")))
  #Set dataframe names
  rownames(increment_counts) <- paste0(#increment_counts$Status, " + ", 
    gsub("_", " ", increment_counts$Substrate))
  names(increment_counts) <- gsub("_", " ", names(increment_counts))
  
  #Get heatmap labels
  counts_labels <- selected_longitudinal %>% mutate_all(funs(str_replace_all(., ".*Time: ", "T:")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "Individuals: ", "I:")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "Increment: ", "Incr:")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "Percentage: ", "P:")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., " and ", ",")))
  #Remove standard deviations
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., " \\+ SD\\:[:space:].*[:space:]\\+[:space:]P", " \\+ P")))
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., " \\+ SD\\:[:space:].*", "")))
  #Remove separators (replace + by ;)
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "[:space:]\\+[:space:]", "; ")))
  #Remove - from empty cells
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "-", "")))
  #Set carriage return
  counts_labels <- counts_labels %>% mutate_all(funs(str_replace_all(., "P:", "\nP:")))
  
  
  #Set dataframe names
  names(counts_labels) <- names(increment_counts)
  rownames(counts_labels) <- rownames(increment_counts)
  
  #COUNTS HEATMAP
  #Prepare heatmap data
  heatmap_data <- t(increment_counts)[3:ncol(increment_counts),]
  heatmap_data <- mutate_all(as.data.frame(heatmap_data), function(x) as.numeric(as.character(x)))
  #Prepare heatmap labels
  heatmap_labels <- t(counts_labels)[3:ncol(counts_labels),]
  
  #Remove special characters from microbial taxa in order to generate clean plots
  rownames(heatmap_data) <- str_replace_all(rownames(heatmap_data), "[^[:alnum:]- ]", "")
  rownames(heatmap_labels) <- str_replace_all(rownames(heatmap_labels), "[^[:alnum:]- ]", "")
  
  #Re-order alphabetically bacterial clade names
  # heatmap_data <- heatmap_data[sort(rownames(heatmap_data)),]
  # heatmap_labels <- heatmap_labels[sort(rownames(heatmap_labels)),]
  
  #Create a new vector to store maximum values per taxa
  max_increments_per_taxa <- c()
  #Determine different ranges of increments in bacterial values (for high and low-abundant genera)
  for (i in 1:nrow(heatmap_data)) {
    max_increments_per_taxa <- c(max_increments_per_taxa, max(heatmap_data[i,]))
  }
  
  #Generate new dataframes containing low, medium and high increments
  #Low increments
  # rownames(heatmap_data) <- sub(".*? ", "", rownames(heatmap_data)) #Remove humann numbers
  data <- heatmap_data
  labels <- heatmap_labels
  
  #Format data for heatmap plots
  data_format <- melt(t(data), id = colnames(t(data)))
  data_labels_format <- melt(t(labels), id = colnames(t(labels)))
  #Remove increments from label (will be added to the plot separately)
  data_labels_format <- data_labels_format %>% mutate_all(funs(str_replace_all(., " Incr:.*;", "")))
  
  
  
  
  #Plot heatmap
  #Increments
  plot <- ggplot(data = data_format, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(aes(fill = value)) + scale_y_discrete(limits=rev) +
    scale_fill_gradient2(name = "Significant increments \nin read counts", 
                         low = "green", high = "darkgreen") +  
    labs(x = "") +
    geom_tile(colour = "gray50") +
    theme(text=element_text(size=25), 
          axis.title.y = element_blank(),
          legend.position = "top", legend.justification="center",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"),
          legend.title = element_text(size=15), legend.text = element_text(size=15)) #+
    # geom_text(aes(label = gsub(" \\(\\+0\\)", "", paste0(data_labels_format$value, " (+", round(data_format$value), ")"))), size=5.0)
  
  
  
  
  
  
  #Plot multiple heatmaps simultaneously
  dev.new(width=5, height=4)
  gridExtra::grid.arrange(plot)
  low_plot_record <- recordPlot()
  
  #Save heatmap plots
  #g <- gridExtra::arrangeGrob(low_plot) #Generates g
  #ggsave(file="heatmaps.png", g, device='png', dpi=700) #Saves g
  
  #g <- gridExtra::arrangeGrob(medium_plot, high_plot, nrow=1, widths=c(5, 5)) #Generates g
  #ggsave(file="heatmaps.png", g, device='png', dpi=700) #Saves g
  
  #Return longitudinal heatmaps as list
  heatmap_plots <- list(plot)
  names(heatmap_plots) <- c("plot")
  return(heatmap_plots)
  
}

####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3
