
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to calculate statistically significant differences
statistical_tests <- function(classif_meta, classif_group) {
  
  #Combine labels and numeric data
  classif_data <- cbind.data.frame(classif_meta, classif_group)
  names(classif_data)[1] <- c("Group")
  
  #1- Test normality (Kolmogorov-Smirnov for higher than 50 observations, otherwise, Shapiro-Wilk)
  #Normal distribution: p-value above 0.05 - visualize with hist() function
  #Non-normal distribution: p-value below 0.05 - visualize with hist() function
  shapiro_value_df <- data.frame(shapiro_value= numeric(0), stringsAsFactors = FALSE)
  for (i in (1+ncol(classif_meta)):ncol(classif_data)){
    shapiro <- shapiro.test(classif_data[,i])
    shapiro_p_value <- shapiro$p.value
    shapiro_value_df <- rbind(shapiro_value_df, shapiro_p_value[1])
  }
  rownames(shapiro_value_df) <- names(classif_group)
  names(shapiro_value_df) <- c("p_value")
  
  #Select non-normal and normal data
  non_normal <- data.frame(metric = numeric(0), feature = numeric(0), stringsAsFactors = FALSE)
  normal <- data.frame(metric = numeric(0), feature = numeric(0), stringsAsFactors = FALSE)
  for (i in ncol(classif_meta):ncol(classif_data)){
    if (is.na(shapiro_value_df[i,]) == F) {
      if(shapiro_value_df[i,] < 0.05) {
        p_value_and_name <- data.frame(shapiro_value_df[i,], rownames(shapiro_value_df)[i])
        non_normal <- rbind(non_normal, p_value_and_name)
      } else {
        p_value_and_name <- data.frame(shapiro_value_df[i,], rownames(shapiro_value_df)[i])
        normal <- rbind(normal, p_value_and_name)
      }
    }
  }
  #Set names
  names(non_normal) <- c("metric", "feature")
  names(normal) <- c("metric", "feature")
  
  
  
  
  #2- Check variance homogeneity (needed for parametric distribution)
  #F/Fisher test specific for 2 groups, Bartlett test for more groups
  #Levent test for every number of groups
  #p-value (Pr(>F)) above 0.05 - homogeneity
  #p-value (Pr(>F)) below 0.05 - no homogeneity
  levene_value_df <- data.frame(levene_value= numeric(0), stringsAsFactors = FALSE)
  for (i in (1+ncol(classif_meta)):ncol(classif_data)){
    levene <- leveneTest(classif_data[,i], classif_data$Group)
    levene_p_value <- levene$`Pr(>F)`
    levene_value_df <- rbind(levene_value_df, levene_p_value[1])
  }
  rownames(levene_value_df) <- names(classif_group)
  names(levene_value_df) <- c("p_value")
  
  #Select non-homogenous and homogenous data
  no_homogeneity <- data.frame(metric = numeric(0), feature = numeric(0), stringsAsFactors = FALSE)
  homogeneity <- data.frame(metric = numeric(0), feature = numeric(0), stringsAsFactors = FALSE)
  for (i in ncol(classif_meta):ncol(classif_data)){
    if (is.na(levene_value_df[i,]) == F) {
      if(levene_value_df[i,] < 0.05) {
        p_value_and_name <- data.frame(levene_value_df[i,], rownames(levene_value_df)[i])
        no_homogeneity <- rbind(no_homogeneity, p_value_and_name)
      } else {
        p_value_and_name <- data.frame(levene_value_df[i,], rownames(levene_value_df)[i])
        homogeneity <- rbind(homogeneity, p_value_and_name)
      }
    }
  }
  #Set names
  names(no_homogeneity) <- c("metric", "feature")
  names(homogeneity) <- c("metric", "feature") #Select parametric data
  
  
  
  
  #3- Select parametric and non-parametric criteria
  parametric <- intersect(normal[,2], homogeneity[,2])
  parametric <- as.data.frame(parametric)
  names(parametric) <- c("feature")
  #Select parametric and non-parametric data
  non_parametric_data <- classif_data[,! names(classif_data) %in% parametric$feature, drop = F]
  parametric_data <- classif_data[,names(classif_data) %in% parametric$feature, drop = F]
  parametric_data <- cbind(classif_data[1], parametric_data)
  
  
  
  
  #4- NON-PARAMETRIC DATA
  #Mann-Whitney U: for two groups
  #Get p-values: Kruskal Wallis by rank test for more than 2 groups
  kruskal_value_df <- data.frame(p_value= numeric(0), stringsAsFactors = FALSE)
  #If there are any non-parametric variables, the number of culumns should be higher than 1
  if (ncol(non_parametric_data) > 1) {
    for (i in (1+ncol(classif_meta)):ncol(non_parametric_data)){
      kruskal <- kruskal.test(non_parametric_data[,i] ~ non_parametric_data$Group)
      p_value <- kruskal$p.value
      kruskal_value_df <- rbind(kruskal_value_df, p_value)
    }
    rownames(kruskal_value_df) <- names(non_parametric_data)[2:i]
    names(kruskal_value_df) <- c("p_value")
  }
  
  #Calculate q-values
  q_value <- p.adjust(kruskal_value_df[,1], method = "fdr")
  
  #Select significant and non-significant features
  # metric <- as.data.frame(q_value)
  metric <- kruskal_value_df
  significant_features <- data.frame(metric = numeric(0), feature = numeric(0), stringsAsFactors = FALSE)
  non_significant_features <- data.frame(metric = numeric(0), feature = numeric(0), stringsAsFactors = FALSE)
  #Iterate features
  for (i in (ncol(classif_meta)):ncol(non_parametric_data)){
    if (is.na(metric[i,]) == F) {
      if(metric[i,] < 0.05) {
        metric_and_name <- data.frame(metric[i,], rownames(metric)[i])
        significant_features <- rbind(significant_features, metric_and_name)
      } else {
        metric_and_name <- data.frame(metric[i,], rownames(metric)[i])
        non_significant_features <- rbind(non_significant_features, metric_and_name)
      }
    }
  }
  #Set names
  names(significant_features) <- c("metric", "feature")
  names(non_significant_features) <- c("metric", "feature")
  #Generate new objects
  no_parametric_significant_features <- significant_features
  no_parametric_non_significant_features <- non_significant_features
  
  
  
  
  #5- PARAMETRIC DATA
  #Get p-values: t-test only for two groups, ANOVA for multiple groups
  anova_value_df <- data.frame(p_value= numeric(0), stringsAsFactors = FALSE)
  #If there are any parametric variables, the number of culumns should be higher than 1
  if (ncol(parametric_data) > 1) {
    for (i in (1+ncol(classif_meta)):ncol(parametric_data)){
      anova <- aov(parametric_data[,i] ~ parametric_data$Group, data = parametric_data)
      p_value <- summary(anova)[[1]][["Pr(>F)"]][1]
      anova_value_df <- rbind(anova_value_df, p_value)
    }
    rownames(anova_value_df) <- names(parametric_data)[2:i]
    names(anova_value_df) <- c("p_value")
  }
  
  #Calculate q-values
  q_value <- p.adjust(anova_value_df[,1], method = "fdr")
  
  #Select significant and non-significant features
  # metric <- as.data.frame(q_value)
  metric <- anova_value_df
  significant_features <- data.frame(metric = numeric(0), feature = numeric(0), stringsAsFactors = FALSE)
  non_significant_features <- data.frame(metric = numeric(0), feature = numeric(0), stringsAsFactors = FALSE)
  
  #Iterate features
  for (i in (ncol(classif_meta)):ncol(non_parametric_data)){
    if (is.na(metric[i,]) == F) {
      if(metric[i,] < 0.05) {
        metric_and_name <- data.frame(metric[i,], rownames(metric)[i])
        significant_features <- rbind(significant_features, metric_and_name)
      } else {
        metric_and_name <- data.frame(metric[i,], rownames(metric)[i])
        non_significant_features <- rbind(non_significant_features, metric_and_name)
      }
    }
  }
  #Set names
  names(significant_features) <- c("metric", "feature")
  names(non_significant_features) <- c("metric", "feature")
  #Generate new objects
  parametric_significant_features <- significant_features
  parametric_non_significant_features <- non_significant_features
  
  
  
  
  #Export results as list
  #Important: feature number refers to feature number (other than Group) in parametric or 
  #non-parametric datasets, not the original metadata file
  results <- list(normal, non_normal, homogeneity, no_homogeneity, 
                  parametric_data, non_parametric_data,
                  no_parametric_significant_features, no_parametric_non_significant_features,
                  parametric_significant_features, parametric_non_significant_features)
  
  #Set list names
  names(results) <- c("normal", "non_normal", "homogeneity", "no_homogeneity", 
                                  "parametric_data", "non_parametric_data",
                                  "no_parametric_significant_features", "no_parametric_non_significant_features",
                                  "parametric_significant_features", "parametric_non_significant_features")
  
  #Export list
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to perform post-hoc tests
post_hoc <- function(sig_param_data, sig_no_param_data) {
  
  #Parametric features
  #Perform post-hoc test if there are any significant features
  if (length(sig_param_data) > 0) {
    #Create a new list to store the results
    param_group_list <- list()
    #Iterate parametric dataset
    for (i in 2:ncol(sig_param_data)) {
      #Run ANOVA
      i_aov <- aov(lm(sig_param_data[[names(sig_param_data)[i]]] ~ sig_param_data[["Group"]], 
                      sig_param_data))
      #Run tukey test
      i_tukey <- HSD.test(i_aov, 'sig_param_data[["Group"]]', group=TRUE)
      
      #Get groups
      i_tukey_means <- i_tukey$means
      i_tukey_groups <- i_tukey$groups
      names(i_tukey_means)[1] <- names(sig_param_data)[i]
      names(i_tukey_groups)[1] <- names(sig_param_data)[i]
      
      #Round mean values
      i_tukey_groups[1] <- format(round(i_tukey_groups[1],2), nsmall=2)
      #Add standard deviations
      i_tukey_groups$SD <- format(round(i_tukey_means$std,2), nsmall=2)
      
      #Concatenate mean and SD
      i_tukey_groups$Mean_SD <- paste(i_tukey_groups[[1]], i_tukey_groups[[3]], sep=" ± ")
      
      #Append tukey results to the list
      param_group_list[[length(param_group_list)+1]] <- i_tukey_groups[c(2,4)]
      names(param_group_list)[length(param_group_list)] <- names(sig_param_data)[i]
    }
  } else {
    param_group_list <- NULL
  }
  
  
  
  
  #Non-parametric features
  #Perform post-hoc test if there are any significant features
  if (length(sig_no_param_data) > 0) {
    #Format labels to avoid errors
    sig_no_param_data$Group <- gsub(" - ", "_", sig_no_param_data$Group)
    #Create a new list to store the results
    no_param_group_list <- list()
    #Iterate non-parametric dataset
    for (i in 2:ncol(sig_no_param_data)) {
      #Run Dunn post-hoc test
      dunn <- dunnTest(sig_no_param_data[[names(sig_no_param_data)[i]]] ~ sig_no_param_data[["Group"]],
                       data=sig_no_param_data, method="bonferroni")
      #Generate letters
      dunn_letter <- cldList(P.adj ~ Comparison, data = dunn$res, threshold = 0.05)
      
      #Run ANOVA (only to get means and sd)
      i_aov <- aov(lm(sig_no_param_data[[names(sig_no_param_data)[i]]] ~ sig_no_param_data[["Group"]], 
                      sig_no_param_data))
      #Run tukey test (only to get means and sd)
      i_tukey <- HSD.test(i_aov, 'sig_no_param_data[["Group"]]', group=TRUE)
      
      #Get groups
      i_tukey_means <- i_tukey$means
      names(i_tukey_means)[1] <- names(sig_no_param_data)[i]
      names(dunn_letter)[1] <- names(sig_no_param_data)[i]
      
      #Round mean values
      i_tukey_means[1] <- format(round(i_tukey_means[1],2), nsmall=2)
      #Add standard deviations
      i_tukey_means$std <- format(round(i_tukey_means$std,2), nsmall=2)
      
      #Concatenate mean and SD
      dunn_letter$Mean_SD <- paste(i_tukey_means[[1]], i_tukey_means$std, sep=" ± ")
      
      #Format group label
      names(dunn_letter)[3] <- c("groups")
      
      #Append tukey results to the list
      no_param_group_list[[length(no_param_group_list)+1]] <- dunn_letter[c(3,4)]
      names(no_param_group_list)[length(no_param_group_list)] <- names(sig_no_param_data)[i]
    }
  } else {
    no_param_group_list <- NULL
  }
  
  
  
  
  #Export results as list
  #Important: feature number refers to feature number (other than Group) in parametric or 
  #non-parametric datasets, not the original metadata file
  results <- list(param_group_list, no_param_group_list)
  
  #Set list names
  names(results) <- c("Parametric_post_hoc", "Non_parametric_post_hoc")
  
  #Export list
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to analyse multiple variables
multiple_stats <- function(estimator_list, pseq_taxa, metadata_class_no_formula, beta_diversity, gene_count_samples) {
  
  #Create a list to store the results
  results <- list()
  #Iterate the list of estimators
  for (i in 1:length(estimator_list)) {
    #Select numerical data to assess significant differences
    classif_group <- estimator_list[i]
    #Alpha diversity
    if (length(grep("alpha_div", names(estimator_list)[i])) > 0) {
      #Select metadata label to assess significant differences
      classif_meta <- pseq_taxa@sam_data[[metadata_class_no_formula]]
      # classif_meta <- prune_samples(pseq_taxa@sam_data$Time == "1" | pseq_taxa@sam_data$Time == "2", pseq_taxa)@sam_data[[metadata_class_no_formula]]
    }
    #Beta diversity
    if (length(grep("beta_div", names(estimator_list)[i])) > 0) {
      #Select metadata label to assess significant differences
      classif_meta <- beta_diversity$Taxa_diversity[[metadata_class_no_formula]]
    }
    #Gene count
    if (length(grep("gene_count", names(estimator_list)[i])) > 0) {
      #Select metadata label to assess significant differences
      classif_meta <- bind_rows(gene_count_samples[grep("count_data", names(gene_count_samples))], .id = "column_label")[,2]
    }
    
    #Perform statistics
    #Set objects as dataframes
    classif_group <- as.data.frame(classif_group)
    classif_meta <- as.data.frame(classif_meta)
    #Calculate statistically significant differences
    stat_dif <- statistical_tests(classif_meta, classif_group)
    #Get significant features
    #Significant parametric data
    if (length(stat_dif$parametric_significant_features$feature) > 0) {
      sig_param_data <- stat_dif$parametric_data[,c(1,grep(stat_dif$parametric_significant_features$feature,
                                                           names(stat_dif$parametric_data)))]
    } else {
      sig_param_data <- NULL
    }
    #Significant non-parametric data
    if (length(stat_dif$no_parametric_significant_features$feature) > 0) {
      sig_no_param_data <- stat_dif$non_parametric_data[,c(1,grep(stat_dif$no_parametric_significant_features$feature,
                                                                  names(stat_dif$non_parametric_data)))]
    } else {
      sig_no_param_data <- NULL
    }
    
    #Perform post-hoc tests
    letters <- post_hoc(sig_param_data, sig_no_param_data)
    
    #Append the results to a single list
    results[[length(results)+1]] <- stat_dif
    names(results)[[length(results)]] <- paste0(names(estimator_list)[i], "_statistics")
    results[[length(results)+1]] <- letters
    names(results)[[length(results)]] <- paste0(names(estimator_list)[i], "_post_hoc")
  }
  
  #Export results
  return(results)
  
}
  
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to analyse multiple taxonomic clades
multiple_stats_taxa <- function(selected_data) {
  
  #Create a list to store the results
  results <- list()
  #Iterate each taxonomic clade
  for (i in 2:ncol(selected_data)) {
    #Select i numerical data and metadata
    classif_group <- as.data.frame(selected_data[,i,drop=F])
    classif_meta <- as.data.frame(selected_data[,1,drop=F])
    classif_data <- cbind.data.frame(classif_meta, classif_group)
    names(classif_data)[1] <- c("Group")
    #Perform post-hoc tests
    letters <- post_hoc(classif_data, classif_data)
    #Append the results to a single list
    results[[length(results)+1]] <- letters
    names(results)[[length(results)]] <- paste0(names(classif_group), "_post_hoc")
  }
  
  #Export results
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
