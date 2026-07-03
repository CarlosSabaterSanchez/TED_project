
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

#Load libraries
library(dplyr)
library(ggraph)
library(igraph)

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

#Custom function to perform machine learning classification
ml_classification <- function(ml_dataset, compute_models,
                              compute_models_tune, r_objects_path, r_object_name) {
  
  #Create an empty list to summarise modelling results
  model_list <- list()
  #Iterate label levels
  for (i in 1:length(ml_dataset)) {
    #Create a new list of models that will be computed for i dataset
    i_model_list <- list()
    #Get i level
    i_level <- names(ml_dataset)[i]
    #Get i dataset
    i_dataset <- ml_dataset[[i]]
    
    #Remove variables showing near zero variances
    i_dataset_numeric_var <- i_dataset[,2:ncol(i_dataset)]
    # i_dataset_numeric_var <- i_dataset_numeric_var[,-nearZeroVar(i_dataset_numeric_var)]
    # i_dataset <- cbind.data.frame(i_dataset$Label, i_dataset_numeric_var)
    # names(i_dataset)[1] <- c("Label")
    
    #Perform a principal components analysis (PCA) on i dataset
    i_pca_data <- i_dataset
    #Prepare data for PCA model
    i_pca_variables <- i_pca_data[,2:ncol(i_pca_data)]
    names(i_pca_variables) <- gsub("_", " ", names(i_pca_variables))
    i_pca_labels <- gsub("_", " ", i_pca_data$Label)
    #PCA plot showing variables names
    i_pca_plot_with_var_names <- fviz_pca_biplot(prcomp(i_pca_variables), label="var", 
                                                 habillage=i_pca_labels,
                                                 addEllipses=T, ellipse.level=0.95, 
                                                 pointsize=4.4, labelsize=5, arrowsize=1) + 
      theme(text=element_text(size=20))
    #PCA plot without variables names
    i_pca_plot_without_var_names <- fviz_pca_biplot(prcomp(i_pca_variables), label="var", 
                                                    habillage=i_pca_labels,
                                                    addEllipses=T, ellipse.level=0.95, 
                                                    pointsize=4.4, labelsize=0, arrowsize=1) + 
      theme(text=element_text(size=20))
    #Format axis titles
    i_pca_plot_with_var_names$labels$x <- gsub("Dim1", "PC 1", i_pca_plot_with_var_names$labels$x)
    i_pca_plot_with_var_names$labels$y <- gsub("Dim2", "PC 2", i_pca_plot_with_var_names$labels$y)
    i_pca_plot_without_var_names$labels$x <- gsub("Dim1", "PC 1", i_pca_plot_without_var_names$labels$x)
    i_pca_plot_without_var_names$labels$y <- gsub("Dim2", "PC 2", i_pca_plot_without_var_names$labels$y)
    #Append i PCA to the model list
    i_model_list[[length(i_model_list)+1]] <- i_pca_variables
    names(i_model_list)[[length(i_model_list)]] <- paste0("pca_variables")
    i_model_list[[length(i_model_list)+1]] <- i_pca_plot_with_var_names
    names(i_model_list)[[length(i_model_list)]] <- paste0("pca_with_var_names")
    i_model_list[[length(i_model_list)+1]] <- i_pca_plot_without_var_names
    names(i_model_list)[[length(i_model_list)]] <- paste0("pca_without_var_names")
    
    #Prepare the different inputs for supervised machine learning models
    #Generate empty data frames to store train and test variables and classes
    i_train_variables <- data.frame()
    i_train_labels <- data.frame()
    i_test_variables <- data.frame()
    i_test_labels <- data.frame()
    
    #Check class imbalance
    set.seed(600)
    smote_ml_data <- SMOTE(i_dataset[,2:ncol(i_dataset)], target=i_dataset$Label, K=1, dup_size=6)
    #Format SMOTE output
    colnames(smote_ml_data$data)[ncol(smote_ml_data$data)] <- names(i_dataset)[1]
    smote_ml_data$data <- smote_ml_data$data[,names(i_dataset)]
    i_dataset <- rbind.data.frame(i_dataset, smote_ml_data$data)
    
    #Iterate all levels from the factor variable (i.e. microbial taxa)
    for (j in 1:length(levels(i_dataset$Label))) {
      #Subset samples for i level
      j_level <- i_dataset[which(i_dataset$Label == levels(i_dataset$Label)[j]),]
      #Randomly select 70% samples from i_level
      set.seed(600)
      j_random_selection <- sample(1:nrow(j_level),round(0.3*nrow(j_level)))
      #Generate train and test sets from 70 and 30% of samples, respectively
      j_train <- j_level[-j_random_selection,]
      j_test <- j_level[j_random_selection,]
      #Select train and test variables
      j_train_variables <- j_train[,2:ncol(j_train)]
      j_test_variables <- j_test[,2:ncol(j_test)]
      #Select train and test label
      j_train_label <- factor(j_train[,1], levels=unique(j_train[,1]))
      j_test_label <- factor(j_test[,1], levels=unique(j_test[,1]))
      #Append j train and test sets to the previously created dataframes
      #Numeric variables
      i_train_variables <- rbind(i_train_variables, j_train_variables)
      i_test_variables <- rbind(i_test_variables, j_test_variables)
      #Labels
      i_train_labels <- unlist(list(i_train_labels, j_train_label))
      i_test_labels <- unlist(list(i_test_labels, j_test_label))
    }
    
    #Compute different models for i dataset
    for (k in 1:length(compute_models)) {
      #Perform supervised learning
      print(paste0("Generating ", compute_models[k], " model for ", i_level, " dataset"))
      #Train a partial least squares (pls) model
      set.seed(600) 
      k_model <- train(i_train_variables, i_train_labels, method=compute_models[k],
                       preProcess=c("range"),
                       tuneLength=compute_models_tune[i,k], importance=T)
      #Calculate the train accuracy rate
      k_model_predictions_train <- predict(k_model, i_train_variables)
      k_model_confusion_matrix_train <- confusionMatrix(as.factor(i_train_labels), k_model_predictions_train)
      #Calculate the test accuracy rate
      k_model_predictions_test <- predict(k_model, i_test_variables)
      k_model_confusion_matrix_test <- confusionMatrix(as.factor(i_test_labels), k_model_predictions_test)
      #Append k model to the model list
      i_model_list[[length(i_model_list)+1]] <- k_model
      names(i_model_list)[[length(i_model_list)]] <- paste0(compute_models[k], "_model")
      #Append confussion matrix (train phase) of k model to the model list
      i_model_list[[length(i_model_list)+1]] <- k_model_confusion_matrix_train
      names(i_model_list)[[length(i_model_list)]] <- paste0(compute_models[k], "_model_conf_matrix_train")
      #Append confussion matrix (test phase) of k model to the model list
      i_model_list[[length(i_model_list)+1]] <- k_model_confusion_matrix_test
      names(i_model_list)[[length(i_model_list)]] <- paste0(compute_models[k], "_model_conf_matrix_test")
    }
    
    #Compare the performance metrics of different machine learning models studied
    #Selected parameters: Accuracy and Kappa (instead of ROC, sensitivity and specificity which are only for binary classifications)
    #Make statistical statements about models performance differences by first collect the resampling results
    print(paste0("Calculating model performance for ", i_level, " dataset"))
    i_resample_list <- i_model_list[grep("_model$", names(i_model_list))]
    names(i_resample_list) <- toupper(gsub("_model", "", names(i_resample_list)))
    i_resamples <- resamples(i_resample_list)
    #Summarise resampling distributions
    summary(i_resamples)
    #Visualize model comparison
    #Boxplot of model performance metrics
    theme1 <- trellis.par.get()
    theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
    theme1$plot.symbol$pch = 30
    theme1$plot.line$col = rgb(1, 0, 0, .7)
    theme1$plot.line$lwd <- 2
    trellis.par.set(theme1)
    i_model_comparison_boxplot <- bwplot(i_resamples, layout=c(3, 1), cex=3, cex.axis=20)
    #Confidence intervals for accuracy and kappa metrics
    trellis.par.set(caretTheme())
    i_model_comparison_accuracy_plot <- dotplot(i_resamples, metric="Accuracy")
    i_model_comparison_kappa_plot <- dotplot(i_resamples, metric="Kappa")
    #Calculate statistical differences between models, simple t-test
    i_model_comparison_differences <- diff(i_resamples)
    i_model_t_test <- summary(i_model_comparison_differences)
    
    #Append i model boxplot to the model list
    i_model_list[[length(i_model_list)+1]] <- i_model_comparison_boxplot
    names(i_model_list)[[length(i_model_list)]] <- paste0("model_boxplot")
    #Append i model accuracy plot to the model list
    i_model_list[[length(i_model_list)+1]] <- i_model_comparison_accuracy_plot
    names(i_model_list)[[length(i_model_list)]] <- paste0("model_accuracy_plot")
    #Append i model kappa plot to the model list
    i_model_list[[length(i_model_list)+1]] <- i_model_comparison_kappa_plot
    names(i_model_list)[[length(i_model_list)]] <- paste0("model_kappa_plot")
    #Append i model t-test to the model list
    i_model_list[[length(i_model_list)+1]] <- i_model_t_test
    names(i_model_list)[[length(i_model_list)]] <- paste0("model_t_test")
    
    #Create a new dataframe to summarize the performance metrics of all models generated
    model_metrics <- data.frame()
    #Create a new dataframe to summarize variable importance coefficients of all models generated
    model_variable_importance <- data.frame(matrix(nrow=ncol(i_train_variables), ncol=0))
    
    #Get the performance metrics and variable importance coefficients of all models generated
    for (k in 1:length(compute_models)) {
      #k model stats - train phase
      k_train_stats <- i_model_list[[grep(paste0(compute_models[k], "_model_conf_matrix_train"), 
                                          names(i_model_list))]]$byClass
      #Set rownames of k model stats - train phase
      rownames(k_train_stats) <- paste0(toupper(compute_models[k]), " - Train phase: ", gsub("_", " ", gsub("Class: ", "", rownames(k_train_stats))))
      
      #k model stats - test phase
      k_test_stats <- i_model_list[[grep(paste0(compute_models[k], "_model_conf_matrix_test"), 
                                         names(i_model_list))]]$byClass
      #Set rownames of k model stats - test phase
      rownames(k_test_stats) <- paste0(toupper(compute_models[k]), " - Test phase: ", gsub("_", " ", gsub("Class: ", "", rownames(k_test_stats))))
      #Get the performance metrics of k model (including both train and test phases)
      k_stats <- rbind.data.frame(k_train_stats, k_test_stats)
      k_stats <- k_stats[order(gsub(".*: ", "", rownames(k_stats))),]
      #Append k model stats to a single dataframe
      model_metrics <- rbind.data.frame(model_metrics, k_stats)
      
      #Perform a variable importance analysis
      if (compute_models[k] != "rf") {
        if (compute_models[k] != "nnet") {
          print(paste0("Calculating variable importance coefficients for ", compute_models[k], " and ", i_level, " dataset"))
          #Model-specific variable importance metrics are available for several models 
          #General variable importance method implemented in caret package will be used to analyse the rest of models
          #PLS variable importance
          k_importance <- abs(data.frame(varImp(i_model_list[[grep(paste0(compute_models[k], "_model$"), names(i_model_list))]])$importance))
          #Transform variable importance coefficients to relative percentages (%)
          k_importance <- round(k_importance*100/max(k_importance), digits=2)
          #Format row and col names
          rownames(k_importance) <- gsub("_", " ", rownames(k_importance))
          names(k_importance) <- gsub("_", " ", names(k_importance))
        }
      }
      
      #Perform a variable importance analysis specific of random forest models
      if (compute_models[k] == "rf") {
        #Get k random forest final model
        print(paste0("Analyzing ", compute_models[k], " model for ", i_level, " dataset"))
        k_rf <- i_model_list[[grep(paste0(compute_models[k], "_model$"), names(i_model_list))]]$finalModel
        #Random forest variable importance
        k_importance <- as.data.frame(importance(k_rf))
        #Select the global mean decrease accuracy indicator as variable importance estimator
        k_importance <- abs(k_importance[,1:(ncol(k_importance)-2), drop=F])
        #Transform variable importance coefficients to relative percentages (%)
        k_importance <- round(k_importance*100/max(k_importance), digits=2)
        #Format row names
        rownames(k_importance) <- gsub("_", " ", rownames(k_importance))
        
        #Plot random trees
        #Plot the tree with the biggest number of nodes
        k_tree_num <- which(k_rf$forest$ndbigtree == max(k_rf$forest$ndbigtree))
        #Plot the tree with the smallest number of nodes
        # k_tree_num <- which(k_rf$forest$ndbigtree == min(k_rf$forest$ndbigtree))
        #Plot random forest
        k_rf_plot <- tree_func(final_model=k_rf, k_tree_num[1])
        #Append model plot to the model list
        i_model_list[[length(i_model_list)+1]] <- k_rf_plot
        names(i_model_list)[[length(i_model_list)]] <- paste0(compute_models[k], "_plot")
      }
      
      #Perform a variable importance analysis specific of neural network models
      if (compute_models[k] == "nnet") {
        # if (compute_models[k] == "mlp") {
        #Get k neural network final model
        print(paste0("Analyzing ", compute_models[k], " model for ", i_level, " dataset"))
        k_nnet <- i_model_list[[grep(paste0(compute_models[k], "_model$"), names(i_model_list))]]$finalModel
        #Get output labels of k neural network model
        k_nnet_levels <- i_model_list[[grep(paste0(compute_models[k], "_model$"), names(i_model_list))]]$levels
        #Generate an empty list to store variable importance coefficients
        k_importance <- list()
        #Iterate output labels of k neural network model
        for (l in 1:length(k_nnet_levels)) {
          #Calculate neural network variable importance for l output
          l_importance <- olden(k_nnet, bar_plot=F, out_var=paste0("", k_nnet_levels[l]))
          # l_importance <- olden(k_nnet, bar_plot=F, out_var=paste0("Output_", k_nnet_levels[l]))
          #Append neural network variable importance for l output to a single list
          k_importance <- c(k_importance, l_importance)
          #Convert neural network variable importance for l output to a dataframe and set row and col names
          k_importance[[l]] <- as.data.frame(k_importance[[l]])
          names(k_importance[[l]]) <- gsub("_", " ", k_nnet_levels[l])
          rownames(k_importance[[l]]) <- gsub("_", " ", gsub("Input_", "", rownames(l_importance)))
        }
        #Convert the list of variable importance coefficients to a dataframe
        k_importance <- as.data.frame(k_importance)
        names(k_importance) <- gsub("[.]", " ", names(k_importance))
        #Transform variable importance coefficients to relative percentages (%)
        k_importance <- round(k_importance*100/max(k_importance), digits=2)
        
        #Plot neural network
        #This plot takes time, it should be generated for networks computed with a reduced number of input variables
        if (ncol(i_pca_variables) <= 50) {
          #Set graphical parameters
          x11()
          par(oma=c(0.001, 0.001, 0.001, 0.001))
          par(mar=c(0.001, 0.001, 0.001, 0.001))
          par(mar = numeric(4), family = 'serif')
          plotnet(k_nnet, alpha=0.3, cex_val=1.0, max_sp=T, pad_x = 0.7,
                  circle_cex=2, node_labs=F, bias=F)
          k_nnet_plot <- recordPlot()
          dev.off()
          #Append model plot to the model list
          i_model_list[[length(i_model_list)+1]] <- k_nnet_plot
          names(i_model_list)[[length(i_model_list)]] <- paste0(compute_models[k], "_plot")
        }
      }
      
      #Append k model variable importance results to a single dataframe
      names(k_importance) <- paste0(toupper(compute_models[k]), " (%): ", names(k_importance))
      model_variable_importance <- cbind.data.frame(model_variable_importance, k_importance)
      
    }
    
    #Set column names
    names(model_metrics) <- paste0(names(model_metrics), " (%)")
    names(model_metrics) <- gsub("[.]", " ", names(model_metrics))
    #Round performance metrics values
    model_metrics <- format(round(model_metrics*100, digits=2), nsmall=2)
    #Append model performance metrics to the model list
    i_model_list[[length(i_model_list)+1]] <- model_metrics
    names(i_model_list)[[length(i_model_list)]] <- paste0("model_metrics")
    
    #Append model variable importance to the model list
    model_variable_importance <- model_variable_importance[,order(gsub(".*:", "", names(model_variable_importance)))]
    i_model_list[[length(i_model_list)+1]] <- model_variable_importance
    names(i_model_list)[[length(i_model_list)]] <- paste0("model_variable_importance")
    
    #Generate a correlation network of predictor variables
    print(paste0("Generating a correlation network for ", i_level, " dataset"))
    #Calculate Pearson correlation coefficients for enzyme domains
    correlation_data <- i_dataset[,2:ncol(i_dataset)]
    names(correlation_data) <- gsub("_", " ", names(correlation_data))
    #Calculate correlations
    i_correlations <- correlate(correlation_data)
    #Add a norm noise to prevent errors during plotting
    i_correlations[, -1] <- i_correlations[, -1] + rnorm(cumprod(dim(i_correlations[, -1]))[-1], sd = 0.001)
    #Generate correlation plots
    i_correlation_plot <- network_plot(i_correlations, min_cor=0.5)
    i_correlation_plot <- i_correlation_plot + theme(text=element_text(size=25), legend.text=element_text(face="bold")) 
    #Append correlation network to the model list
    i_model_list[[length(i_model_list)+1]] <- i_correlations
    names(i_model_list)[[length(i_model_list)]] <- paste0("correlation_matrix")
    i_model_list[[length(i_model_list)+1]] <- i_correlation_plot
    names(i_model_list)[[length(i_model_list)]] <- paste0("correlation_network")
    
    #Append i model list to a general list summarising all datasets
    model_list[[length(model_list)+1]] <- i_model_list
    names(model_list)[[length(model_list)]] <- i_level
  }
  
  #Export R objects
  saveRDS(model_list, paste0(r_objects_path, r_object_name, ".rds"))
  #Return results
  return(model_list)
  
}

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

#Custom function to plot random forest trees
tree_func <- function(final_model, 
                      tree_num) {
  #Get tree by index
  tree <- randomForest::getTree(final_model, 
                                k=tree_num, 
                                labelVar=TRUE) %>%
    tibble::rownames_to_column() %>%
    #Make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point`=ifelse(is.na(prediction), `split point`, NA))
  
  #Prepare data frame for graph
  graph_frame <- data.frame(from=rep(tree$rowname, 2),
                            to=c(tree$`left daughter`, tree$`right daughter`))
  
  #Convert to graph and delete the last node that we don't want to plot
  graph <- graph_from_data_frame(graph_frame) %>%
    delete_vertices("0")
  
  #Set node labels
  V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  V(graph)$leaf_label <- as.character(tree$prediction)
  V(graph)$leaf_label <- gsub("_", " ", V(graph)$leaf_label)
  V(graph)$split <- as.character(round(tree$`split point`, digits=2))
  
  #Plot
  plot <- ggraph(graph, 'dendrogram') + 
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label=node_label), na.rm=TRUE, repel=TRUE) +
    geom_node_label(aes(label=split), vjust=2.5, na.rm=TRUE, fill="white") +
    geom_node_label(aes(label=leaf_label, fill=leaf_label), na.rm=TRUE, 
                    repel=TRUE, colour="white", fontface="bold", show.legend=FALSE) +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_blank(),
          plot.background=element_rect(fill="white"),
          panel.border=element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_text(size=30))
  #Print plot
  print(plot)
  
}

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
