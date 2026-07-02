
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to calculate and plot significant correlation coefficients from two datasets
correlation_pvalue_plot <- function(compare_group_1, compare_group_2, coef_threshold, axis_size) {
  
  #Prune rows containing NAs
  #Replace NAs by 0
  compare_group_1[is.na(compare_group_1)] <- 0
  compare_group_2 <- compare_group_2[rownames(compare_group_1),]
  #Generate an empty dataframe to store correlation coefficients and pvalues
  coefficients_and_pvalues <- c()
  #Iterate dataframes
  for (i in 1:ncol(compare_group_1)) {
    for (j in 1:ncol(compare_group_2)) {
      #Calculate correlation coefficient
      i_j_correlation <- cor.test(as.matrix(compare_group_1[,i]), 
                                  as.matrix(compare_group_2[,j]), 
                                  method = "pearson")
      #Get p_value
      i_j_pvalue <- i_j_correlation$p.value
      #Get coefficient
      i_j_coefficient <- i_j_correlation$estimate
      #Get correlation info
      i_j_info <- cbind.data.frame(names(compare_group_1)[i], names(compare_group_2)[j], i_j_pvalue, i_j_coefficient)
      names(i_j_info) <- c("Dataset_1", "Dataset_2", "P_value", "Coefficient")
      #Append correlation info to the complete correlation dataframe
      coefficients_and_pvalues <- rbind(coefficients_and_pvalues, i_j_info)
    }
  }
  
  #Subset significant correlations
  significant_correlations <- subset(coefficients_and_pvalues, coefficients_and_pvalues$P_value < 0.05)
  #Prune coefficients according to a threshold
  significant_correlations <- subset(significant_correlations, abs(significant_correlations$Coefficient) >= coef_threshold)
  
  #Calculate correlation matrix considering those variables showing significant correlations
  correlation_matrix <- cor(compare_group_1[,unique(significant_correlations$Dataset_1)], 
                            compare_group_2[,unique(significant_correlations$Dataset_2)], 
                            method = c("pearson"))
  
  #Delete empty rows (NA) from correlation matrix
  correlation_matrix <- na.omit(correlation_matrix)
  
  #Order correlation matrix alphabetically
  correlation_matrix <- correlation_matrix[order(rownames(correlation_matrix)),]
  correlation_matrix <- correlation_matrix[,order(colnames(correlation_matrix))]
  
  #Plot correlation heatmap
  corrplot::corrplot(correlation_matrix, 
                     tl.cex = 0.7, cl.cex = 1.3, number.cex = 1.0, pch.cex =1.0)
  #significant_correlations_plot <- recordPlot()
  
  #Melt this matrix and remove NA values:
  correlation_matrix_melt <- melt(correlation_matrix, na.rm = TRUE)
  
  #Prune small correlation coefficients lower than -0.3 or +0.3 (label) 
  cell_label <- correlation_matrix_melt
  cell_label$value[which(abs(cell_label$value) < coef_threshold)] <- NA
  
  #Plot correlation matrix using ggplot2
  #Save ggplot2 graph
  significant_correlations_plot <- ggplot(data=correlation_matrix_melt, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() +
    theme(plot.title=element_text(hjust=0.5, face="bold", colour="black"), 
          text=element_text(size=axis_size), 
          axis.text.y.left=element_markdown(),
          axis.text.x=element_markdown(angle = 90, vjust = 0.5, hjust=1, face="bold"),
          legend.title=element_text(face="bold", colour="black"),
          #panel.grid.major=element_line(color="#8ccde3", size=0.75, linetype=2),
          panel.background=element_blank(),
          panel.border=element_blank()) +
    scale_fill_gradient2(limits=c(-1, 1), breaks=seq(-1, 1, by=0.5), 
                         low = "red", mid="white", high = "blue", name="Correlations") +
    guides(size=guide_legend(), colour=guide_legend(reverse=T)) +
    geom_tile(color="black",
              lwd=1.5,
              linetype=1) +
    labs(title="Correlation analysis", x="", y ="") + scale_y_discrete(limits=rev)
  
  #Return correlation results as a list
  correlation_results <- list(significant_correlations, correlation_matrix, significant_correlations_plot, correlation_matrix_melt)
  names(correlation_results) <- c("Significant_correlations", "Correlation_matrix", "Significant_correlations_plot", "correlation_matrix_melt")
  return(correlation_results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to perform a multiomic analysis
multiomics <- function(multiomic_input, selected_levels, label_names, quantile_keep_x, circos_cutoff, zero_var, color_scale, top_var) {
  
  #Create a list to store the results
  results_list <- list()
  
  #Iterate input datasets
  for (i in 1:length(multiomic_input)) {
    #Subset i study
    i_study <- multiomic_input[[i]]
    #Create new lists to store X and Y inputs for multiomic methods
    i_X <- list()
    i_Y <- list()
    #Iterate each dataset within each study
    for (j in 1:length(i_study)) {
      #Subset j dataset
      j_dataset <- i_study[[j]]
      #Get j Y labels
      i_Y[[length(i_Y)+1]] <- j_dataset[,1]
      names(i_Y)[[length(i_Y)]] <- names(i_study)[j]
      #Get j X numeric data
      i_X[[length(i_X)+1]] <- j_dataset[,2:ncol(j_dataset)]
      names(i_X)[[length(i_X)]] <- names(i_study)[j]
    }
    #Check all labels in each dataset are equal
    allSame <- function(x) length(unique(x)) == 1
    print(paste0("Labels from each dataset match: ", allSame(i_Y)))
    #Convert Y to factor
    i_Y <- factor(i_Y[[1]], levels=selected_levels[[i]])
    
    #1- Multiblock PLS-DA
    results_diablo <- block.plsda(i_X, i_Y, scale=T, near.zero.var=zero_var) #Run the method
    
    #2- Multiblock sPLS-DA
    #Set the number of features to use for the X datasets
    #Prepare Keep X list to compute the model
    keep_X_text <- paste0("list(",
                          paste0(paste0(names(i_X), "=round(quantile(1:ncol(i_X$",names(i_X), "),",
                                        paste0("c(", paste0(quantile_keep_x, collapse=","), ")"),"))"), collapse=","),
                          ")")
    keep_X <- eval(parse(text=keep_X_text))
    #Run the method
    results_sparse_diablo <- block.splsda(i_X, i_Y, keepX=keep_X, scale=T, near.zero.var=zero_var) 
    
    #Create a list of different models
    model_list <- list(Diablo=results_diablo, Sparse_diablo=results_sparse_diablo)
    i_results <- list()
    
    #Iterate the list of different models to generate plots
    for (j in 1:length(model_list)) {
      
      #Print a message
      print(paste0("Analyzing model: ",  names(model_list)[j]))
      
      #Plot samples from each dataset
      sample_plot <- plotIndiv(model_list[[j]],
                               ind.names=F, legend=T,
                               ellipse=F, star=T, group=i_Y,
                               pch=i_Y,
                               X.label='Component 1', Y.label='Component 2',
                               subtitle=names(i_X),
                               legend.title=label_names[i],
                               cex=10, size.subtitle=30, size.xlabel=40, size.ylabel=40,
                               size.axis=40, size.legend=20, size.legend.title=20, point.lwd=2, col.per.group=color_scale)

      #Plot arrows
      arrow_plot <- plotArrow(model_list[[j]],
                              ind.names=F, legend=T,
                              ellipse=F, star=T, group=i_Y,
                              pch=i_Y,
                              X.label='Component 1', Y.label='Component 2',
                              subtitle=names(i_X),
                              legend.title=label_names[i],
                              cex=5, size.subtitle=20, size.xlabel=25, size.ylabel=25,
                              size.axis=25, size.legend=25, size.legend.title=25, col.per.group=color_scale) #Plot the samples
      
      #Plot all variables
      x11(width=14,height=7.5) #Must be full screen
      plotVar(model_list[[j]], legend=TRUE, pch=rep(16, length(i_X))) #Plot the variables
      all_var_plot <- recordPlot()
      dev.off()
      
      #Diablo plot
      x11(width=14,height=7.5) #Must be full screen
      plotDiablo(model_list[[j]], ncomp=1, legend=T, col.per.group=color_scale)
      diablo_plot <- recordPlot()
      dev.off()
      
      #Variable selection (importance)
      var_importance_comp_1 <- selectVar(model_list[[j]], comp=1)
      var_importance_comp_2 <- selectVar(model_list[[j]], comp=2)
      
      #Circos plot - Correlations between datasets
      x11(width=14,height=7.5) #Must be full screen
      circosPlot(model_list[[j]], cutoff=circos_cutoff, line=F, legend=F,
                 color.blocks=color_scale,
                 color.cor=c("darkblue","red"), size.labels=0.01, size.variables=0.6,
                 size.legend=1.0)
      circos_plot <- recordPlot()
      dev.off()
      
      #Relevance network
      # x11(width=14,height=7.5)
      # network(model_list[[j]], blocks=c(1:length(i_X)),
      #         # color.node=c('darkorchid', 'brown1', 'lightgreen', 'lightblue'),
      #         cutoff=0.5)
      # relevance_network <- recordPlot()
      # dev.off()
      
      #Clustered image map specifically implemented to represent the multi-omics molecular signature expression for each sample
      if (length(grep("block.splsda", class(model_list[[j]]))) != 0) {
        x11(width=14,height=7.5)
        cimDiablo(model_list[[j]], margins=c(10,20), legend.position="right", size.legend=1.0)
        cluster_plot <- recordPlot()
        dev.off()
      } else {
        cluster_plot <- NULL
      }
      
      #Generate individual figures per dataset
      k_results <- list()
      for (k in 1:length(i_X)) {
        #Prepare a vector to select k figure
        selection_vector <- rep(F, length(i_X))
        selection_vector[k] <- T
        #Plot k variables
        x11(width=14,height=7.5) #Must be full screen
        plotVar(model_list[[j]], var.names=selection_vector,
                legend=T, pch=rep(16, length(i_X)))
        k_var_plot <- recordPlot()
        dev.off()
        
        #Plot the contributions of each feature to each dimension
        #Most influent variables
        x11(width=14,height=10) #Must be full screen
        par(oma=c(0,12,0,0))
        plotLoadings(model_list[[j]], ncomp=1, contrib="max",
                     legend.title=label_names[i], block=k, #name.var=colnames(i_X[[k]]),
                     title=names(i_X)[k], size.title=2, size.legend=1.3, size.name=1.5, legend.color=color_scale,
                     ndisplay=top_var) #Variable size: size.name=0.7
        k_importance_plot <- recordPlot()
        dev.off()
        
        #ROC curves
        roc_comp_1 <- auroc(model_list[[j]], roc.block=names(i_X)[k], roc.comp=1)
        roc_comp_2 <- auroc(model_list[[j]], roc.block=names(i_X)[k], roc.comp=2)
        
        #Append results to a single list
        k_results[[length(k_results)+1]] <- k_var_plot
        names(k_results)[[length(k_results)]] <- paste0(names(i_X)[k], "_var_plot")
        k_results[[length(k_results)+1]] <- k_importance_plot
        names(k_results)[[length(k_results)]] <- paste0(names(i_X)[k], "_importance_plot")
        k_results[[length(k_results)+1]] <- roc_comp_1
        names(k_results)[[length(k_results)]] <- paste0(names(i_X)[k], "_ROC_comp_1")
        k_results[[length(k_results)+1]] <- roc_comp_2
        names(k_results)[[length(k_results)]] <- paste0(names(i_X)[k], "_ROC_comp_2")
      }
      
      #Append results from j model to a single list
      j_results <- list(Model=model_list[[j]],
                        Sample_plot=sample_plot, Arrow_plot=arrow_plot, Diablo_plot=diablo_plot,
                        All_var_plot=all_var_plot, 
                        Individual_plots=k_results,
                        Var_importance_comp_1=var_importance_comp_1, Var_importance_comp_2=var_importance_comp_2,
                        Circos_plot=circos_plot, 
                        # Relevance_network=relevance_network, 
                        Cluster_plot=cluster_plot)
      i_results[[length(i_results)+1]] <- j_results
      names(i_results)[[length(i_results)]] <- names(model_list)[j]
    }
    
    #Append all results to a single list
    results_list[[length(results_list)+1]] <- i_results
    names(results_list)[[length(results_list)]] <- names(multiomic_input)[i]
    
  }
  
  #Export results list
  return(results_list)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
