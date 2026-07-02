
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to generate heatmap and dendogram
heatmap_dendogram <- function(data_to_plot, x_size, y_size, 
                              palette_to_use, 
                              strain_levels) {
  
  #Load libraries
  library(tidyverse)
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  library(reshape2)
  library(grid)
  library(stringi)
  library(gtools)
  library(Cairo)
  library(ComplexHeatmap)
  
  #Generate heatmap and dendogram plots
  #Replace NAs by 0
  data_to_plot[is.na(data_to_plot)] <- 0
  #Remove categorical variable
  data_to_plot_no_class <- data_to_plot[,2:ncol(data_to_plot),drop=F]
  
  #Format microbial classification colors
  #Format dendogram leaves
  #https://jokergoo.github.io/2023/07/07/what-is-a-hclust-object/
  dendrogram <- as.dendrogram(hclust(d=dist(x=data_to_plot_no_class)))
  dendrogram <- rev(dendrogram)
  dendrogram_structure <- capture.output(str(dendrogram))
  dendrogram_structure <- rev(dendrogram_structure)
  dendrogram_leaves <- grep("leaf", dendrogram_structure)
  #Remove special characters to align dendogram and microbial metadata
  data_to_plot_no_special_chars <- data_to_plot
  rownames(data_to_plot_no_special_chars) <- gsub("[[:punct:]]", "", rownames(data_to_plot_no_special_chars))
  data_to_plot_no_special_chars <- data_to_plot_no_special_chars[gsub("[[:punct:]]", "", dendro_data(dendrogram)$labels$label),]
  all(rownames(data_to_plot_no_special_chars) == gsub("[[:punct:]]", "", dendro_data(dendrogram)$labels$label))
  #Assign colors to each level
  microbial_classification_colors <- factor(data_to_plot_no_special_chars$Class, levels=strain_levels)
  microbial_classification_levels <- levels(microbial_classification_colors)
  for (i in 1:length(microbial_classification_levels)) {
    #Assign i color to i level
    microbial_classification_colors <- gsub(paste0("^", microbial_classification_levels[i], "$"), palette_to_use_dendogram[i], 
                                            microbial_classification_colors)
  }
  
  #Plot heatmap
  #https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap
  #https://yunuuuu.github.io/ggalign/articles/plot-options.html
  heatmap_plot <- ggheatmap(data_to_plot_no_class) + 
    geom_tile(aes(fill=value), width=1, height=1, color="lightgrey") +
    scale_fill_gradient2(low="white", high="black") +
    scale_y_continuous(position="right") +
    theme(axis.text.y.right=element_markdown(), 
          axis.text.x=element_text(size=x_size, angle=90, vjust=0.5, hjust=1, face="bold"),
          axis.text.y=element_text(size=y_size),
          legend.position="none",
          panel.border=element_rect(colour="black", fill=NA, size=1)) + 
    hmanno("r", size=unit(3, "cm")) +
    # align_dendro(distance="euclidean", method="complete", k=3L) + #Discriminates three leaves
    align_dendro(data=data_to_plot[,2:ncol(data_to_plot)], distance="euclidean", method="complete") + #Does not discriminates three leaves
    geom_point(size=4, data=~.x[dendrogram_leaves,], color=rev(microbial_classification_colors))
  
  #Append results
  plots <- list(heatmap_plot)
  names(plots) <- c("heatmap_plot")
  #Return results
  return(plots)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to generate heatmap and dendogram with legend
heatmap_dendogram_legend <- function(data_to_plot, x_size, y_size, 
                              palette_to_use, 
                              strain_levels, scale_text) {
  
  #Load libraries
  library(tidyverse)
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  library(reshape2)
  library(grid)
  library(stringi)
  library(gtools)
  library(Cairo)
  library(ComplexHeatmap)
  
  #Generate heatmap and dendogram plots
  #Replace NAs by 0
  data_to_plot[is.na(data_to_plot)] <- 0
  #Remove categorical variable
  data_to_plot_no_class <- data_to_plot[,2:ncol(data_to_plot),drop=F]
  
  #Format microbial classification colors
  #Format dendogram leaves
  #https://jokergoo.github.io/2023/07/07/what-is-a-hclust-object/
  dendrogram <- as.dendrogram(hclust(d=dist(x=data_to_plot_no_class)))
  dendrogram <- rev(dendrogram)
  dendrogram_structure <- capture.output(str(dendrogram))
  dendrogram_structure <- rev(dendrogram_structure)
  dendrogram_leaves <- grep("leaf", dendrogram_structure)
  #Remove special characters to align dendogram and microbial metadata
  data_to_plot_no_special_chars <- data_to_plot
  rownames(data_to_plot_no_special_chars) <- gsub("[[:punct:]]", "", rownames(data_to_plot_no_special_chars))
  data_to_plot_no_special_chars <- data_to_plot_no_special_chars[gsub("[[:punct:]]", "", dendro_data(dendrogram)$labels$label),]
  all(rownames(data_to_plot_no_special_chars) == gsub("[[:punct:]]", "", dendro_data(dendrogram)$labels$label))
  #Assign colors to each level
  microbial_classification_colors <- factor(data_to_plot_no_special_chars$Class, levels=strain_levels)
  microbial_classification_levels <- levels(microbial_classification_colors)
  for (i in 1:length(microbial_classification_levels)) {
    #Assign i color to i level
    microbial_classification_colors <- gsub(paste0("^", microbial_classification_levels[i], "$"), palette_to_use_dendogram[i], 
                                            microbial_classification_colors)
  }
  
  #Plot heatmap
  #https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap
  #https://yunuuuu.github.io/ggalign/articles/plot-options.html
  heatmap_plot <- ggheatmap(data_to_plot_no_class) + 
    geom_tile(aes(fill=value), width=1, height=1, color="lightgrey") +
    scale_fill_gradient2(low="white", high="black", name=scale_text) +
    scale_y_continuous(position="right") +
    theme(axis.text.y.right=element_markdown(), 
          axis.text.x=element_text(size=x_size, angle=90, vjust=0.5, hjust=1, face="bold"),
          axis.text.y=element_text(size=y_size),
          legend.position="top", legend.text=element_text(size=20, face = "bold"), legend.title=element_text(size=20, face="bold", vjust=1),
          panel.border=element_rect(colour="black", fill=NA, size=1)) + 
    hmanno("r", size=unit(3, "cm")) +
    # align_dendro(distance="euclidean", method="complete", k=3L) + #Discriminates three leaves
    align_dendro(data=data_to_plot[,2:ncol(data_to_plot)], distance="euclidean", method="complete") + #Does not discriminates three leaves
    geom_point(size=4, data=~.x[dendrogram_leaves,], color=rev(microbial_classification_colors))
  
  #Append results
  plots <- list(heatmap_plot)
  names(plots) <- c("heatmap_plot")
  #Return results
  return(plots)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to generate heatmap and dendogram using two color scales
heatmap_dendogram_two_color_scales <- function(data_to_plot, x_size, y_size,
                                               palette_to_use_dendogram, palette_to_use_heatmap,
                                               strain_levels,
                                               enzyme_classification, x_axis_legend, x_axis_legend_rows) {
  
  #Load libraries
  library(tidyverse)
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  library(reshape2)
  library(grid)
  library(stringi)
  library(gtools)
  library(Cairo)
  library(ComplexHeatmap)
  
  #Generate heatmap and dendogram plots
  #Replace NAs by 0
  data_to_plot[is.na(data_to_plot)] <- 0
  #Remove categorical variable
  data_to_plot_no_class <- data_to_plot[,2:ncol(data_to_plot),drop=F]
  
  #Order enzyme classification data
  enzyme_classification <- enzyme_classification[colnames(data_to_plot)[2:length(colnames(data_to_plot))],]
  
  #Add a new color scale to x axis
  #Careful: do not subset a dataframe using a factor, convert to character first
  #Format enzyme classification colors
  enzyme_classification_colors <- factor(enzyme_classification$Group, levels=levels(enzyme_classification$Group))
  enzyme_classification_levels <- levels(enzyme_classification_colors)
  for (i in 1:length(enzyme_classification_levels)) {
    #Assign i color to i level
    enzyme_classification_colors <- gsub(paste0("^", enzyme_classification_levels[i], "$"), palette_to_use_heatmap[i], 
                                         enzyme_classification_colors)
  }
  #Format microbial classification colors
  #Format dendogram leaves
  #https://jokergoo.github.io/2023/07/07/what-is-a-hclust-object/
  dendrogram <- as.dendrogram(hclust(d=dist(x=data_to_plot_no_class)))
  dendrogram <- rev(dendrogram)
  dendrogram_structure <- capture.output(str(dendrogram))
  dendrogram_structure <- rev(dendrogram_structure)
  dendrogram_leaves <- grep("leaf", dendrogram_structure)
  #Remove special characters to align dendogram and microbial metadata
  data_to_plot_no_special_chars <- data_to_plot
  rownames(data_to_plot_no_special_chars) <- gsub("[[:punct:]]", "", rownames(data_to_plot_no_special_chars))
  data_to_plot_no_special_chars <- data_to_plot_no_special_chars[gsub("[[:punct:]]", "", dendro_data(dendrogram)$labels$label),]
  all(rownames(data_to_plot_no_special_chars) == gsub("[[:punct:]]", "", dendro_data(dendrogram)$labels$label))
  #Assign colors to each level
  microbial_classification_colors <- factor(data_to_plot_no_special_chars$Class, levels=strain_levels)
  microbial_classification_levels <- levels(microbial_classification_colors)
  for (i in 1:length(microbial_classification_levels)) {
    #Assign i color to i level
    microbial_classification_colors <- gsub(paste0("^", microbial_classification_levels[i], "$"), palette_to_use_dendogram[i], 
                                            microbial_classification_colors)
  }
  
  #Plot heatmap
  #https://stackoverflow.com/questions/42047896/joining-a-dendrogram-and-a-heatmap
  #https://yunuuuu.github.io/ggalign/articles/plot-options.html
  heatmap_plot <- ggheatmap(data_to_plot_no_class) + theme(axis.text.x=element_text(colour=enzyme_classification_colors))  +
    geom_tile(aes(fill=value), width=1, height=1, color="lightgrey") +
    scale_fill_gradient2(low="white", high="black") +
    scale_y_continuous(position="right") +
    theme(axis.text.y.right=element_markdown(), 
          axis.text.x=element_text(size=x_size, angle=90, vjust=0.5, hjust=1, face="bold"),
          axis.text.y=element_text(size=y_size),
          legend.position="none",
          panel.border=element_rect(colour="black", fill=NA, size=1)) + 
    hmanno("r", size=unit(3, "cm")) +
    # align_dendro(distance="euclidean", method="complete", k=3L) + #Discriminates three leaves
    align_dendro(data=data_to_plot[,2:ncol(data_to_plot)], distance="euclidean", method="complete") + #Does not discriminates three leaves
    geom_point(size=4, data=~.x[dendrogram_leaves,], color=rev(microbial_classification_colors))
  
  #Add a custom legend
  lgd <- Legend(labels=enzyme_classification_levels, title="", labels_gp=gpar(fontsize=x_axis_legend), 
               nrow=x_axis_legend_rows, legend_gp=gpar(fill=palette_to_use_heatmap))
  # draw(lgd, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
  # heatmap_plot_legend <- recordPlot()
  heatmap_plot_legend <- lgd
  
  #Append results
  plots <- list(heatmap_plot, heatmap_plot_legend)
  names(plots) <- c("heatmap_plot", "heatmap_plot_legend")
  #Return results
  return(plots)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to generate heatmap and dendogram using two color scales and a intensity scale
heatmap_dendogram_two_color_scales_no_dendo_with_scale <- function(data_to_plot, x_size, y_size, 
                                                                   palette_to_use_dendogram, palette_to_use_heatmap,
                                                                   strain_levels,
                                                                   enzyme_classification, x_axis_legend, x_axis_legend_rows,
                                                                   y_axis_legend, y_axis_legend_rows, scale_text, scale_size) {
  
  #Load libraries
  library(tidyverse)
  library(ggplot2)
  library(ggdendro)
  library(dendextend)
  library(reshape2)
  library(grid)
  library(stringi)
  library(gtools)
  library(Cairo)
  library(ComplexHeatmap)
  
  #Generate heatmap and dendogram plots
  #Replace NAs by 0
  data_to_plot[is.na(data_to_plot)] <- 0
  #Remove categorical variable
  data_to_plot_no_class <- data_to_plot[,2:ncol(data_to_plot)]
  #Transpose data to avoid errors during heatmap plotting
  # data_to_plot_t <- t(data_to_plot) #Transpose function does not work with dataframes
  data_to_plot_t_no_class <- t(data_to_plot[,2:ncol(data_to_plot),drop=F])
  
  #Format data for heatmap plotting
  heatmap_data <- melt(t(data_to_plot_t_no_class), id=colnames(t(data_to_plot_t_no_class)))
  #Set additional heatmap labels
  #Careful: do not subset a dataframe using a factor, convert to character first
  heatmap_data$group <- enzyme_classification[as.character(heatmap_data$Var2),"Group"]
  heatmap_data$group <- factor(heatmap_data$group, levels=levels(enzyme_classification$Group))
  #Create the heatmap plot
  heatmap_plot <- ggplot(data=heatmap_data, aes(x=Var2, y=Var1, fill=value, label=group)) + theme_classic() + 
    geom_tile(aes(fill=value), width=1, height=1, color="lightgrey") +
    scale_fill_gradient2(name=scale_text, 
                         low="white", high="black") + labs(x="") +
    scale_y_discrete(position="right", limits=rev, expand=c(0,0)) + 
    scale_x_discrete(expand=c(0,0)) +
    # geom_tile(colour="gray50") +
    theme(text=element_text(size=x_size), legend.text=element_text(size=scale_size), legend.title=element_text(size=scale_size),
          # axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          legend.position="top", 
          plot.margin=unit(c(1,1,1,1), 'cm'),
          axis.text.x=element_text(angle=90, vjust=0.5, hjust=1, face="bold"),
          axis.text.y=element_text(vjust=0.5, hjust=1, face="bold"))
    
  #Add a new color scale to x axis
  #Careful: do not subset a dataframe using a factor, convert to character first
  enzyme_classification_colors <- enzyme_classification[as.character(unique(heatmap_data$Var2)), "Group"]
  enzyme_classification_levels <- levels(enzyme_classification_colors)
  for (i in 1:length(enzyme_classification_levels)) {
    #Assign i color to i level
    enzyme_classification_colors <- gsub(paste0("^", enzyme_classification_levels[i], "$"), palette_to_use_heatmap[i], 
                                         enzyme_classification_colors)
  }
  
  #Add a custom legend
  heatmap_plot <- heatmap_plot + theme(axis.text.x=element_text(colour=enzyme_classification_colors)) 
  grid.newpage()
  # plot(heatmap_plot)
  lgd <- Legend(labels=enzyme_classification_levels, title="", labels_gp=gpar(fontsize=x_axis_legend), 
                nrow=x_axis_legend_rows, legend_gp=gpar(fill=palette_to_use_heatmap))
  # draw(lgd, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
  # heatmap_plot_legend <- recordPlot()
  heatmap_plot_legend_x <- lgd
  
  #Add a new color scale to y axis
  #Careful: do not subset a dataframe using a factor, convert to character first
  y_classification_colors <- data_to_plot$Class
  y_classification_levels <- levels(y_classification_colors)
  for (i in 1:length(y_classification_levels)) {
    #Assign i color to i level
    y_classification_colors <- gsub(paste0("^", y_classification_levels[i], "$"), palette_to_use_dendogram[i], 
                                    y_classification_colors)
  }
  
  #Add a custom legend
  heatmap_plot <- heatmap_plot + theme(axis.text.y=element_text(colour=rev(y_classification_colors))) 
  grid.newpage()
  # plot(heatmap_plot)
  lgd <- Legend(labels=y_classification_levels, title="", labels_gp=gpar(fontsize=y_axis_legend), 
                nrow=y_axis_legend_rows, legend_gp=gpar(fill=palette_to_use_dendogram))
  # draw(lgd, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
  # heatmap_plot_legend <- recordPlot()
  heatmap_plot_legend_y <- lgd
  
  #Append results
  plots <- list(heatmap_plot, heatmap_plot_legend_x, heatmap_plot_legend_y)
  names(plots) <- c("heatmap_plot", "heatmap_plot_legend_x", "heatmap_plot_legend_y")
  #Return results
  return(plots)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Finish
