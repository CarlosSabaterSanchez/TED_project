#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to calculate alpha-diversity in pseq objects contaning functional metagenomics data
pseq_alpha_diversity <- function(pseq_div, metadata_class_no_formula, selected_gene_count, results_name, selected_levels) {
  
  #Alpha diversity of taxonomic data
  #Multiply percentages * reads to get counts and round the result
  # for (i in 1:ncol(pseq_div@otu_table)) {
  #   pseq_div@otu_table[,i] <- pseq_div@otu_table[,i]*metadata$`Reads Knead`[i]
  # }
  #Scale OTU data: multiply relative proportion per 1e+06
  # pseq_div@otu_table <- pseq_div@otu_table*1e+06/max(pseq_div@otu_table)
  # pseq_div@otu_table <- round(pseq_div@otu_table)
  
  #Transform class to factor
  pseq_div@sam_data[[metadata_class_no_formula]] <- factor(pseq_div@sam_data[[metadata_class_no_formula]], levels=selected_levels)
  
  #Alpha diversity (especially for all taxonomic classes combined and numerous taxa)
  alpha_div <- alpha(pseq_div, index=c("Chao1", "Shannon", "diversity_gini_simpson", "diversity_inverse_simpson"))
  #Chao1
  global_div <- paste(round(mean(alpha_div$chao1),1), "±", round(sd(alpha_div$chao1),1))
  #Get chao1 values for each level
  diversity_per_level <- list()
  for (i in 1:length(levels(pseq_div@sam_data[[metadata_class_no_formula]]))) {
    i_level <- levels(pseq_div@sam_data[[metadata_class_no_formula]])[i]
    i_mean <- mean(alpha(prune_samples(pseq_div@sam_data[[metadata_class_no_formula]]==i_level, pseq_div), index=c("Chao1", "Shannon", "diversity_gini_simpson", "diversity_inverse_simpson"))$chao1)
    i_sd <- sd(alpha(prune_samples(pseq_div@sam_data[[metadata_class_no_formula]]==i_level, pseq_div), index=c("Chao1", "Shannon", "diversity_gini_simpson", "diversity_inverse_simpson"))$chao1)
    diversity_per_level <- c(diversity_per_level, paste(round(i_mean,1), "±", round(i_sd,1)))
    names(diversity_per_level)[length(diversity_per_level)] <- i_level
  }
  #Inspect Chao1 values
  div_summary <- as.data.frame(unlist(diversity_per_level))
  names(div_summary) <- c("Chao1")
  xlsx::write.xlsx(div_summary, paste0("Results/", results_name, "_chao1_div_summary.xlsx"))
  
  #Alpha diversity plots
  alpha_list <- c("Chao1", "Shannon", "Simpson", "InvSimpson")
  # alpha_list <- c("Chao1", "Shannon")
  
  #Diversity plot
  if (is.null(pseq_div@sam_data$Donor) == F) {
    #With donors
    alpha_plot_with_donor <- plot_richness(pseq_div, metadata_class_no_formula, measures=alpha_list) + 
      theme(legend.position="none", text=element_text(size=20)) +
      theme(panel.background=element_rect(fill='white', colour='black')) + labs(x="") + 
      ggrepel:::geom_label_repel(aes(label=Donor), label.size=1, size=5) 
  } else {
    alpha_plot_with_donor <- ""
  }
  #Without donors
  alpha_plot_without_donor <- plot_richness(pseq_div, metadata_class_no_formula, measures=alpha_list) + 
    theme(legend.position="none", text=element_text(size=20)) +
    theme(panel.background=element_rect(fill='white', colour='black')) + labs(x="") 
  #Compare estimators
  # pairs(estimate_richness(pseq_div, measures=alpha_list)[,c(1,3,4,5)], cex.axis=3, cex.labels=3)
  pairs(estimate_richness(pseq_div, measures=alpha_list)[,c(1,3,4,5)], cex.axis=3, cex.labels=3)
  alpha_comparison_plot <- recordPlot()
  
  #Manually calculate alpha diversity and append gene count
  alpha_data <- alpha(pseq_div, index=c("Chao1", "Shannon", "diversity_gini_simpson", "diversity_inverse_simpson"))
  gene_count_data <- bind_rows(selected_gene_count[grep("count_data", names(selected_gene_count))])
  alpha_data <- alpha_data[rownames(gene_count_data),]
  alpha_data <- cbind(gene_count_data, alpha_data)
  #Get statistics of each estimator
  count_stats <- aggregate(alpha_data$Count, list(alpha_data[[metadata_class_no_formula]]), function(x) c(mean=mean(x), sd=sd(x), var=var(x), range=range(x)))
  chao1_stats <- aggregate(alpha_data$chao1, list(alpha_data[[metadata_class_no_formula]]), function(x) c(mean=mean(x), sd=sd(x), var=var(x), range=range(x)))
  shannon_stats <- aggregate(alpha_data$diversity_shannon, list(alpha_data[[metadata_class_no_formula]]), function(x) c(mean=mean(x), sd=sd(x), var=var(x), range=range(x)))
  simpson_stats <- aggregate(alpha_data$diversity_gini_simpson, list(alpha_data[[metadata_class_no_formula]]), function(x) c(mean=mean(x), sd=sd(x), var=var(x), range=range(x)))
  invsimpson_stats <- aggregate(alpha_data$diversity_inverse_simpson, list(alpha_data[[metadata_class_no_formula]]), function(x) c(mean=mean(x), sd=sd(x), var=var(x), range=range(x)))

  #Generate legends
  count_legend <- paste0(count_stats[,1], "\n", format(round(count_stats$x[,1],1), nsmall=1), " ± ", format(round(count_stats$x[,2],1), nsmall=1),
                         "\n", "Range: ", format(round(count_stats$x[,4],1), nsmall=1), " - ", format(round(count_stats$x[,5],1), nsmall=1))
  chao1_legend <- paste0(chao1_stats[,1], "\n", format(round(chao1_stats$x[,1],1), nsmall=1), " ± ", format(round(chao1_stats$x[,2],1), nsmall=1),
                         "\n", "Range: ", format(round(chao1_stats$x[,4],1), nsmall=1), " - ", format(round(chao1_stats$x[,5],1), nsmall=1))
  shannon_legend <- paste0(shannon_stats[,1], "\n", format(round(shannon_stats$x[,1],1), nsmall=1), " ± ", format(round(shannon_stats$x[,2],1), nsmall=1),
                           "\n", "Range: ", format(round(shannon_stats$x[,4],1), nsmall=1), " - ", format(round(shannon_stats$x[,5],1), nsmall=1))
  simpson_legend <- paste0(simpson_stats[,1], "\n", format(round(simpson_stats$x[,1],1), nsmall=1), " ± ", format(round(simpson_stats$x[,2],1), nsmall=1),
                           "\n", "Range: ", format(round(simpson_stats$x[,4],1), nsmall=1), " - ", format(round(simpson_stats$x[,5],1), nsmall=1))
  invsimpson_legend <- paste0(invsimpson_stats[,1], "\n", format(round(invsimpson_stats$x[,1],1), nsmall=1), " ± ", format(round(invsimpson_stats$x[,2],1), nsmall=1),
                           "\n", "Range: ", format(round(invsimpson_stats$x[,4],1), nsmall=1), " - ", format(round(invsimpson_stats$x[,5],1), nsmall=1))
  
  #Generate individual plots per estimator
  count_plot <- ggplot(alpha_data, aes(x=alpha_data[[metadata_class_no_formula]], y=Count)) + geom_point() + 
    theme(legend.position="none", text=element_text(size=30, face="bold"), panel.background=element_rect(fill='white', colour='black'), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + labs(x="", y="Gene count") +
    scale_x_discrete(labels=count_legend)
  chao1_plot <- ggplot(alpha_data, aes(x=alpha_data[[metadata_class_no_formula]], y=chao1)) + geom_point() + 
    theme(legend.position="none", text=element_text(size=30, face="bold"), panel.background=element_rect(fill='white', colour='black'), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + labs(x="", y="Chao1") +
    scale_x_discrete(labels=chao1_legend)
  shannon_plot <- ggplot(alpha_data, aes(x=alpha_data[[metadata_class_no_formula]], y=diversity_shannon)) + geom_point() + 
    theme(legend.position="none", text=element_text(size=30, face="bold"), panel.background=element_rect(fill='white', colour='black'), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + labs(x="", y="Shannon") +
    scale_x_discrete(labels=shannon_legend)
  simpson_plot <- ggplot(alpha_data, aes(x=alpha_data[[metadata_class_no_formula]], y=diversity_gini_simpson )) + geom_point() + 
    theme(legend.position="none", text=element_text(size=30, face="bold"), panel.background=element_rect(fill='white', colour='black'), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + labs(x="", y="Simpson") +
    scale_x_discrete(labels=simpson_legend)
  invsimpson_plot <- ggplot(alpha_data, aes(x=alpha_data[[metadata_class_no_formula]], y=diversity_inverse_simpson )) + geom_point() + 
    theme(legend.position="none", text=element_text(size=30, face="bold"), panel.background=element_rect(fill='white', colour='black'), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + labs(x="", y="InvSimpson") +
    scale_x_discrete(labels=invsimpson_legend)
  
  #Combine individual plots
  gridExtra:::grid.arrange(count_plot, chao1_plot, shannon_plot, simpson_plot, invsimpson_plot, ncol=5)
  alpha_count_plot <- recordPlot()
  
  #Return alpha diversity results as a list
  alpha_diversity_list <- list(global_div, div_summary, alpha_plot_with_donor, alpha_plot_without_donor, alpha_comparison_plot, alpha(pseq_div, index=c("Chao1", "Shannon", "diversity_gini_simpson", "diversity_inverse_simpson")), alpha_count_plot)
  names(alpha_diversity_list) <- c("global_div", "div_summary", "alpha_plot_with_donor", "alpha_plot_without_donor", "alpha_comparison_plot", "estimators", "alpha_count_plot")
  return(alpha_diversity_list)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to calculate beta-diversity in pseq objects contaning functional metagenomics data
pseq_beta_diversity <- function(pseq_taxa, pseq_families, pseq_pathways, metadata_class_no_formula, cluster_labels, variety_class_no_formula, selected_levels) {
  
  #In this version of the script, beta diversity will be calculated considering all microbiome data (not only the significant ones)
  #To perform diversity analysis on significant data, add:
  #significant_taxa, significant_genes, significant_path
  
  #Transform class to factor
  #Beta diversity labels
  pseq_taxa@sam_data[[metadata_class_no_formula]] <- factor(pseq_taxa@sam_data[[metadata_class_no_formula]], levels=selected_levels)
  pseq_families@sam_data[[metadata_class_no_formula]] <- factor(pseq_families@sam_data[[metadata_class_no_formula]], levels=selected_levels)
  pseq_pathways@sam_data[[metadata_class_no_formula]] <- factor(pseq_pathways@sam_data[[metadata_class_no_formula]], levels=selected_levels)
  #Cluster labels
  pseq_taxa@sam_data[[cluster_labels]] <- factor(pseq_taxa@sam_data[[cluster_labels]], levels=selected_levels)
  pseq_families@sam_data[[cluster_labels]] <- factor(pseq_families@sam_data[[cluster_labels]], levels=selected_levels)
  pseq_pathways@sam_data[[cluster_labels]] <- factor(pseq_pathways@sam_data[[cluster_labels]], levels=selected_levels)
  
  #TAXA
  #Prepare pseq object with significant taxa
  #significant_pseq_taxa <- prune_taxa(significant_taxa$taxon, pseq_taxa)
  significant_pseq_taxa <- pseq_taxa
  #Determine all levels for metadata class
  class_levels <- levels(significant_pseq_taxa@sam_data[[metadata_class_no_formula]])
  #Generate and empty dataframe to store mean values per class level
  diversity_levels <- data.frame(matrix(, nrow=nrow(significant_pseq_taxa@otu_table), ncol=0))
  #Subset significant class data to get the means of each group
  for (i in 1:length(class_levels)) {
    #Generate a pseq object for each level
    significant_i <- prune_samples(significant_pseq_taxa@sam_data[[metadata_class_no_formula]] == class_levels[i], significant_pseq_taxa)
    #Calculate beta-diversity for all functions
    i_diversity <- divergence(significant_i,
                              apply(abundances(significant_i), 1, median),
                              method = "bray")
    #Generate i box
    i_box <- data.frame(cbind.data.frame(significant_i@sam_data[[metadata_class_no_formula]], i_diversity))
    names(i_box) <- c(metadata_class_no_formula, "Diversity")
    #Append i box to common dataframe
    diversity_levels <- rbind(diversity_levels, i_box)
  }
  
  #Transform class to factor
  diversity_levels[[metadata_class_no_formula]] <- factor(diversity_levels[[metadata_class_no_formula]])
  
  #Transform diversity to numeric
  diversity_levels$Diversity <- as.numeric(diversity_levels$Diversity)
  
  #Get statistics
  beta_stats <- aggregate(diversity_levels$Diversity, list(diversity_levels[[metadata_class_no_formula]]), function(x) c(mean=mean(x), sd=sd(x), var=var(x), range=range(x)))
  #Generate legend
  beta_legend <- paste0(beta_stats[,1], "\n", format(round(beta_stats$x[,1],2), nsmall=2), " ± ", format(round(beta_stats$x[,2],2), nsmall=2),
                         "\n", "Range: ", format(round(beta_stats$x[,4],2), nsmall=2), " - ", format(round(beta_stats$x[,5],2), nsmall=2))
  
  #Boxplot beta diversity
  taxa_boxplot <- ggplot(diversity_levels, aes(x=.data[[metadata_class_no_formula]], y=Diversity, fill=diversity_levels[,1])) +
    geom_boxplot(alpha=0.3) +
    theme(legend.position="none", text=element_text(size=30), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
    xlab("") + ylab("Taxonomy")
  taxa_boxplot <- taxa_boxplot + 
    theme(panel.background=element_rect(fill='white', colour='black')) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 2)) + ylim(0, 1.0) + scale_x_discrete(labels=beta_legend)
  
  #Beta diversity - clustering
  mydist <- phyloseq::distance(significant_pseq_taxa, method="bray")
  myhclust <- try(hclust(mydist))
  if (class(myhclust) == "try-error") {
    mydist <- phyloseq::distance(significant_pseq_taxa, method="gower")
    myhclust <- hclust(mydist)
  }
  #Cluster showing classes
  myhclust$labels <- significant_pseq_taxa@sam_data[[cluster_labels]]
  plot(myhclust, main="Bray-Curtis Dissimilarity",
       xlab=metadata_class_no_formula, ylab="Beta diversity - Taxonomy", 
       cex=1.1, cex.axis=2, cex.lab=2, cex.main=2)
  
  
  #Plot cluster using ggplot2
  ggdendrogram(myhclust, rotate = TRUE, size = 4, theme_dendro = FALSE) + 
    theme(legend.position="none", text = element_text(size=25),
          panel.background=element_rect(fill='white', colour='black'),
          axis.text.y = element_text(size=11, face="italic")) +
    labs(x = "Bray-Curtis dissimilarity",
         y = metadata_class_no_formula) +
    ggtitle("Beta diversity - Taxonomy")
  
  
  #Generate an advanced dendogram plot
  dend <- as.dendrogram(myhclust)
  dend <- dend %>%
    # set("by_labels_branches_col", value=levels(labels(dend))[2]) %>% #Color branches by the presence of one sample from a given label
    # set("by_labels_branches_col", type = "any", value=levels(labels(dend))[1]) %>% #Color branches by the presence of samples from the same given label
    # set("by_labels_branches_lwd", value=levels(labels(dend))[2]) %>% #Set the line width of branches with specific labels
    # set("by_labels_branches_lty", value=levels(labels(dend))[2]) %>% #Set the line type of branches with specific labels
    set("nodes_pch", 19) %>% #Set the nodes' point type
    set("nodes_cex", 2) %>%  #Set the nodes' point size
    # color_branches(k=4) #Color by the number of clusters
    set("labels_cex", 2.1) %>% 
    set("branches_lwd", 3)
  
  #Set label colors by class 
  colors_to_use <- as.numeric(labels(dend))
  labels_colors(dend) <- colors_to_use
  #Set the same colors than in the rest of Figures
  labels_colors(dend) <- gsub(1, "#FDD6D3", labels_colors(dend))
  labels_colors(dend) <- gsub(2, "#D8E7B3", labels_colors(dend))
  labels_colors(dend) <- gsub(1, "#B3ECED", labels_colors(dend))
  labels_colors(dend) <- gsub(1, "#EED8FF", labels_colors(dend))
  #Plot dendogram
  plot(dend, main="Beta diversity - Taxonomy", ylab="Bray-Curtis dissimilarity",
       xlab="", cex.axis=1.5, cex.lab=1.5, cex.main=2, ylim=c(-0.1,1))
  
  taxa_cluster_plot <- recordPlot()
  
  #Export diversity values
  taxa_diversity <- diversity_levels
  
  
  #PCoA
  #Calculate strength of association
  #Assess differences through adonis and permanova test (show p-value)
  #Pick relative abundances (compositional) and sample metadata 
  pseq.rel <- microbiome::transform(significant_pseq_taxa, "compositional")
  otu <- abundances(pseq.rel)
  meta <- meta(pseq.rel)
  set.seed(600)
  permanova <- adonis(eval(parse(text=paste0("t(otu) ~ ",  metadata_class_no_formula))), data=meta, permutations=999, method="bray")
  permanova_p <- as.data.frame(permanova$aov.tab)[metadata_class_no_formula, "Pr(>F)"]
  permanova_coef <- coefficients(permanova)
  #Ordinate the data
  set.seed(600)
  pcoa_ord <- try(ordinate(significant_pseq_taxa, "PCoA", "bray"))
  if (class(pcoa_ord) == "try-error") {
    set.seed(600)
    pcoa_ord <- try(ordinate(significant_pseq_taxa, "PCoA", "gower"))
  }
  #Plot PCoA
  taxa_pcoa <- plot_ordination(significant_pseq_taxa, pcoa_ord, color = metadata_class_no_formula, shape = variety_class_no_formula) +
    geom_point(size = 8) + theme(text = element_text(size=35)) + ggtitle("PCoA - Taxonomy") + 
  geom_line(aes(group = significant_pseq_taxa@sam_data[[metadata_class_no_formula]]), color="grey",
            arrow = arrow(type = "closed", length=unit(0.075, "inches")))

  #Clean axis labels
  taxa_pcoa$labels$x <- gsub("Axis.", "PC ",taxa_pcoa$labels$x)
  taxa_pcoa$labels$y <- gsub("Axis.", "PC ",taxa_pcoa$labels$y)
  taxa_pcoa$labels$colour <- "Group"
  #Get variance percentages
  x_perc <- format(as.numeric(stringr::str_match(taxa_pcoa$labels$x, "\\[(.*?)%")[2]), nsmall=1)
  y_perc <- format(as.numeric(stringr::str_match(taxa_pcoa$labels$y, "\\[(.*?)%")[2]), nsmall=1)
  #Add formatted variance percentages
  taxa_pcoa$labels$x <- paste0(gsub("\\[(.*?)%]", paste0("(", x_perc, "%)"), taxa_pcoa$labels$x), " - PERMANOVA (p = " , permanova_p, ")")
  taxa_pcoa$labels$y <- gsub("\\[(.*?)%]", paste0("(", y_perc, "%)"), taxa_pcoa$labels$y)
  taxa_pcoa <- taxa_pcoa  + theme(panel.background=element_rect(fill='white', colour='black')) +
    scale_x_continuous(labels = function(x) format(x, nsmall = 2)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 2)) 
  #PCoA plot with donors
  if (is.null(significant_pseq_taxa@sam_data$Donor) == F) {
    taxa_pcoa_with_donor <- taxa_pcoa + ggrepel:::geom_label_repel(aes(label=Donor), label.size=1, size=5) 
  } else {
    taxa_pcoa_with_donor <- ""
  }
  
  
  
  
  #GENE FAMILIES
  #Prepare pseq object with significant gene families
  #significant_pseq_families <- prune_taxa(significant_genes$taxon, pseq_families)
  significant_pseq_families <- pseq_families
  #Determine all levels for metadata class
  class_levels <- levels(significant_pseq_families@sam_data[[metadata_class_no_formula]])
  #Generate and empty dataframe to store mean values per class level
  diversity_levels <- data.frame(matrix(, nrow=nrow(significant_pseq_families@otu_table), ncol=0))
  #Subset significant class data to get the means of each group
  for (i in 1:length(class_levels)) {
    #Generate a pseq object for each level
    significant_i <- prune_samples(significant_pseq_families@sam_data[[metadata_class_no_formula]] == class_levels[i], significant_pseq_families)
    #Calculate beta-diversity for all functions
    i_diversity <- divergence(significant_i,
                               apply(abundances(significant_i), 1, median),
                               method = "bray")
    #Generate i box
    i_box <- data.frame(cbind.data.frame(significant_i@sam_data[[metadata_class_no_formula]], i_diversity))
    names(i_box) <- c(metadata_class_no_formula, "Diversity")
    #Append i box to common dataframe
    diversity_levels <- rbind(diversity_levels, i_box)
  }

  #Transform class to factor
  diversity_levels[[metadata_class_no_formula]] <- factor(diversity_levels[[metadata_class_no_formula]])
  
  #Transform diversity to numeric
  diversity_levels$Diversity <- as.numeric(diversity_levels$Diversity)
  
  #Get statistics
  beta_stats <- aggregate(diversity_levels$Diversity, list(diversity_levels[[metadata_class_no_formula]]), function(x) c(mean=mean(x), sd=sd(x), var=var(x), range=range(x)))
  #Generate legend
  beta_legend <- paste0(beta_stats[,1], "\n", format(round(beta_stats$x[,1],2), nsmall=2), " ± ", format(round(beta_stats$x[,2],2), nsmall=2),
                        "\n", "Range: ", format(round(beta_stats$x[,4],2), nsmall=2), " - ", format(round(beta_stats$x[,5],2), nsmall=2))
  
  #Boxplot beta diversity
  families_boxplot <- ggplot(diversity_levels, aes(x=.data[[metadata_class_no_formula]], y=Diversity, fill=diversity_levels[,1])) +
                    geom_boxplot(alpha=0.3) +
    theme(legend.position="none", text=element_text(size=30), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
                    xlab("") + ylab("Gene families")
  families_boxplot <- families_boxplot + 
    theme(panel.background=element_rect(fill='white', colour='black')) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 2)) + ylim(0, 1.0) + scale_x_discrete(labels=beta_legend)
  
  #Beta diversity - clustering
  mydist <- phyloseq::distance(significant_pseq_families, method="bray")
  myhclust <- try(hclust(mydist))
  if (class(myhclust) == "try-error") {
    mydist <- phyloseq::distance(significant_pseq_families, method="gower")
    myhclust <- hclust(mydist)
  }
  #Cluster showing classes
  myhclust$labels <- significant_pseq_families@sam_data[[cluster_labels]]
  plot(myhclust, main="Bray-Curtis Dissimilarity",
                       xlab=metadata_class_no_formula, ylab="Beta diversity - Gene families", 
                       cex=1.1, cex.axis=2, cex.lab=2, cex.main=2)
  
  
  #Plot cluster using ggplot2
  ggdendrogram(myhclust, rotate = TRUE, size = 4, theme_dendro = FALSE) + 
    theme(legend.position="none", text = element_text(size=25),
          panel.background=element_rect(fill='white', colour='black'),
          axis.text.y = element_text(size=11, face="italic")) +
    labs(x = "Bray-Curtis dissimilarity",
         y = metadata_class_no_formula) +
    ggtitle("Beta diversity - Gene families")
  
  
  #Generate an advanced dendogram plot
  dend <- as.dendrogram(myhclust)
  dend <- dend %>%
    # set("by_labels_branches_col", value=levels(labels(dend))[2]) %>% #Color branches by the presence of one sample from a given label
    # set("by_labels_branches_col", type = "any", value=levels(labels(dend))[1]) %>% #Color branches by the presence of samples from the same given label
    # set("by_labels_branches_lwd", value=levels(labels(dend))[2]) %>% #Set the line width of branches with specific labels
    # set("by_labels_branches_lty", value=levels(labels(dend))[2]) %>% #Set the line type of branches with specific labels
    set("nodes_pch", 19) %>% #Set the nodes' point type
    set("nodes_cex", 2) %>%  #Set the nodes' point size
    # color_branches(k=4) #Color by the number of clusters
    set("labels_cex", 2.1) %>% 
    set("branches_lwd", 3)
  
  #Set label colors by class 
  colors_to_use <- as.numeric(labels(dend))
  labels_colors(dend) <- colors_to_use
  #Set the same colors than in the rest of Figures
  labels_colors(dend) <- gsub(1, "#FDD6D3", labels_colors(dend))
  labels_colors(dend) <- gsub(2, "#D8E7B3", labels_colors(dend))
  labels_colors(dend) <- gsub(1, "#B3ECED", labels_colors(dend))
  labels_colors(dend) <- gsub(1, "#EED8FF", labels_colors(dend))
  #Plot dendogram
  plot(dend, main="Beta diversity - Gene families", ylab="Bray-Curtis dissimilarity",
       xlab="", cex.axis=1.5, cex.lab=1.5, cex.main=2, ylim=c(-0.1,1))
  
  families_cluster_plot <- recordPlot()
  
  #Export diversity values
  families_diversity <- diversity_levels
  
  
  #PCoA
  #Calculate strength of association
  #Assess differences through adonis and permanova test (show p-value)
  #Pick relative abundances (compositional) and sample metadata 
  pseq.rel <- microbiome::transform(significant_pseq_families, "compositional")
  otu <- abundances(pseq.rel)
  meta <- meta(pseq.rel)
  set.seed(600)
  permanova <- adonis(eval(parse(text=paste0("t(otu) ~ ",  metadata_class_no_formula))), data=meta, permutations=999, method="bray")
  permanova_p <- as.data.frame(permanova$aov.tab)[metadata_class_no_formula, "Pr(>F)"]
  permanova_coef <- coefficients(permanova)
  #Ordinate the data
  set.seed(600)
  pcoa_ord <- try(ordinate(significant_pseq_families, "PCoA", "bray"))
  if (class(pcoa_ord) == "try-error") {
    set.seed(600)
    pcoa_ord <- try(ordinate(significant_pseq_families, "PCoA", "gower"))
  }
  #Plot PCoA
  families_pcoa <- plot_ordination(significant_pseq_families, pcoa_ord, color = metadata_class_no_formula, shape = variety_class_no_formula) +
    geom_point(size = 8) + theme(text = element_text(size=35)) + ggtitle("PCoA - Gene families") + 
  geom_line(aes(group = significant_pseq_families@sam_data[[metadata_class_no_formula]]), color="grey",
            arrow = arrow(type = "closed", length=unit(0.075, "inches")))
  
  #Clean axis labels
  families_pcoa$labels$x <- gsub("Axis.", "PC ",families_pcoa$labels$x)
  families_pcoa$labels$y <- gsub("Axis.", "PC ",families_pcoa$labels$y)
  families_pcoa$labels$colour <- "Group"
  #Get variance percentages
  x_perc <- format(as.numeric(stringr::str_match(families_pcoa$labels$x, "\\[(.*?)%")[2]), nsmall=1)
  y_perc <- format(as.numeric(stringr::str_match(families_pcoa$labels$y, "\\[(.*?)%")[2]), nsmall=1)
  #Add formatted variance percentages
  families_pcoa$labels$x <- paste0(gsub("\\[(.*?)%]", paste0("(", x_perc, "%)"), families_pcoa$labels$x), " - PERMANOVA (p = " , permanova_p, ")")
  families_pcoa$labels$y <- gsub("\\[(.*?)%]", paste0("(", y_perc, "%)"), families_pcoa$labels$y)
  families_pcoa <- families_pcoa  + theme(panel.background=element_rect(fill='white', colour='black')) +
    scale_x_continuous(labels = function(x) format(x, nsmall = 2)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 2))
  #PCoA plot with donors
  if (is.null(significant_pseq_families@sam_data$Donor) == F) {
    families_pcoa_with_donor <- families_pcoa + ggrepel:::geom_label_repel(aes(label=Donor), label.size=1, size=5) 
  } else {
    families_pcoa_with_donor <- ""
  }
  
  
  
  
  #GENE PATH
  #Prepare pseq object with significant gene path
  #significant_pseq_path <- prune_taxa(significant_path$taxon, pseq_pathways)
  significant_pseq_path <- pseq_pathways
  #Determine all levels for metadata class
  class_levels <- levels(significant_pseq_path@sam_data[[metadata_class_no_formula]])
  #Generate and empty dataframe to store mean values per class level
  diversity_levels <- data.frame(matrix(, nrow=nrow(significant_pseq_path@otu_table), ncol=0))
  #Subset significant class data to get the means of each group
  for (i in 1:length(class_levels)) {
    #Generate a pseq object for each level
    significant_i <- prune_samples(significant_pseq_path@sam_data[[metadata_class_no_formula]] == class_levels[i], significant_pseq_path)
    #Calculate beta-diversity for all functions
    i_diversity <- divergence(significant_i,
                              apply(abundances(significant_i), 1, median),
                              method = "bray")
    #Generate i box
    i_box <- data.frame(cbind.data.frame(significant_i@sam_data[[metadata_class_no_formula]], i_diversity))
    names(i_box) <- c(metadata_class_no_formula, "Diversity")
    #Append i box to common dataframe
    diversity_levels <- rbind(diversity_levels, i_box)
  }
  
  #Transform class to factor
  diversity_levels[[metadata_class_no_formula]] <- factor(diversity_levels[[metadata_class_no_formula]])
  
  #Transform diversity to numeric
  diversity_levels$Diversity <- as.numeric(diversity_levels$Diversity)
  
  #Get statistics
  beta_stats <- aggregate(diversity_levels$Diversity, list(diversity_levels[[metadata_class_no_formula]]), function(x) c(mean=mean(x), sd=sd(x), var=var(x), range=range(x)))
  #Generate legend
  beta_legend <- paste0(beta_stats[,1], "\n", format(round(beta_stats$x[,1],2), nsmall=2), " ± ", format(round(beta_stats$x[,2],2), nsmall=2),
                        "\n", "Range: ", format(round(beta_stats$x[,4],2), nsmall=2), " - ", format(round(beta_stats$x[,5],2), nsmall=2))
  
  #Boxplot beta diversity
  path_boxplot <- ggplot(diversity_levels, aes(x=.data[[metadata_class_no_formula]], y=Diversity, fill=diversity_levels[,1])) +
    geom_boxplot(alpha=0.3) +
    theme(legend.position="none", text=element_text(size=30), axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
    xlab("") + ylab("Pathways")
  path_boxplot <- path_boxplot + 
    theme(panel.background=element_rect(fill='white', colour='black')) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 2)) + ylim(0, 1.0) + scale_x_discrete(labels=beta_legend)
  
  #Beta diversity - clustering
  mydist <- phyloseq::distance(significant_pseq_path, method="bray")
  myhclust <- try(hclust(mydist))
  if (class(myhclust) == "try-error") {
    mydist <- phyloseq::distance(significant_pseq_path, method="gower")
    myhclust <- hclust(mydist)
  }
  #Cluster showing classes
  myhclust$labels <- significant_pseq_path@sam_data[[cluster_labels]]
  plot(myhclust, main="Bray-Curtis Dissimilarity",
       xlab=metadata_class_no_formula, ylab="Beta diversity - Metabolic pathways", 
       cex=1.1, cex.axis=2, cex.lab=2, cex.main=2)
  
  
  #Plot cluster using ggplot2
  ggdendrogram(myhclust, rotate = TRUE, size = 4, theme_dendro = FALSE) + 
    theme(legend.position="none", text = element_text(size=25),
          panel.background=element_rect(fill='white', colour='black'),
          axis.text.y = element_text(size=11, face="italic")) +
    labs(x = "Bray-Curtis dissimilarity",
         y = metadata_class_no_formula) +
    ggtitle("Beta diversity - Metabolic pathways")
  
  
  #Generate an advanced dendogram plot
  dend <- as.dendrogram(myhclust)
  dend <- dend %>%
    # set("by_labels_branches_col", value=levels(labels(dend))[2]) %>% #Color branches by the presence of one sample from a given label
    # set("by_labels_branches_col", type = "any", value=levels(labels(dend))[1]) %>% #Color branches by the presence of samples from the same given label
    # set("by_labels_branches_lwd", value=levels(labels(dend))[2]) %>% #Set the line width of branches with specific labels
    # set("by_labels_branches_lty", value=levels(labels(dend))[2]) %>% #Set the line type of branches with specific labels
    set("nodes_pch", 19) %>% #Set the nodes' point type
    set("nodes_cex", 2) %>%  #Set the nodes' point size
    # color_branches(k=4) #Color by the number of clusters
    set("labels_cex", 2.1) %>% 
    set("branches_lwd", 3)
  
  #Set label colors by class 
  colors_to_use <- as.numeric(labels(dend))
  labels_colors(dend) <- colors_to_use
  #Set the same colors than in the rest of Figures
  labels_colors(dend) <- gsub(1, "#FDD6D3", labels_colors(dend))
  labels_colors(dend) <- gsub(2, "#D8E7B3", labels_colors(dend))
  labels_colors(dend) <- gsub(1, "#B3ECED", labels_colors(dend))
  labels_colors(dend) <- gsub(1, "#EED8FF", labels_colors(dend))
  #Plot dendogram
  plot(dend, main="Metabolic pathways", ylab="Bray-Curtis dissimilarity",
       xlab="", cex.axis=1.5, cex.lab=1.5, cex.main=2, ylim=c(-0.1,1))
  
  path_cluster_plot <- recordPlot()
    
  #Export diversity values
  path_diversity <- diversity_levels
  
  
  #PCoA
  #Calculate strength of association
  #Assess differences through adonis and permanova test (show p-value)
  #Pick relative abundances (compositional) and sample metadata 
  pseq.rel <- microbiome::transform(significant_pseq_families, "compositional")
  otu <- abundances(pseq.rel)
  meta <- meta(pseq.rel)
  set.seed(600)
  permanova <- adonis(eval(parse(text=paste0("t(otu) ~ ",  metadata_class_no_formula))), data=meta, permutations=999, method="bray")
  permanova_p <- as.data.frame(permanova$aov.tab)[metadata_class_no_formula, "Pr(>F)"]
  permanova_coef <- coefficients(permanova)
  #Ordinate the data
  set.seed(600)
  pcoa_ord <- try(ordinate(significant_pseq_path, "PCoA", "bray"))
  if (class(pcoa_ord) == "try-error") {
    set.seed(600)
    pcoa_ord <- try(ordinate(significant_pseq_path, "PCoA", "gower"))
  }
  #Plot PCoA
  path_pcoa <- plot_ordination(significant_pseq_path, pcoa_ord, color = metadata_class_no_formula, shape = variety_class_no_formula) +
    geom_point(size = 8) + theme(text = element_text(size=35)) + ggtitle("PCoA - Metabolic pathways") + 
  geom_line(aes(group = significant_pseq_path@sam_data[[metadata_class_no_formula]]), color="grey",
            arrow = arrow(type = "closed", length=unit(0.075, "inches")))

  #Clean axis labels
  path_pcoa$labels$x <- gsub("Axis.", "PC ",path_pcoa$labels$x)
  path_pcoa$labels$y <- gsub("Axis.", "PC ",path_pcoa$labels$y)
  path_pcoa$labels$colour <- "Group"
  #Get variance percentages
  x_perc <- format(as.numeric(stringr::str_match(path_pcoa$labels$x, "\\[(.*?)%")[2]), nsmall=1)
  y_perc <- format(as.numeric(stringr::str_match(path_pcoa$labels$y, "\\[(.*?)%")[2]), nsmall=1)
  #Add formatted variance percentages
  path_pcoa$labels$x <- paste0(gsub("\\[(.*?)%]", paste0("(", x_perc, "%)"), path_pcoa$labels$x), " - PERMANOVA (p = " , permanova_p, ")")
  path_pcoa$labels$y <- gsub("\\[(.*?)%]", paste0("(", y_perc, "%)"), path_pcoa$labels$y)
  path_pcoa <- path_pcoa  + theme(panel.background=element_rect(fill='white', colour='black')) +
    scale_x_continuous(labels = function(x) format(x, nsmall = 2)) +
    scale_y_continuous(labels = function(x) format(x, nsmall = 2))
  #PCoA plot with donors
  if (is.null(significant_pseq_path@sam_data$Donor) == F) {
    path_pcoa_with_donor <- path_pcoa + ggrepel:::geom_label_repel(aes(label=Donor), label.size=1, size=5) 
  } else {
    path_pcoa_with_donor <- ""
  }
  
  
  
  
  #Return beta diversity results as a list
  beta_diversity_list <- list(taxa_boxplot, taxa_cluster_plot, 
                              families_boxplot, families_cluster_plot, 
                              path_boxplot, path_cluster_plot,
                              taxa_diversity, families_diversity, path_diversity,
                              taxa_pcoa, families_pcoa, path_pcoa,
                              taxa_pcoa_with_donor, families_pcoa_with_donor, path_pcoa_with_donor)
  names(beta_diversity_list) <- c("Taxa_boxplot", "Taxa_cluster", 
                                  "Families_boxplot", "Families_cluster", 
                                  "Path_boxplot", "Path_cluster_plot",
                                  "Taxa_diversity", "Families_diversity", "Path_diversity",
                                  "Taxa_PCoA", "Families_PCoA", "Path_PCoA",
                                  "Taxa_PCoA_with_donor", "Families_PCoA_with_donor", "Path_PCoA_with_donor")
  return(beta_diversity_list)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
