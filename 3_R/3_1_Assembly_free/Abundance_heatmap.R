
#LONGITUDINAL MICROBIOME ANALYSIS
#update R
#library(installr)
#updateR()

#Reading Data, check getwd, header = TRUE tells the name of the variables
getwd()
library(microbiome)
library(phyloseq)
library(ggplot2)
library(vegan)
library(caret)
library(Rtsne)
library(plyr)
library(qvalue)
library(Maaslin2)
library(RColorBrewer)
library(ggdendro)
library(dendextend)
library(rcompanion)
library(FSA)
library(dplyr)
library(DESeq2)
library(microbiomeMarker)

#Call custom functions to generate pseq objects from metagenomics functional analyses
source("Input/Functions/1_Pseq_object_kraken_humann.R")
source("Input/Functions/2_Longitudinal_study.R")

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#PREPARE INPUT DATA
#Open a new window to store plots
x11()

#Read metadata
metadata <- openxlsx::read.xlsx("Input/Data/Metadata_TED.xlsx", sheet=1)
rownames(metadata) <- metadata$Shotgun_code

#Read kraken taxonomic assignments
family_kraken <- read.table(file="Input/Data/kraken_family.tsv", sep ='\t', header=T, check.names=F)
family_maxikraken <- read.table(file="Input/Data/maxikraken_family.tsv", sep ='\t', header=T, check.names=F)
genus_kraken <- read.table(file="Input/Data/kraken_genus.tsv", sep ='\t', header=T, check.names=F)
genus_maxikraken <- read.table(file="Input/Data/maxikraken_genus.tsv", sep ='\t', header=T, check.names=F)
species_kraken <- read.table(file="Input/Data/kraken_species.tsv", sep ='\t', header=T, check.names=F)
species_maxikraken <- read.table(file="Input/Data/maxikraken_species.tsv", sep ='\t', header=T, check.names=F)
strain_kraken <- openxlsx::read.xlsx("Input/Data/Metadata_TED.xlsx", sheet=4)
strain_maxikraken <- openxlsx::read.xlsx("Input/Data/Metadata_TED.xlsx", sheet=3)

#Combine kraken outputs
family_taxa <- combine_kraken_taxa(family_kraken, family_maxikraken, metadata)
genus_taxa <- combine_kraken_taxa(genus_kraken, genus_maxikraken, metadata)
species_taxa <- combine_kraken_taxa(species_kraken, species_maxikraken, metadata)
strain_taxa <- combine_kraken_taxa(strain_kraken, strain_maxikraken, metadata)

#Read humann taxonomic assignments
humann_taxa <- read.table(file="Input/Data/metaphlan3_merged.txt", sep ='\t', header=T, check.names=F)
#Set taxonomic clades as rownames
rownames(humann_taxa) <- humann_taxa$clade_name
#Remove first two columns (redundant taxonomy and NCBI codes)
humann_taxa <- humann_taxa[,-c(1,2)]
#Transpose taxonomy table to match metadata format
humann_taxa <- t(humann_taxa)
#Order samples
humann_taxa <- humann_taxa[rownames(strain_taxa$selected_metadata),]

#Read humann metabolic pathways (takes less than 1 min)
humann_pathways <- read.csv("Input/Data/humann3_pathabundance_cpm_stratified.tsv", row.names=1, header=T, sep="\t", check.names=F)
colnames(humann_pathways) <- gsub("_.*", "", colnames(humann_pathways))
#Order samples
humann_pathways <- humann_pathways[,rownames(strain_taxa$selected_metadata)]

#Read gene families (takes approximately 15 min)
humann_families <- read.csv("Input/Data/humann3_genefamilies_cpm_stratified.tsv", row.names=1, header=T, sep="\t", check.names=F)
colnames(humann_families) <- gsub("_Abundance-RPKs", "", colnames(humann_families))
#Order samples
humann_families <- humann_families[,rownames(strain_taxa$selected_metadata)]
#Save R object
saveRDS(humann_families, "R_objects/humann_families.rds")
#Read R object (takes less than 1 min)
humann_families <- readRDS("R_objects/humann_families.rds")

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Prepare a list of pseq inputs
pseq_input <- list(family_taxa$taxa_counts, family_taxa$taxa_percentage,
                   genus_taxa$taxa_counts, genus_taxa$taxa_percentage,
                   species_taxa$taxa_counts, species_taxa$taxa_percentage,
                   strain_taxa$taxa_counts, strain_taxa$taxa_percentage,
                   humann_taxa, #humann_families, 
                   humann_pathways)
names(pseq_input) <- c("kraken_family_counts", "kraken_family_percentage",
                       "kraken_genus_counts", "kraken_genus_percentage",
                       "kraken_species_counts", "kraken_species_percentage",
                       "kraken_strain_counts", "kraken_strain_percentage",
                       "humann_taxa", #"humann_func_families", 
                       "humann_func_pathways")

#Rename sample codes
rownames(strain_taxa$selected_metadata) <- gsub("-", "", rownames(strain_taxa$selected_metadata))
strain_taxa$selected_metadata$Shotgun_code <- gsub("-", "", strain_taxa$selected_metadata$Shotgun_code)
for (i in 1:length(pseq_input)) {
  if (length(grep("humann_func", names(pseq_input)[i])) != 0) {
    colnames(pseq_input[[i]]) <- gsub("-", "", colnames(pseq_input[[i]]))
  } else {
    rownames(pseq_input[[i]]) <- gsub("-", "", rownames(pseq_input[[i]]))
  }
}

#Generate pseq objects
#abundance_filter: 0.1 (%)
#frequency_filter: 0.5 (0-1)
abundance_filter <- 0.1
frequency_filter <- 0.5
selected_metadata <-  strain_taxa$selected_metadata

#Format sample codes
colnames(pseq_input$humann_func_pathways) <- paste0("M", colnames(pseq_input$humann_func_pathways))
rownames(pseq_input$kraken_family_counts) <- paste0("M", rownames(pseq_input$kraken_family_counts))
rownames(pseq_input$kraken_family_percentage) <- paste0("M", rownames(pseq_input$kraken_family_percentage))
rownames(pseq_input$kraken_genus_counts) <- paste0("M", rownames(pseq_input$kraken_genus_counts))
rownames(pseq_input$kraken_genus_percentage) <- paste0("M", rownames(pseq_input$kraken_genus_percentage))
rownames(pseq_input$kraken_species_counts) <- paste0("M", rownames(pseq_input$kraken_species_counts))
rownames(pseq_input$kraken_species_percentage) <- paste0("M", rownames(pseq_input$kraken_species_percentage))
rownames(pseq_input$kraken_strain_counts) <- paste0("M", rownames(pseq_input$kraken_strain_counts))
rownames(pseq_input$kraken_strain_percentage) <- paste0("M", rownames(pseq_input$kraken_strain_percentage))
rownames(pseq_input$humann_taxa) <- paste0("M", rownames(pseq_input$humann_taxa))
rownames(selected_metadata) <- paste0("M", rownames(selected_metadata))
selected_metadata$Shotgun_code <- paste0("M", selected_metadata$Shotgun_code)

#MICROBIOTA STATISTICAL ANALYSIS
#Generate pseq objects of metagenomic functional data
pseq_list <- generate_pseq_kraken_humann(pseq_input, selected_metadata, abundance_filter, frequency_filter)
#Save R object
saveRDS(pseq_list, "R_objects/pseq_list.rds")
#Read R object (takes less than 1 min)
pseq_list <- readRDS("R_objects/pseq_list.rds")

#Subset pseqs list
pseq_taxa <- pseq_list$kraken_species_counts_filtered
pseq_families <- pseq_list$kraken_family_counts_filtered 
pseq_pathways <- pseq_list$humann_func_pathways_filtered

#Load additional libraries
library(ANCOMBC)
library(dplyr)

#Perform longitudinal study
#Select labels for statistical analysis
health_label <- names(pseq_taxa@sam_data)[3]
donor_label <- names(pseq_taxa@sam_data)[6]
time_label <- names(pseq_taxa@sam_data)[5]
substrate_label <- names(pseq_taxa@sam_data)[4]
#Perform longitudinal study using a custom function developed for this purpose
p_value_threshold <- 0.05
selected_pseq <- pseq_taxa
longitudinal <- longitudinal_study_trends(selected_pseq, substrate_label, health_label, time_label, p_value_threshold)
#Export the results of the longitudinal table as excel sheets
xlsx::write.xlsx(longitudinal$significant_all_per_substrate, "Results/longitudinal_tables.xlsx")
#Save R objects
saveRDS(longitudinal, file = "R_objects/longitudinal.rds")

####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3
####################################################################################################################################3

#Read R objects
longitudinal <- readRDS("R_objects/longitudinal.rds")
#Heatmap
#Select the most interesting results
selected_longitudinal <- longitudinal$significant_all_per_substrate
#Determine increment quantiles in order to filter low increments
increment_df <- data.frame(selected_longitudinal[,3:ncol(selected_longitudinal)])
increment_values <- as.data.frame(apply(increment_df, 1, function(x) {gsub(".*Increment: ", "", x)}))
increment_values <- as.data.frame(apply(increment_values, 1, function(x) {gsub(" .*", "", x)}))
increment_values <- as.data.frame(apply(increment_values, 1, function(x) {gsub("-", "0", x)}))
increment_values <- increment_values %>% mutate_all(~ as.numeric(as.character(.)))
quantile(apply(increment_values,1,max))
#Filter increments lower than 5
for (i in 3:ncol(selected_longitudinal)) {
  for (j in 1:nrow(selected_longitudinal)) {
    #Get j i cell
    j_i_cell <- selected_longitudinal[j,i]
    if (length(grep("Increment", j_i_cell))) {
      #Get j i increment
      j_i_increment <- as.numeric(gsub(".*Increment: (.+) \\+.*", "\\1", j_i_cell))
      #Filter increments
      if (j_i_increment < quantile(apply(increment_values,1,max))[2]) {
        selected_longitudinal[j,i] <- "-"
      }
    }
  }
}
#Filter empty columns
remove_cols <- c()
for (i in 3:ncol(selected_longitudinal)) {
  if (length(grep("-", t(selected_longitudinal[,i]))) == nrow(selected_longitudinal)) {
    remove_cols <- c(remove_cols,i)
  }
}
#Remove empty columns
if (is.null(remove_cols) != T) {
  selected_longitudinal <- selected_longitudinal[,-remove_cols]
}

#Sort substrates for heatmap plotting
substrate_order <- c("FAB", "BAB",
                     "BSG-U", "BSG-Vz")
selected_longitudinal <- selected_longitudinal[order(match(selected_longitudinal$Substrate, substrate_order)),]
#Run longitudinal heatmaps
selected_heatmaps <- longitudinal_heatmap(selected_pseq, selected_longitudinal,
                                          time_label, substrate_label, health_label) 
#Low increment heatmap
#Generate a high quality Figure
win.metafile(width=14, height=7) #Must be full screen
png("Figures/low_heatmap.png", width=20, height=10, units="in", res=600, pointsize=50)
selected_heatmaps$low_plot_record
dev.off()
#Medium and high increment heatmaps
#Generate a high quality Figure
win.metafile(width=14, height=7) #Must be full screen
png("Figures/medium_high_heatmap.png", width=20, height=10, units="in", res=600, pointsize=50)
selected_heatmaps$high_medium_plot_record
dev.off()

#Alphabetic order
selected_longitudinal <- selected_longitudinal[,c(1,2, order(sub(".*?_", "", names(selected_longitudinal)[3:length(selected_longitudinal)]))+2)]
#Run a single longitudinal heatmap
selected_heatmaps <- longitudinal_single_heatmap(selected_pseq, selected_longitudinal,
                                          time_label, substrate_label, health_label) 

#Generate a high quality Figure
win.metafile(width=14, height=7) #Must be full screen
png("Figures/single_heatmap.png", width=20, height=10, units="in", res=600, pointsize=50)
selected_heatmaps$plot
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Finish
