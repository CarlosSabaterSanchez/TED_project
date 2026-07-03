
#Reading Data, check getwd, header = TRUE tells the name of the variables
getwd()
library(microbiome) 
library(phyloseq)
library(plyr) 
library(dplyr) 
library(DESeq2) 
library(microbiomeMarker) 
library(Maaslin2)
library(ggdendro) 
library(dendextend) 
library(RColorBrewer) 
library(reshape2) 
library(caret)
library(factoextra) 
library(randomForest) 
library(NeuralNetTools) 
library(corrr) 
library(smotefamily) 

#Call custom functions to generate pseq objects from metagenomics functional analyses
source("Input/Functions/1_Pseq_object_kraken_humann.R")
source("Input/Functions/2_Pseq_stats_kraken_humann.R")
source("Input/Functions/3_Pseq_gene_count.R")
source("Input/Functions/4_Pseq_beta_tax_functional.R")
source("Input/Functions/5_Pseq_barplot_tax_functional.R")
source("Input/Functions/6_Statistical_differences.R")
source("Input/Functions/7_ML.R")
source("Input/Functions/2_Longitudinal_study.R")

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Open a new window to store plots
x11()

#Read metadata
metadata <- openxlsx::read.xlsx("Input/Data/Metadata.xlsx", sheet=1)
rownames(metadata) <- metadata$Shotgun_code

#Read kraken taxonomic assignments
family_kraken <- read.table(file="Input/Data/kraken_family.tsv", sep ='\t', header=T, check.names=F)
family_maxikraken <- read.table(file="Input/Data/maxikraken_family.tsv", sep ='\t', header=T, check.names=F)
genus_kraken <- read.table(file="Input/Data/kraken_genus.tsv", sep ='\t', header=T, check.names=F)
genus_maxikraken <- read.table(file="Input/Data/maxikraken_genus.tsv", sep ='\t', header=T, check.names=F)
species_kraken <- read.table(file="Input/Data/kraken_species.tsv", sep ='\t', header=T, check.names=F)
species_maxikraken <- read.table(file="Input/Data/maxikraken_species.tsv", sep ='\t', header=T, check.names=F)
strain_kraken <- openxlsx::read.xlsx("Input/Data/Metadata.xlsx", sheet=4)
strain_maxikraken <- openxlsx::read.xlsx("Input/Data/Metadata.xlsx", sheet=3)

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
pseq_list <- generate_pseq_kraken_humann(pseq_input, selected_metadata, abundance_filter, frequency_filter)
#Save R object
saveRDS(pseq_list, "R_objects/pseq_list.rds")
#Read R object (takes less than 1 min)
pseq_list <- readRDS("R_objects/pseq_list.rds")

#Check humann taxonomy abundance data##################################################################
pseq_composition <- aggregate_taxa(pseq_list$humann_taxa_abundance, "Species")
#Remove unknown taxa
pseq_composition <- prune_taxa(grep("unknown", rownames(pseq_composition@otu_table),
                                    ignore.case=T, value=T, invert=T), pseq_composition)
#Remove duplicated taxa generated by aggregate
#Remove unknown taxa
pseq_composition <- prune_taxa(grep("[|]", rownames(pseq_composition@otu_table),
                                    ignore.case=T, value=T, invert=T), pseq_composition)
#Inspect colSums
colSums(pseq_composition@otu_table)

#Get those species determined by both humann and kraken#################################################
select_humann <- aggregate_taxa(pseq_list$humann_taxa_counts_filtered, "Species")
select_kraken <- pseq_list$kraken_species_counts_filtered
#Find common taxa-
kraken_and_human_species <- intersect(gsub("_", " ", rownames(select_humann@otu_table)),
                                      rownames(select_kraken@otu_table))
#Get those species determined by humann or kraken
only_humann_species <- setdiff(gsub("_", " ", rownames(select_humann@otu_table)), kraken_and_human_species)
only_kraken_species <- setdiff(rownames(select_kraken@otu_table), kraken_and_human_species)

#Get all taxonomic clades determined
all_clades <- sort(unique(c(
  rownames(pseq_list$kraken_family_counts_filtered@otu_table),
  rownames(pseq_list$kraken_genus_counts_filtered@otu_table),
  rownames(pseq_list$kraken_species_counts_filtered@otu_table),
  rownames(pseq_list$kraken_strain_counts_filtered@otu_table),
  gsub("_", " ", gsub(".*remove", "", sub("_", "remove",  rownames(pseq_list$humann_taxa_counts_filtered@otu_table))))
)))

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Statistical analyses:
#Select filtered pseq objects
pseq_list_filtered <- pseq_list[grep("filtered", names(pseq_list))]
pseq_list_time <- list()

#Subset Time
times <- unique(pseq_list_filtered$kraken_family_counts_filtered@sam_data$Time)
#Iterate pseq list
for (i in 1:length(pseq_list_filtered)) {
  i_pseq <- pseq_list_filtered[[i]]
  #Iterate manufacturers
  for (j in 1:length(times)) {
    j_time <- times[j]
    #Subset i and j
    i_j_pseq <- subset_samples(i_pseq, i_pseq@sam_data$Time == j_time)
    #Append i and j subset to a single list
    pseq_list_time[[length(pseq_list_time)+1]] <- i_j_pseq
    names(pseq_list_time)[[length(pseq_list_time)]] <- paste0(names(pseq_list_filtered)[i], "_", j_time)
  }
}

#Generate the final subset list
subset_list <- list()
for (i in 1:length(times)) {
  subset_list[[length(subset_list)+1]] <- pseq_list_time[grep(paste0("_",times[i]), names(pseq_list_time))]
  names(subset_list)[[length(subset_list)]] <- paste0("Time_", times[i])
}

#Microbiota statistical analyses
subset_list <- list(pseq_list_filtered, pseq_list_filtered) 
names(subset_list) <- c("Time", "Substrate")
#Prepare a vector of factors to be used as categorical variables in the statistical analyses
categorical_var <- c("Time", "Substrate")
#Run multiple statistical tests
increment_threshold <- 10
p_threshold <- 0.05
padj_threshold <- 0.25
#Run the custom function
multiple_microbiota_stats <- microbiota_multiple_tests(subset_list, categorical_var, 10, 0.05, 0.25, all_clades)
#Save R object
saveRDS(multiple_microbiota_stats, "R_objects/multiple_microbiota_stats.rds")
#Read R object (takes less than 1 min)
multiple_microbiota_stats <- readRDS("R_objects/multiple_microbiota_stats.rds")

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Generate a composition barplot
barplot_time_0 <- abundance_barplot(pseq_list_time$kraken_species_percentage_filtered_0, c("Species"), "Substrate", unique(pseq_list_time$kraken_species_percentage_filtered_0@sam_data$Substrate))

#Plot composition data
win.metafile(width=15, height=7) #Must be full screen
ggsave(file="Figures/Origin_barplot_time_0.png", barplot_time_0$composition_barplot, device='png', dpi=600) #+ expand_limits(x=c(-4, 4)) #Saves g
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Generate pseq objects
#abundance_filter: 0.1 (%)
#frequency_filter: 0.5 (0-1)
abundance_filter <- 0.45
frequency_filter <- 0.13
selected_metadata <-  strain_taxa$selected_metadata
pseq_list <- generate_pseq_kraken_humann(pseq_input, selected_metadata, abundance_filter, frequency_filter)

#Abundance analysis
#Prepare inputs
#Select filtered pseq objects
pseq_list_filtered <- pseq_list[grep("filtered", names(pseq_list))]
pseq_list_substrate <- list()

#Subset Substrate
substrates <- unique(pseq_list_filtered$kraken_family_counts_filtered@sam_data$Substrate)
#Iterate pseq list
for (i in 1:length(pseq_list_filtered)) {
  i_pseq <- pseq_list_filtered[[i]]
  #Iterate substrates
  for (j in 1:length(substrates)) {
    j_substrate <- substrates[j]
    #Subset i and j
    i_j_pseq <- subset_samples(i_pseq, i_pseq@sam_data$Substrate == j_substrate)
    #Append i and j subset to a single list
    pseq_list_substrate[[length(pseq_list_substrate)+1]] <- i_j_pseq
    names(pseq_list_substrate)[[length(pseq_list_substrate)]] <- paste0(names(pseq_list_filtered)[i], "_", j_substrate)
  }
}

#Barley barplots
barplot_BSG_U <- abundance_barplot(pseq_list_substrate$`kraken_family_counts_filtered_BSG-U`, c("Species"), "Time", 
                                    unique(pseq_list_substrate$`kraken_species_percentage_filtered_BSG-U`@sam_data$Time))

barplot_BSG_Vz <- abundance_barplot(pseq_list_substrate$`kraken_species_percentage_filtered_BSG-Vz`, c("Species"), "Time", 
                                unique(pseq_list_substrate$`kraken_species_percentage_filtered_BSG-Vz`@sam_data$Time))



#Artichocke barplots
barplot_FAB <- abundance_barplot(pseq_list_substrate$kraken_species_percentage_filtered_FAB, c("Species"), "Time", 
                                unique(pseq_list_substrate$kraken_species_percentage_filtered_FAB@sam_data$Time))

barplot_BAB <- abundance_barplot(pseq_list_substrate$kraken_species_percentage_filtered_BAB, c("Species"), "Time", 
                                unique(pseq_list_substrate$kraken_species_percentage_filtered_BAB@sam_data$Time))
                    
#Basal samples (Time 0)
#Prepare inputs
#Select filtered pseq objects
pseq_list_time <- list()

#Subset time
times <- unique(pseq_list_filtered$kraken_family_counts_filtered@sam_data$Time)
#Iterate pseq list
for (i in 1:length(pseq_list_filtered)) {
  i_pseq <- pseq_list_filtered[[i]]
  #Iterate substrates
  for (j in 1:length(times)) {
    j_time <- times[j]
    #Subset i and j
    i_j_pseq <- subset_samples(i_pseq, i_pseq@sam_data$Time == j_time)
    #Append i and j subset to a single list
    pseq_list_time[[length(pseq_list_time)+1]] <- i_j_pseq
    names(pseq_list_time)[[length(pseq_list_time)]] <- paste0(names(pseq_list_filtered)[i], "_", j_time)
  }
}

#basal barplots
barplot_basal_samples <- abundance_barplot(pseq_list_time$kraken_species_percentage_filtered_0, c("Species"), "Substrate", 
                                    unique(pseq_list_time$kraken_species_percentage_filtered_0@sam_data$Substrate))

#Plot basal abundance data
win.metafile(width=15, height=7) #Must be full screen
ggsave(file="Figures/Abundance_basal_samples.png", barplot_basal_samples$composition_barplot, device='png', dpi=600) #+ expand_limits(x=c(-4, 4)) #Saves g
dev.off()

#Final barplots
barplot_final_samples <- abundance_barplot(pseq_list_time$kraken_species_percentage_filtered_72, c("Species"), "Substrate", 
                                           unique(pseq_list_time$kraken_species_percentage_filtered_72@sam_data$Substrate))

#Plot abundance data
win.metafile(width=15, height=7) #Must be full screen
ggsave(file="Figures/Abundance_final_samples.png", barplot_final_samples$composition_barplot, device='png', dpi=600) #+ expand_limits(x=c(-4, 4)) #Saves g
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

library(phyloseq)
library(ggplot2)

df <- psmelt(pseq_BSG)

ggplot(df, aes(x = Sample, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Substrate, scales = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Genus level
#basal barplots
barplot_basal_samples_genus <- abundance_barplot(pseq_list_time$kraken_genus_percentage_filtered_0, c("Species"), "Substrate", 
                                           unique(pseq_list_time$kraken_genus_percentage_filtered_0@sam_data$Substrate))

#Plot basal abundance data
win.metafile(width=15, height=7) #Must be full screen
ggsave(file="Figures/Abundance_basal_samples_genus.png", barplot_basal_samples_genus$composition_barplot, device='png', dpi=600) #+ expand_limits(x=c(-4, 4)) #Saves g
dev.off()

#Final barplots
barplot_final_samples_genus <- abundance_barplot(pseq_list_time$kraken_genus_percentage_filtered_72, c("Species"), "Substrate", 
                                           unique(pseq_list_time$kraken_genus_percentage_filtered_72@sam_data$Substrate))

#Plot abundance data
win.metafile(width=15, height=7) #Must be full screen
ggsave(file="Figures/Abundance_final_samples_genus.png", barplot_final_samples_genus$composition_barplot, device='png', dpi=600) #+ expand_limits(x=c(-4, 4)) #Saves g
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Finish
