
#Reading Data, check getwd, header = TRUE tells the name of the variables
getwd()

#Load libraries
library(readr)
library(xml2)
library(htmltools)
library(stringr)
library(purrr)
library(openxlsx)
library(readxl)
library(gtools)
library(dplyr)
library(stringi)
#Extra libraries to mix italics and normal fonts in Figure axis
library(ggtree)
library(ggtext)
#Extra libraries to automatically generate heatmaps with dendograms
library(ggalign)
#Extra libraries to plot correlations
library(ggpubr)

#Invoke custom functions 
source("Functions/1_Plots.R")
source("Functions/3_CAZy_domains.R")

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Select the file list and coverage threshold
file_list <- list.files(path="cazy_functional_domains/")
coverage_threshold <- 0.8
#Remove (or not) microbial genomes not showing any domain of interest
remove_taxa_without_domains <- F
#Select latest dbCAN substrate metadata (does not contain GT domains and cohesin)
latest_dbcan_metadata <- "https://bcb.unl.edu/dbCAN2/download/dbsub_data/Substrate_05-03-2022.txt"
#Inspect CAZy substrates available
cazy_metadata <- read.csv(latest_dbcan_metadata, sep="\t", skip=1, encoding="UTF-8")
unique(cazy_metadata$Substrate_Simple_V2)
#Select the search terms of interest for each substrate
#Pectin
# terms_interest <- c("pecti", "rhamnogalacturonan", "galacturon")
# exclude_terms <- c("amylopectin")
#Xylan/Arabinoxylan
# terms_interest <- c("xylan", "xylooligosaccharide")
# exclude_terms <- NULL
#Mannan
# terms_interest <- c("mann")
# exclude_terms <- c("mannosylglycerate", "alpha-D-mannoside", "beta-D-mannoside")
#Cellulose
# terms_interest <- c("cellul", "cello")
# exclude_terms <- NULL
#Lignin
# terms_interest <- c("lign")
# exclude_terms <- NULL
#Glucan
# terms_interest <- c("glucan")
# exclude_terms <- NULL
#Fructan
# terms_interest <- c("fructan", "oligofructose", "fructooligosaccharide", "fructo-oligosaccharide", "inulin")
# exclude_terms <- NULL

#Format strain and variable labels
#Read metadata file
metadata <- openxlsx::read.xlsx("Metadata.xlsx", sheet=2)
rownames(metadata) <- metadata$ID

#Subset artichoke samples
metadata <- metadata[(metadata$Substrate=="BSG-U") | (metadata$Substrate=="BSG-Vz"), ]
file_list <- grep(paste0("^", metadata$ID, "_hmmer", collapse="|"), file_list, value=T)

#Get strain sequences
strain_seq <- metadata$ID
#Get dendogram class
dendo_class <- "Substrate_Time"

# #Add abbreviated names to the metada file
# #Add strain abbreviations
metadata$Taxonomy_metadata <- gsub("_", " ", gsub("_.? ", " ", gsub(".*__", "", gsub(";s__$", "", metadata$Taxonomy_metadata))))
metadata$Taxonomy_metadata <- sub("(\\S+ \\S+) [A-Za-z]+(?!\\S)", "\\1", metadata$Taxonomy_metadata, perl=TRUE)
metadata$Taxonomy_metadata <- paste0(metadata$Taxonomy_metadata, " MAG ", metadata$MAG, " ", metadata$Substrate, " ", metadata$Time, "h")
metadata_abbreviated <- metadata$Taxonomy_metadata
metadata$Abbreviated <- metadata_abbreviated

#No duplicated sequences were included in this study
duplicated_sequences <- NULL

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Process CAZy domains
#All domains
terms_interest <- c("")
exclude_terms <- NULL
#CAZy domains
all_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                            latest_dbcan_metadata, terms_interest, exclude_terms,
                            metadata, strain_seq, dendo_class)

#Select numeric data
all_domain_data <- all_cazy$selected_domain_data_numeric
all_domain_data <- all_domain_data[,2:ncol(all_domain_data)]

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Process CAZy domains
#Pectin
terms_interest <- c("pecti", "rhamnogalacturonan", "galacturon")
exclude_terms <- c("amylopectin")
#CAZy domains
pectin_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                            latest_dbcan_metadata, terms_interest, exclude_terms,
                            metadata, strain_seq, dendo_class)

#Xylan/Arabinoxylan
terms_interest <- c("xylan", "xylooligosaccharide")
exclude_terms <- NULL
#CAZy domains
xylan_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                            latest_dbcan_metadata, terms_interest, exclude_terms,
                            metadata, strain_seq, dendo_class)

#Mannan
terms_interest <- c("mann")
exclude_terms <- c("mannosylglycerate", "alpha-D-mannoside", "beta-D-mannoside")
#CAZy domains
mannan_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                           latest_dbcan_metadata, terms_interest, exclude_terms,
                           metadata, strain_seq, dendo_class)

#Cellulose
terms_interest <- c("cellul", "cello")
exclude_terms <- NULL
#CAZy domains
cellulose_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                            latest_dbcan_metadata, terms_interest, exclude_terms,
                            metadata, strain_seq, dendo_class)

#Glucan
terms_interest <- c("glucan")
exclude_terms <- NULL
#CAZy domains
glucan_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                            latest_dbcan_metadata, terms_interest, exclude_terms,
                            metadata, strain_seq, dendo_class)

#Fructan
terms_interest <- c("fructan", "oligofructose", "fructooligosaccharide", "fructo-oligosaccharide", "inulin")
exclude_terms <- NULL
#CAZy domains
# fructan_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
#                             latest_dbcan_metadata, terms_interest, exclude_terms,
#                             metadata, strain_seq, dendo_class)

#Lignin
terms_interest <- c("lign")
exclude_terms <- NULL
#CAZy domains
lignin_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                               latest_dbcan_metadata, terms_interest, exclude_terms,
                               metadata, strain_seq, dendo_class)
#Pullulan and cyclodextrin
terms_interest <- c("pullul", "cyclodextrin")
exclude_terms <- NULL
#CAZy domains
pullulan_cyclodextrin_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                                           latest_dbcan_metadata, terms_interest, exclude_terms,
                                           metadata, strain_seq, dendo_class)

#Maltodextrin and malto-oligosaccharides
terms_interest <- c("malto-oligosaccharide", "maltooligosaccharide", "maltodextrin")
exclude_terms <- NULL
#CAZy domains
mos_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                         latest_dbcan_metadata, terms_interest, exclude_terms,
                         metadata, strain_seq, dendo_class)

#All fibre structures
terms_interest <- c("pecti", "rhamnogalacturonan", "galacturon",
                    "xylan", "xylooligosaccharide",
                    "mann",
                    "cellul", "cello",
                    "glucan",
                    "fructan", "oligofructose", "fructooligosaccharide", "fructo-oligosaccharide", "inulin",
                    "lign",
                    "pullul", "cyclodextrin",
                    "malto-oligosaccharide", "maltooligosaccharide", "maltodextrin")
exclude_terms <- c("amylopectin",
                   "mannosylglycerate", "alpha-D-mannoside", "beta-D-mannoside")
#CAZy domains
fibre_cazy <- process_cazy(file_list, coverage_threshold, remove_taxa_without_domains,
                           latest_dbcan_metadata, terms_interest, exclude_terms,
                           metadata, strain_seq, dendo_class)

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#PECTIN
sel_data <- pectin_cazy$selected_domain_data_numeric
sel_description <- pectin_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
#Remove domains not found in the genomes as well as unspecific domains 
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
if (length(domains_not_found) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
  sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
}
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
                                        c("#ff0000", "#00ff00", "#0000ff")) 
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
pectin_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                        strain_levels, sel_description, 20, 1)
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#XYLAN
sel_data <- xylan_cazy$selected_domain_data_numeric
sel_description <- xylan_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
#Remove domains not found in the genomes as well as unspecific domains 
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
if (length(domains_not_found) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
  sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
}
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
                                        c("#ff0000", "#00ff00", "#0000ff")) 
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
xylan_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                        strain_levels, sel_description, 20, 1)
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#MANNAN
sel_data <- mannan_cazy$selected_domain_data_numeric
sel_description <- mannan_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
#Remove domains not found in the genomes as well as unspecific domains 
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
# if (length(domains_not_found) != 0) {
#   sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
#   sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
# }
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
                                        c("#ff0000", "#00ff00", "#0000ff")) 
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
mannan_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                       strain_levels, sel_description, 20, 1)
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#CELLULOSE
sel_data <- cellulose_cazy$selected_domain_data_numeric
sel_description <- cellulose_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
# unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
domains_not_found <- NULL
#GH5 is very important for cellulose metabolism and is considered an unspecific domain
unspecific_domains <- NULL
#Remove domains not found in the genomes as well as unspecific domains 
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
if (length(domains_not_found) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
  sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
}
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
                                        c("#ff0000", "#00ff00", "#0000ff")) 
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
cellulose_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                       strain_levels, sel_description, 20, 1)
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#GLUCAN
sel_data <- glucan_cazy$selected_domain_data_numeric
sel_description <- glucan_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
#Remove domains not found in the genomes as well as unspecific domains 
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
if (length(domains_not_found) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
  sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
}
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
                                        c("#ff0000", "#00ff00", "#0000ff")) 
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
glucan_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                           strain_levels, sel_description, 20, 2)
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

# #FRUCTAN
# sel_data <- fructan_cazy$selected_domain_data_numeric
# sel_description <- fructan_cazy$selected_cazy_metadata_sorted
# #Get missing and unspecific domains
# domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
# unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
# #Remove domains not found in the genomes as well as unspecific domains 
# if (length(unspecific_domains) != 0) {
#   sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
#   sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
# }
# if (length(domains_not_found) != 0) {
#   sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
#   sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
# }
# sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))
# 
# #Load extra libraries
# library(paletteer)
# library(Polychrome)
# sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
# palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
# palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
#                                         c("#ff0000", "#00ff00", "#0000ff")) 
# #Set group levels
# strain_levels <- unique(sel_data$Class)
# #Generate heatmaps and dendograms
# x11()
# fructan_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
#                                                         strain_levels, sel_description, 20, 1)
# dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#LIGNIN
sel_data <- lignin_cazy$selected_domain_data_numeric
sel_description <- lignin_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
#Remove domains not found in the genomes as well as unspecific domains
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
if (length(domains_not_found) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
  sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
}
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),
                                        c("#ff0000", "#00ff00", "#0000ff"))
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
lignin_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                       strain_levels, sel_description, 20, 1)
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Pullulan and cyclodextrin
sel_data <- pullulan_cyclodextrin_cazy$selected_domain_data_numeric
sel_description <- pullulan_cyclodextrin_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
#Remove domains not found in the genomes as well as unspecific domains
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
if (length(domains_not_found) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
  sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
}
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
                                        c("#ff0000", "#00ff00", "#0000ff")) 
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
pullulan_cyclodextrin_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                                       strain_levels, sel_description, 20, 1)
dev.off()

#Export data for correlation analysis
saveRDS(sel_data, "R_objects/strains_pullulan.rds")

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Maltodextrin and malto-oligosaccharides
sel_data <- mos_cazy$selected_domain_data_numeric
sel_description <- mos_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
#Remove domains not found in the genomes as well as unspecific domains
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
if (length(domains_not_found) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
  sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
}
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
                                        c("#ff0000", "#00ff00", "#0000ff")) 
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
mos_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                     strain_levels, sel_description, 20, 1)
dev.off()

#Export data for correlation analysis
saveRDS(sel_data, "R_objects/strains_mos.rds")

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#All fibre structures
sel_data <- fibre_cazy$selected_domain_data_numeric
sel_description <- fibre_cazy$selected_cazy_metadata_sorted
#Get missing and unspecific domains
domains_not_found <- names(which(colSums(sel_data[,2:ncol(sel_data)]) == 0))
unspecific_domains <- rownames(sel_description)[grep("/", sel_description$Group)]
#Remove domains not found in the genomes as well as unspecific domains
if (length(unspecific_domains) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% unspecific_domains)]
  sel_description <- sel_description[!(row.names(sel_description) %in% unspecific_domains),]
}
if (length(domains_not_found) != 0) {
  sel_data <- sel_data[,-which(names(sel_data) %in% domains_not_found)]
  sel_description <- sel_description[!(row.names(sel_description) %in% domains_not_found),]
}
sel_description$Group <- factor(sel_description$Group, levels=unique(sel_description$Group))

#Load extra libraries
library(paletteer)
library(Polychrome)
sel_data$Class <- factor(sel_data$Class, levels=unique(sel_data$Class))
palette_to_use_dendogram <- createPalette(length(unique(metadata$Substrate_Time)), c("#ff0000", "#00ff00", "#0000ff")) 
palette_to_use_heatmap <- createPalette(length(levels(sel_description$Group)),  
                                        c("#ff0000", "#00ff00", "#0000ff")) 
#Set group levels
strain_levels <- unique(sel_data$Class)
#Generate heatmaps and dendograms
x11()
fibre_cazy_plots <- heatmap_dendogram_two_color_scales(sel_data, 20, 20, palette_to_use_dendogram, palette_to_use_heatmap,
                                                       strain_levels, sel_description, 20, 2)
dev.off()

#Export data for correlation analysis
saveRDS(sel_data, "R_objects/strains_fibre.rds")

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Print figures
#Pectin
x11()
grid.newpage()
pectin_cazy_plots$heatmap_plot + theme(plot.margin = unit(c(0,0,1,0), "cm"))
draw(pectin_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
pectin_cazy_heatmap <- recordPlot()
dev.off()

#Xylan
x11()
grid.newpage()
xylan_cazy_plots$heatmap_plot
draw(xylan_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
xylan_cazy_heatmap <- recordPlot()
dev.off()

#Mannan
x11()
grid.newpage()
mannan_cazy_plots$heatmap_plot
draw(mannan_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
mannan_cazy_heatmap <- recordPlot()
dev.off()

#Cellulose
x11()
grid.newpage() 
cellulose_cazy_plots$heatmap_plot + theme(plot.margin = unit(c(0,0,1,0), "cm"))
draw(cellulose_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
cellulose_cazy_heatmap <- recordPlot()
dev.off()

#Glucan
x11()
grid.newpage()
glucan_cazy_plots$heatmap_plot + theme(plot.margin = unit(c(0,0,1,0), "cm"))
draw(glucan_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
glucan_cazy_heatmap <- recordPlot()
dev.off()

# #Fructan
# x11()
# grid.newpage()
# fructan_cazy_plots$heatmap_plot
# draw(fructan_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
# fructan_cazy_heatmap <- recordPlot()
# dev.off()

#Lignin
x11()
grid.newpage()
lignin_cazy_plots$heatmap_plot
draw(lignin_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
lignin_cazy_heatmap <- recordPlot()
dev.off()

#Pullulan and cyclodextrin
x11()
grid.newpage()
pullulan_cyclodextrin_cazy_plots$heatmap_plot
draw(pullulan_cyclodextrin_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
pullulan_cyclodextrin_cazy_heatmap <- recordPlot()
dev.off()

#Maltodextrin and malto-oligosaccharides
x11()
grid.newpage()
mos_cazy_plots$heatmap_plot
draw(mos_cazy_plots$heatmap_plot_legend, x=unit(1.8, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
mos_cazy_heatmap <- recordPlot()
dev.off()

#All fibre structures
x11()
grid.newpage()
fibre_cazy_plots$heatmap_plot + theme(plot.margin = unit(c(0,0,2.6,0), "cm"))
draw(fibre_cazy_plots$heatmap_plot_legend, x=unit(1.4, "cm"), y=unit(0.1, "cm"), just=c("left", "bottom"))
fibre_cazy_heatmap <- recordPlot()
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Export high-quality Figures
#Pectin
win.metafile() #Must be full screen
png("Figures/CAZy/pectin_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
pectin_cazy_heatmap
dev.off()

#Xylan
win.metafile() #Must be full screen
png("Figures/CAZy/xylan_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
xylan_cazy_heatmap
dev.off()

#Mannan
win.metafile() #Must be full screen
png("Figures/CAZy/mannan_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
xylan_cazy_heatmap
dev.off()

#Cellulose
win.metafile() #Must be full screen
png("Figures/CAZy/cellulose_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
cellulose_cazy_heatmap
dev.off()

#Glucan
win.metafile() #Must be full screen
png("Figures/CAZy/glucan_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
glucan_cazy_heatmap
dev.off()

# #Fructan
# win.metafile() #Must be full screen
# png("Figures/CAZy/fructan_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
# fructan_cazy_heatmap
# dev.off()

#Lignin
win.metafile() #Must be full screen
png("Figures/CAZy/lignin_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
lignin_cazy_heatmap
dev.off()

#Pullulan and cyclodextrin
win.metafile() #Must be full screen
png("Figures/CAZy/pullulan_cyclodextrin_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
pullulan_cyclodextrin_cazy_heatmap
dev.off()

#Maltodextrin and malto-oligosaccharides
win.metafile() #Must be full screen
png("Figures/CAZy/mos_cazy_heatmap_bsg.png", width=28, height=14, units="in", res=600)
mos_cazy_heatmap
dev.off()

#All fibre structures
win.metafile() #Must be full screen
png("Figures/CAZy/fibre_cazy_heatmap_bsg.png", width=55, height=16, units="in", res=600)
fibre_cazy_heatmap
dev.off()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Finish
