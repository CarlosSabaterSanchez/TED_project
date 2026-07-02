
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Load libraries
library(agricolae)
library(FSA)
library(rcompanion)
library(car)
library(reshape2)
library(ggplot2)
library(ggfortify)
library(ggtext)
library(dplyr)

#Load custom functions
source("Functions/1_Statistical_differences.R")

#Figures and Tables showing statistical differences
#Figure 1c (Supplementary material Table S2)
#Figure 1d (Supplementary material Table S5)
#Supplementary material Table S2
table_s2 <- openxlsx::read.xlsx("Tables_S2_S5.xlsx", sheet=1, rowNames=T)
names(table_s2)[1] <- "Class"
table(table_s2$Class)
#Supplementary material Table S5
table_s5 <- openxlsx::read.xlsx("Tables_S2_S5.xlsx", sheet=2, rowNames=T)
names(table_s5)[1] <- "Class"
table(table_s5$Class)

#Reformat labels
#Table S2
table_s2$Class <- gsub(1, "BSG-U", table_s2$Class)
table_s2$Class <- gsub(2, "BSG-F-2d", table_s2$Class)
table_s2$Class <- gsub(3, "BSG-F-6d", table_s2$Class)
table_s2$Class <- gsub(4, "BSG-Vz-U", table_s2$Class)
table_s2$Class <- gsub(5, "BSG-Vz-F-2d", table_s2$Class)
table_s2$Class <- gsub(6, "BSG-Vz-F-6d", table_s2$Class)

#Table S5
table_s5$Class <- gsub(1, "BAB-U", table_s5$Class)
table_s5$Class <- gsub(2, "BAB-F-2d", table_s5$Class)
table_s5$Class <- gsub(3, "BAB-F-7d", table_s5$Class)
table_s5$Class <- gsub(4, "FAB-U", table_s5$Class)
table_s5$Class <- gsub(5, "FAB-F-2d", table_s5$Class)
table_s5$Class <- gsub(6, "FAB-F-7d", table_s5$Class)

#Replace NAs by 0
table_s2[is.na(table_s2)] <- 0
table_s5[is.na(table_s5)] <- 0

#Run statistical tests - Table S2
stats_table_s2 <- multiple_stats_no_pseq(estimator_list=list(Table_S2=table_s2), metadata_class_no_formula="Class")
#Inspect statistical differences
#Significant
stats_table_s2$Table_S2_significant_post_hoc$Parametric_post_hoc
stats_table_s2$Table_S2_significant_post_hoc$Parametric_effect_size
stats_table_s2$Table_S2_parametric_significant_p_values
stats_table_s2$Table_S2_significant_post_hoc$Non_parametric_post_hoc
stats_table_s2$Table_S2_significant_post_hoc$Non_parametric_effect_size
stats_table_s2$Table_S2_non_parametric_significant_p_values
#Non-significant
stats_table_s2$Table_S2_non_significant_post_hoc$Parametric_post_hoc
stats_table_s2$Table_S2_non_significant_post_hoc$Parametric_effect_size
stats_table_s2$Table_S2_parametric_non_significant_p_values
stats_table_s2$Table_S2_non_significant_post_hoc$Non_parametric_post_hoc
stats_table_s2$Table_S2_non_significant_post_hoc$Non_parametric_effect_size
stats_table_s2$Table_S2_non_parametric_non_significant_p_values

#Run statistical tests - Table S5
stats_table_s5 <- multiple_stats_no_pseq(list(Table_S5=table_s5), "Class")
#Inspect statistical differences
#Significant
stats_table_s5$Table_S5_significant_post_hoc$Parametric_post_hoc
stats_table_s5$Table_S5_significant_post_hoc$Parametric_effect_size
stats_table_s5$Table_S5_parametric_significant_p_values
stats_table_s5$Table_S5_significant_post_hoc$Non_parametric_post_hoc
stats_table_s5$Table_S5_significant_post_hoc$Non_parametric_effect_size
stats_table_s5$Table_S5_non_parametric_significant_p_values
#Non-significant
stats_table_s5$Table_S5_non_significant_post_hoc$Parametric_post_hoc
stats_table_s5$Table_S5_non_significant_post_hoc$Parametric_effect_size
stats_table_s5$Table_S5_parametric_non_significant_p_values
stats_table_s5$Table_S5_non_significant_post_hoc$Non_parametric_post_hoc
stats_table_s5$Table_S5_non_significant_post_hoc$Non_parametric_effect_size
stats_table_s5$Table_S5_non_parametric_non_significant_p_values

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Finish
