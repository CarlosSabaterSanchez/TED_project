
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to process CAZy domains
process_cazy <- function(file_list, coverage_threshold, remove_taxa_without_domains,
                         latest_dbcan_metadata, terms_interest, exclude_terms,
                         metadata, strain_seq, dendo_class) {
  
  #Compute heatmaps
  #Read individual csv files from the file list and store them in a raw file list
  raw_files <- list()
  for (i in 1:length(file_list)) {
    raw_files[[i]] <- read.csv(paste0("cazy_functional_domains/", file_list[i]), header=T, sep="\t")
  }
  #Export raw data as an R object
  saveRDS(raw_files, "R_objects/raw_files.rds")
  #Read R object containing the raw data
  raw_files <- readRDS("R_objects/raw_files.rds")
  #Iterate individual files to select high-quality functional domains
  filtered_domains <- list()
  for (i in 1:length(raw_files)) {
    filtered_domains[[i]] <- raw_files[[i]][which(raw_files[[i]]$Coverage >= coverage_threshold), ]
  }
  #Overwrite previous list
  raw_files <- filtered_domains
  #Set the names of the list containing raw data
  names(raw_files) <- sub("(.)", "\\U\\1", gsub("_hmmer.*", "", file_list), perl=T)
  #Gather codes indicating domain families from all samples present in the list containing filtered domain data
  domain_data <- data.frame()
  for (i in 1:length(raw_files)) {
    domain_data <- rbind(domain_data, as.data.frame(raw_files[[i]]$HMM.Profile))
  }
  #Remove duplicated domain codes from the domain data frame
  domain_data <- unique(domain_data)
  names(domain_data) <- c("Domain")
  #Export domain data as an R object
  saveRDS(domain_data, "R_objects/domain_data.rds")
  #Read R object containing domain data
  domain_data <- readRDS("R_objects/domain_data.rds")
  
  #Search for interesting functional domains using manually-curated pfam codes
  # pfam_codes <- openxlsx::read.xlsx("../Input/Selected_Pfam.xlsx", sheet=1)
  # pfam_selected_codes <- pfam_codes$Pfam
  #Prune domain data to select interesting functional domains
  # domain_data$Domain <- gsub("[.].*", "", domain_data$Domain)
  rownames(domain_data) <- domain_data$Domain
  # domain_data <- domain_data[intersect(pfam_selected_codes, rownames(domain_data)),, drop=F]
  
  #Generate a new data frame summarising the number of glycosidase families per genome sequence
  for (i in 1:length(raw_files)) {
    
    #Determine glycosidase families present in i genome sequence
    i_domains <- intersect(domain_data$Domain, raw_files[[i]]$HMM.Profile)
    # i_domains <- intersect(domain_data$Domain, gsub("[.].*", "", raw_files[[i]]$HMM.Profile))
    i_domains <- data.frame(i_domains)
    
    #Continue the analysis if any glycosidases are found in i genome sequence
    if (length(i_domains) > 0) {
      names(i_domains) <- names(raw_files)[i]
      
      #If the number of domains found in i genome is lower than total number of domains found in the complete dataset, 
      #Then, add NAs corresponding to the missing domains
      if (nrow(i_domains) < nrow(domain_data)) {
        i_add_nas <- nrow(domain_data) - nrow(i_domains)
        i_domains[nrow(i_domains)+i_add_nas,] <- NA
      }
      
      #Align domains found in i genome to match the total number of domains found in the complete dataset
      i_match_index <- match(domain_data$Domain, i_domains[,1])
      i_domains_match  <- i_domains[i_match_index,]
      
      #glycosidases found in i genomes to the complete dataset
      domain_data <- cbind(domain_data, i_domains_match)
      #Set the column name corresponding to i genome
      names(domain_data)[ncol(domain_data)] <- names(i_domains)
      
    }
  }
  
  #Format the the data frame summarising the number of domains per genome sequence
  #Set row names (i.e. functional domain) 
  rownames(domain_data) <- gsub("[.].*", "", domain_data$Domain)
  #Remove column containing domain names 
  domain_data <- domain_data[,2:ncol(domain_data)]
  #Replace characters in cells by 1
  domain_data_numeric <- domain_data
  domain_data_numeric <- domain_data_numeric %>% mutate_all(funs(str_replace_all(., ".*hmm", "1")))
  #Replace NAs by 0
  domain_data_numeric[is.na(domain_data_numeric)] <- 0
  #Convert the data frame to numeric
  domain_data_numeric <- mutate_all(domain_data_numeric, function(x) as.numeric(as.character(x)))
  
  #Select specific CAZy domains (i.e. pectinases)
  #Read CAZy database metadata summarising the enzymatic activities of CAZy functional domains and their substrates
  cazy_metadata <- read.csv(latest_dbcan_metadata, sep="\t", skip=1)
  cazy_metadata$Substrate_Simple_V2 <- iconv(cazy_metadata$Substrate_Simple_V2, from = "ISO-8859-1", to = "UTF-8")
  selected_cazy_metadata <- cazy_metadata[grep(paste0(terms_interest, collapse="|"), cazy_metadata$Substrate_Simple_V2, ignore.case=T),]
  #Exclude search terms
  if (length(exclude_terms) != 0) {
    selected_cazy_metadata <- selected_cazy_metadata[-grep(paste0(exclude_terms, collapse="|"), selected_cazy_metadata$Substrate_Simple_V2, ignore.case=T),]
  }
  #Remove duplicates
  selected_cazy_metadata <- selected_cazy_metadata[which(duplicated(selected_cazy_metadata$Family) == F),]
  #Order enzyme domains of interest
  selected_cazy_metadata <- selected_cazy_metadata[order(factor(selected_cazy_metadata$Family, 
                                                                levels=unique(mixedsort(selected_cazy_metadata$Family)))),]
  
  #Format selected activities (update this section depending on the substrate)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #Abbreviate activities
  #Pectin format####################################################################################################################################################
  selected_cazy_metadata$Group <- iconv(selected_cazy_metadata$Name, from = "ISO-8859-1", to = "UTF-8")
  selected_cazy_metadata$Group[grep("pectin methylesterase", selected_cazy_metadata$Group, ignore.case=T)] <- c("pectin methylesterase")
  selected_cazy_metadata$Group[grep("Galacturonidase", selected_cazy_metadata$Group, ignore.case=T)] <- c("Galacturonase")
  selected_cazy_metadata$Group[grep("rhamnogalacturonan alpha-1,2-galacturonohydrolase", selected_cazy_metadata$Group, ignore.case=T)] <- c("Rhamnogalacturonase")
  selected_cazy_metadata$Group[grep("rhamnogalacturonan alpha-L-rhamnohydrolase", selected_cazy_metadata$Group, ignore.case=T)] <- c("Rhamnogalacturonase")
  selected_cazy_metadata$Group[grep("aceric acid", selected_cazy_metadata$Group, ignore.case=T)] <- c("Aceric acid hydrolase")
  selected_cazy_metadata$Group[grep("exo-polygalacturonate lyase", selected_cazy_metadata$Group, ignore.case=T)] <- c("Exopolygalacturonate lyase")
  #Xylan format####################################################################################################################################################
  selected_cazy_metadata$Group[grep("copper-dependent lytic polysaccharide monooxygenases", selected_cazy_metadata$Group, ignore.case=T)] <- c("Polysaccharide monooxygenase")
  selected_cazy_metadata$Group[grep("lytic xylan monooxygenase", selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylan monooxygenase")
  selected_cazy_metadata$Group[grep("^CBM85$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucuronoxylan/glucan/glucomannan/cellulose active")
  selected_cazy_metadata$Group[grep("Binding to mannan, xylan, and cellulose", selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylan/mannan/cellulose active")
  selected_cazy_metadata$Group[grep("Binding to xylan, yeast cell wall glucan and chitin", selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylan/glucan/chitin active")
  selected_cazy_metadata$Group[c(intersect(grep("decorated soluble mannans", selected_cazy_metadata$Group, ignore.case=T),
                                           grep("xylan-degrading enzymes", selected_cazy_metadata$Group, ignore.case=T)))] <- c("Xylan/mannan/galactan active")
  selected_cazy_metadata$Group[grep("cellulose-binding function has been demonstrated in one case on amorphous cellulose and beta-1,4-xylan ", 
                                    selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylan/cellulose active")
  selected_cazy_metadata$Group[grep("The CBM72 found on an endoglucanase ", selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylan/glucan/mannan/cellulose active")
  selected_cazy_metadata$Group[grep("Several members of CBM37", selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylan/cellulose/chitin active")
  selected_cazy_metadata$Group[grep("The cellulose-binding function has been demonstrated in many cases. Several of these modules have been shown to also bind chitin or xylan", 
                                    selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylan/cellulose/chitin active")
  selected_cazy_metadata$Group[grep("The cellulose-binding function has been demonstrated in one case on amorphous cellulose and beta-1,4-xylan. Some of these modules also bind beta-1,3-glucan, beta-1,3-1,4-glucan, and beta-1,4-glucan", 
                                    selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylan/glucan/cellulose active")
  selected_cazy_metadata$Group[grep("beechwood xylan and rye arabinoxylan binding", selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylanase")
  selected_cazy_metadata$Group[grep("exo-oligoxylanase", selected_cazy_metadata$Group, ignore.case=T)] <- c("Exooligoxylanase")
  selected_cazy_metadata$Group[grep("Binding to arabinofuranose", selected_cazy_metadata$Group, ignore.case=T)] <- c("Arabinofuranosidase")
  selected_cazy_metadata$Group[grep("Leth et al", selected_cazy_metadata$Group, ignore.case=T)] <- c("Xylanase")
  selected_cazy_metadata$Group[grep("^CBM36$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylanase")
  selected_cazy_metadata$Group[grep("^CBM89$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylanase")
  selected_cazy_metadata$Group[grep("^CBM60$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylanase")
  selected_cazy_metadata$Group[grep("^CBM31$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylanase")
  selected_cazy_metadata$Group[grep("^CBM32$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylan/glucan active")
  selected_cazy_metadata$Group[grep("^CBM15$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylanase")
  selected_cazy_metadata$Group[grep("^CBM13$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylan/lectin active")
  selected_cazy_metadata$Group[grep("^CBM4$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylan/glucan/cellulose active")
  selected_cazy_metadata$Group[grep("^CBM22$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylan/glucan active")
  #Mannan format####################################################################################################################################################
  selected_cazy_metadata$Group[grep("^CBM88$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Galactoxyloglucan/galactomannan active")
  selected_cazy_metadata$Group[grep("^CBM76$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/Xyloglucan/glucomannan active")
  selected_cazy_metadata$Group[grep("^CBM80$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/Xyloglucan/glucomannan/galactomannan active")
  selected_cazy_metadata$Group[grep("^CBM62$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xyloglucan/galactomannan/arabinogalactan active")
  selected_cazy_metadata$Group[grep("^CBM29$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Mannanase")
  selected_cazy_metadata$Group[grep("^CBM23$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Mannanase")
  selected_cazy_metadata$Group[grep("^CBM27$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Mannanase")
  selected_cazy_metadata$Group[grep("^CBM16$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucomannan/cellulose active")
  #Lignin format####################################################################################################################################################
  selected_cazy_metadata$Group[grep("vanillyl-alcohol oxidase", selected_cazy_metadata$Group, ignore.case=T)] <- c("Vanillyl alcohol oxidase")
  selected_cazy_metadata$Group[grep("cytochrome-c peroxidase", selected_cazy_metadata$Group, ignore.case=T)] <- c("Cytochrome C peroxidase")
  #Cellulose format####################################################################################################################################################
  selected_cazy_metadata$Group[grep("^CBM1$", selected_cazy_metadata$Family, ignore.case=T)] <- c("cellulose/chitin active")
  selected_cazy_metadata$Group[grep("^CBM3$", selected_cazy_metadata$Family, ignore.case=T)] <- c("cellulose/chitin active")
  selected_cazy_metadata$Group[grep("^CBM8$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM9$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM10$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM11$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/cellulose active")
  selected_cazy_metadata$Group[grep("^CBM17$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM28$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/cellulose active")
  selected_cazy_metadata$Group[grep("^CBM30$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM44$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xyloglucan/cellulose active")
  selected_cazy_metadata$Group[grep("^CBM46$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM49$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM63$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM64$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellulase")
  selected_cazy_metadata$Group[grep("^CBM81$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/xyloglucan/cello-oligosaccharide active")
  selected_cazy_metadata$Group[grep("^GH5$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/cellulose active")
  selected_cazy_metadata$Group[grep("^GH7$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellobiohydrolase")
  selected_cazy_metadata$Group[grep("^GH48$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellobiohydrolase")
  selected_cazy_metadata$Group[grep("^GH74$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Cellobiohydrolase")
  selected_cazy_metadata$Group[grep("^GH131$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/cellulose active")
  #Glucan format####################################################################################################################################################
  selected_cazy_metadata$Group[grep("^CBM39$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/lipopolysaccharide/lipoteichoic acid active")
  selected_cazy_metadata$Group[grep("^CBM41$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan/amylose/amylopectin/pullulan active")
  selected_cazy_metadata$Group[grep("^CBM43$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^CBM52$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^CBM54$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylan/glucan/chitin active")
  selected_cazy_metadata$Group[grep("^CBM56$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^CBM65$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xyloglucanase")
  selected_cazy_metadata$Group[grep("^CBM75$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xyloglucanase")
  selected_cazy_metadata$Group[grep("^CBM78$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^CBM79$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH1$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glycosidase")
  selected_cazy_metadata$Group[grep("^GH3$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH9$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xyloglucanase")
  selected_cazy_metadata$Group[grep("^GH10$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH16$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Transglycosidase")
  selected_cazy_metadata$Group[grep("^GH17$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanosyltransglycosylase")
  selected_cazy_metadata$Group[grep("^GH26$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH30$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucosidase")
  selected_cazy_metadata$Group[grep("^GH31$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xylosidase")
  selected_cazy_metadata$Group[grep("^GH45$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Xyloglucanase")
  selected_cazy_metadata$Group[grep("^GH51", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH64$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH65$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucosidase/phosphorylase")
  selected_cazy_metadata$Group[grep("^GH71$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH72$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanosyltransglycosylase")
  selected_cazy_metadata$Group[grep("^GH81$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH87$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH94$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Oligoglucan phosphorylase")
  selected_cazy_metadata$Group[grep("^GH128$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucosidase")
  selected_cazy_metadata$Group[grep("^GH132$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucosidase")
  selected_cazy_metadata$Group[grep("^GH144$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucooligosaccharide sophorohydrolase")
  selected_cazy_metadata$Group[grep("^GH148$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH149$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan phosphorylase")
  selected_cazy_metadata$Group[grep("^GH152$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH157$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH158$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  selected_cazy_metadata$Group[grep("^GH161$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucan phosphorylase")
  selected_cazy_metadata$Group[grep("^GH162$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Glucanase")
  #Fructan format
  selected_cazy_metadata$Group[grep("^CBM38$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Inulin active")
  selected_cazy_metadata$Group[grep("^CBM66$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Fructan active")
  selected_cazy_metadata$Group[grep("^GH91$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Difructofuranose hydrolase")
  selected_cazy_metadata$Group[grep("^GH172$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Difructose synthase")
  #Pullulan and cyclodextrin format
  selected_cazy_metadata$Group[grep("^CBM20$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Starch/cyclodextrin active")
  #Maltodextrin and malto-oligosaccharides format
  selected_cazy_metadata$Group[grep("^GH4$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Alpha glucosidase")
  selected_cazy_metadata$Group[grep("^GH13$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Oligosaccharide alpha glucosyltransferase")
  selected_cazy_metadata$Group[grep("^GH27$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Isomalto dextranase")
  selected_cazy_metadata$Group[grep("^GH57$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Alpha glucanotransferase")
  selected_cazy_metadata$Group[grep("^GH63$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Alpha glucosidase")
  selected_cazy_metadata$Group[grep("^GH126$", selected_cazy_metadata$Family, ignore.case=T)] <- c("Alpha amylase")
  #General format####################################################################################################################################################
  #Remove bonds and parentheses
  selected_cazy_metadata$Group <- gsub(" /.*", "", selected_cazy_metadata$Group)
  selected_cazy_metadata$Group <- gsub(" \\(.*", "", selected_cazy_metadata$Group)
  selected_cazy_metadata$Group <- gsub(".*-", "", selected_cazy_metadata$Group)
  selected_cazy_metadata$Group <- str_to_sentence(selected_cazy_metadata$Group)
  
  #Find common enzyme domains between CAZy database metadata and genome annotations (be careful with subfamilies "_", i.e. GH13_* or GH43_*)
  common_domains <- grep(paste0("^", selected_cazy_metadata$Family, "$", collapse="|"), rownames(domain_data_numeric))
  common_domains_with_subfamilies <- grep(paste0("^", selected_cazy_metadata$Family, "_", collapse="|"), rownames(domain_data_numeric))
  #Prune genome annotations to select only those specific CAZy domains (i.e. pectinases)
  selected_domain_data_numeric <- domain_data_numeric[unique(c(common_domains, common_domains_with_subfamilies)),]
  #Add subfamilies to CAZy metadata
  common_domains_with_subfamilies_codes <- rownames(domain_data_numeric)[common_domains_with_subfamilies]
  #Iterate common domains with subfamilies
  for (i in common_domains_with_subfamilies_codes) {
    i_no_subfamily <- gsub("_.*", "", i)
    i_cazy_metadata <- subset(selected_cazy_metadata, selected_cazy_metadata$Family == i_no_subfamily)
    i_cazy_metadata$Family <- i
    selected_cazy_metadata <- rbind(selected_cazy_metadata, i_cazy_metadata)
  }
  #Remove non-matching domains without subfamilies
  selected_cazy_metadata <- selected_cazy_metadata[match(rownames(selected_domain_data_numeric), selected_cazy_metadata$Family),]
  all(selected_cazy_metadata$Family == rownames(selected_domain_data_numeric))
  
  #Filter genomes not showing any selected CAZy domain
  if (remove_taxa_without_domains == T) {
    if (length(which(colSums(selected_domain_data_numeric) == 0)) != 0) {
      selected_domain_data_numeric <- selected_domain_data_numeric[,-which(colSums(selected_domain_data_numeric) == 0)]
    }
  }
  
  #Order selected enzymes according to their activities 
  selected_cazy_metadata <- selected_cazy_metadata[gtools::mixedorder(selected_cazy_metadata$Group),]
  selected_cazy_metadata_sorted <- data.frame()
  #Order selected enzymes according to their enzyme families 
  for (i in unique(selected_cazy_metadata$Group)) {
    i_subset <- subset(selected_cazy_metadata, selected_cazy_metadata$Group == i)
    i_subset <- i_subset[gtools::mixedorder(i_subset$Family),]
    selected_cazy_metadata_sorted <- rbind(selected_cazy_metadata_sorted, i_subset)
  }
  
  #Order annotation
  selected_cazy_metadata_sorted$Family <- gsub("_", " ", selected_cazy_metadata_sorted$Family)
  rownames(selected_cazy_metadata_sorted) <- selected_cazy_metadata_sorted$Family
  rownames(selected_domain_data_numeric) <- gsub("_", " ", rownames(selected_domain_data_numeric))
  #Format activities
  selected_domain_data_numeric <- selected_domain_data_numeric[selected_cazy_metadata_sorted$Family,]
  
  #Order annotation
  selected_domain_data_numeric <- as.data.frame(t(selected_domain_data_numeric))
  
  #Select sequence files showing CAZy domains without duplicates
  sequences_annotated_without_duplicates <- intersect(gsub("_cazy", "", rownames(selected_domain_data_numeric)), strain_seq)
  all_removed_strains_are_dupli <- length(grep(paste0("^", setdiff(gsub("_cazy", "", rownames(selected_domain_data_numeric)), strain_seq), "$", collapse="|"),
                                               duplicated_sequences)) == length(setdiff(gsub("_cazy", "", rownames(selected_domain_data_numeric)), strain_seq))
  print(paste0("All srains removed were duplicates: ", all_removed_strains_are_dupli))
  metadata <- metadata[grep(paste0(sequences_annotated_without_duplicates, collapse="|"), metadata$ID),]
  #Get annotated and non-duplicated strain sequences again
  rownames(selected_domain_data_numeric) <- gsub("Cazy_out", "cazy_out", rownames(selected_domain_data_numeric))
  length(intersect(paste0(strain_seq), rownames(selected_domain_data_numeric))) == length(strain_seq)
  selected_domain_data_numeric <- selected_domain_data_numeric[paste0(strain_seq),]
  rownames(selected_domain_data_numeric) <- gsub("cazy_out_", "", rownames(selected_domain_data_numeric))
  #Add strain codes
  rownames(selected_domain_data_numeric) <- paste0(gsub(" = .*", "", metadata$Taxonomy_metadata))
  # rownames(selected_domain_data_numeric) <- paste0(gsub(" = .*", "", metadata$Taxonomy_metadata), " ", gsub("_", " ", rownames(selected_domain_data_numeric)))
  selected_domain_data_numeric$Class <- metadata[[dendo_class]]
  
  #Remove alternative strain names
  rownames(selected_domain_data_numeric) <- gsub(".* / ", "", rownames(selected_domain_data_numeric))
  
  #Add strain abbreviations
  substr(rownames(selected_domain_data_numeric), 4, 4) <- "."
  rownames(selected_domain_data_numeric) <- sub("[.].*? ", ". ", rownames(selected_domain_data_numeric))
  
  #Set factor levels
  selected_cazy_metadata_sorted$Group <- factor(selected_cazy_metadata_sorted$Group, levels=unique(selected_cazy_metadata_sorted$Group))
  
  #Get dendogram species and strains
  #Format dendogram to mix both italics and normal fonts in each label
  #https://stackoverflow.com/questions/63788793/how-to-apply-both-italic-and-normal-fonts-in-the-same-label-in-phylogenetic-tree
  #https://r-graph-gallery.com/package/ggtext.html
  #https://stackoverflow.com/questions/72027474/how-to-format-selected-axis-text-words-with-italic-and-bold-using-ggtext
  selected_df <- selected_domain_data_numeric
  species <- ifelse(is.na(word(rownames(selected_df), 1, 2)), rownames(selected_df), word(rownames(selected_df), 1, 2))
  n_words <- str_count(rownames(selected_df), " ") + 1
  strains <- ifelse(is.na(word(rownames(selected_df), 3, n_words)), rownames(selected_df), word(rownames(selected_df), 3, n_words))
  format_species_strains <- paste0("paste(italic('", species, "'),", " ' ", strains, "')")
  format_species_strains <- paste0("*", species, "* ", "**", strains, "**")
  rownames(selected_domain_data_numeric) <- format_species_strains
  
  #Set strain class as the first column
  selected_domain_data_numeric <- selected_domain_data_numeric[,c(length(selected_domain_data_numeric), 1:(length(selected_domain_data_numeric)-1))]
  
  #Append results
  results <- list(selected_domain_data_numeric, selected_cazy_metadata_sorted)
  names(results) <- c("selected_domain_data_numeric", "selected_cazy_metadata_sorted")
  #Return results
  return(results)
  
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

#Custom function to process CAZy domains by taxa
process_cazy_by_taxa <- function(selected_domain_data_numeric) {
  
  #Get the number of different taxa at this level
  taxa <- unique(selected_domain_data_numeric$Class)
  #Create two new data frames to store domain percentages by taxa
  domain_data_by_taxa <- data.frame(matrix(nrow=0, ncol=0))
  domain_description_by_taxa <- data.frame(matrix(nrow=0, ncol=0))
  
  #Generate a new data frame summarizing the number of functional domains per taxa
  for (i in 1:length(taxa)) {
    #Select genomes corresponding to i taxa
    i_taxa <- taxa[i]
    #Select domains corresponding to genomes from i taxa
    i_domains <- subset(selected_domain_data_numeric, selected_domain_data_numeric$Class == i_taxa,drop=F)
    #Get i taxa domain percentages
    i_percentages <- colSums(i_domains[,2:ncol(i_domains)])*100/nrow(i_domains)
    i_n_genomes <- colSums(i_domains[,2:ncol(i_domains)])
    i_percentages_description <- paste0(taxa[i], " - ", gsub("_", " ", names(i_domains)[2:length(names(i_domains))]), 
                                        " ", format(round(i_percentages, 2), nsmall =2), 
                                        "% (n=", i_n_genomes, "/", 
                                        nrow(i_domains), " genomes)")
    #Format results
    i_percentages <- t(data.frame(i_percentages))
    rownames(i_percentages) <- paste0("*", gsub("_", " ", taxa[i]), "*", " (n=", nrow(i_domains), ")")
    i_percentages_description <- data.frame(i_percentages_description)
    
    #Append results
    domain_data_by_taxa <- rbind(domain_data_by_taxa, i_percentages)
    domain_description_by_taxa <- rbind(domain_description_by_taxa, i_percentages_description)
  }
  
  #Append results
  results <- list(domain_data_by_taxa, domain_description_by_taxa)
  names(results) <- c("domain_data_by_taxa", "domain_description_by_taxa")
  #Return results
  return(results)
}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


#Finish
