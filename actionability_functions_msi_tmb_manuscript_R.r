#!/usr/bin/env Rscript

### Annotate IMPACT files using oncokb-annotator ###

### Chakravarty D, Gao J, Phillips SM, et al. OncoKB: A Precision Oncology Knowledge Base. JCO Precis Oncol. 2017;2017:PO.17.00011. doi:10.1200/PO.17.00011 ###

### Actionability Functions ###
# Collection of functions use to clean, process, and analysis actionability data

### Input parameters

# cna_df: OncoKB annotated IMPACT CNA data  
# mut_df: OncoKB annotated IMPACT mutation data  
# fus_df: OncoKB annotated IMPACT fusion data  
# clin_df: OncoKB annotated IMPACT clinical sample data  
# data_freeze: Sample data, must include *SAMPLE_ID*, group_col, and consent_col
# group_col: Column name for the groups (cancer types)
# consent_col: Columns name for 12-245 Part C consent status (YES/NO/NA)
# path_df: Pathway data, must include gene and correpsonding pathway columns (in that order)  
# tsg_list: List of tumor suppresor genes (no header)  
# fusion_list: List of genes to isolate from fusion partners (ie. NTRK1-LMNA fusion becomes NTRK1 fusion)  
# prop_level_df: Output from action_levels_barplot_fun actionability_levels_barplot_table.txt  
# alt_final_df: Output from action_alterations_barplot_fun actionability_master_alterations_table.txt  
# alt_min: Minimum alteration percentage required in one cancer type to visualize alteration on main plot (default 1)  
# status: Include only somatic mutations, only germline mutations, or both (options: somatic, germline, both)
# gene_order: List of genes for gene order, genes not included will be ordered by pathway following this list (no header) 
# only_highest_level: TRUE/FALSE, If true only visualize the highest level of evidence genes in main plot
# msi_tmb_status: TRUE/FALSE, If true include Level 1 MSI/TMB status in actionability barplot, removes MSI/TMB in all other plots
# msi_tmb_df: MSI/TMB annotated file (atypical alterations), visualizes MSI/TMB level 1 for actionability barplot, 
  # removes all samples in file for all other plots


###


# Load libraries
if (!require('tidyverse')) install.packages('tidyverse'); library(tidyverse)
if (!require('cowplot')) install.packages('cowplot'); library(cowplot)
if (!require('reshape2')) install.packages('reshape2'); library(reshape2)

# Collapse oncogenic alterations
collapse_oncogenic <- function(data_frame, sample_column, alteration_type){
  data_frame[, ] <- lapply(data_frame[, ], as.character)
  #data_frame_samp <- data_frame %>% dplyr::filter(oncogenic == "Oncogenic")  #### TESTING
  data_frame_samp <- data_frame[grepl("Oncogenic", data_frame$oncogenic),]
  colnames(data_frame_samp)[which(names(data_frame_samp) == sample_column)] <- "SAMPLE_ID"
  data_frame_samp <- aggregate(oncogenic ~ SAMPLE_ID, data = data_frame_samp, toString, na.omit = TRUE)
  colnames(data_frame_samp)[2] <- paste0(alteration_type, "_oncogenic")
  return(data_frame_samp)
}

# Create frequency data frame by group for subgroup
freq_dataframe <- function(data_frame, split_group, percentage_group){
  # Split group is the column to group by
  # Percentage group is the column to calculate the percentage for, by group
  df <- data_frame %>%
    dplyr::select(percentage_group, split_group) %>%
    group_by_(split_group, percentage_group, .drop = F) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = n / sum(n)) %>%
    ungroup()
  return(df)
}

# Create actionability level barplot
action_levels_barplot_fun <- function(cna_df, mut_df, fus_df, clin_df, data_freeze,
                                      status = c("somatic", "germline", "both"),
                                      group_col,
                                      consent_col,
                                      msi_tmb_status,
                                      msi_tmb_df){
  # Read in data
  cna_df <- read.delim(cna_df)
  fus_df <- read.delim(fus_df)
  mut_df <- read.delim(mut_df)
  clin_df <- read.delim(clin_df)
  data_freeze <- read.delim(data_freeze)
  data_freeze$SAMPLE_ID <- as.character(data_freeze$SAMPLE_ID)
  
  ######
  
  # Optional MSI/TMB addition
  if (msi_tmb_status == TRUE){
    msi_tmb_df <- read.delim(msi_tmb_path)
    msi_tmb_df <- msi_tmb_df %>%
      dplyr::select(SAMPLE_ID) %>% 
      mutate_if(is.factor, as.character) %>%
      mutate(Highest_level = "LEVEL_1_MSI-H_TMB-H") %>%
      distinct()
  } else {
    msi_tmb_df <- data.frame(SAMPLE_ID = character(), Highest_level = character())
  }
  
  #####
  
  # Clean & filter clinical data
  # Add group column
  clin_df <- clin_df %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID),
           HIGHEST_LEVEL = as.character(HIGHEST_LEVEL)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID) %>%
    left_join(data_freeze[,c("SAMPLE_ID", group_col, consent_col)], by = c("SAMPLE_ID")) %>%
    mutate(HIGHEST_LEVEL = ifelse(SAMPLE_ID %in% msi_tmb_df$SAMPLE_ID, "LEVEL_1_MSI-H_TMB-H", HIGHEST_LEVEL))
  group_col_dup <- paste0(group_col, ".y")
  colnames(clin_df)[which(names(clin_df) == group_col_dup)] <- group_col
  colnames(clin_df)[which(names(clin_df) == consent_col)] <- "consent"
  
  # Clean, filter, rename genomic data 
  # Fix column names if upper
  # Filter for columns of interest
  col_list <- c("SAMPLE_ID", "oncogenic", "LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B", "LEVEL_4", "Highest_level")
  
  # Fusions
  fus_df <- fus_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID),
           Highest_level = as.character(Highest_level)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID) %>%
    mutate(Highest_level = ifelse(SAMPLE_ID %in% msi_tmb_df$SAMPLE_ID, "LEVEL_1_MSI-H_TMB-H", Highest_level)) %>%
    dplyr::select(col_list)
  
  # CNA
  cna_df <- cna_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID),
           Highest_level = as.character(Highest_level)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID)  %>%
    mutate(Highest_level = ifelse(SAMPLE_ID %in% msi_tmb_df$SAMPLE_ID, "LEVEL_1_MSI-H_TMB-H", Highest_level)) %>%
    dplyr::select(col_list)
  
  # Mutations
  mut_df <- mut_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID),
           Highest_level = as.character(Highest_level)) %>%
    mutate(Highest_level = ifelse(SAMPLE_ID %in% msi_tmb_df$SAMPLE_ID, "LEVEL_1_MSI-H_TMB-H", Highest_level)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID)
  
  # Filter for status
  if (status == "somatic") {
    mut_somatic_df <- filter(mut_df, Mutation_Status != "GERMLINE" | is.na(Mutation_Status) == T)
    mut_germ_df <- filter(mut_df, Mutation_Status == "GERMLINE")
    mut_df <- mut_somatic_df[col_list]
  } else if (status == "germline") {
    clin_df <- filter(clin_df, consent == "YES")
    mut_germ_df <- filter(mut_df, Mutation_Status == "GERMLINE")
    mut_germ_df <- mut_germ_df[mut_germ_df$SAMPLE_ID %in% clin_df$SAMPLE_ID,]
    mut_df <- mut_germ_df[col_list]
  } else {
    mut_df <- mut_df[col_list]
    clin_germ_df <- filter(clin_df, consent == "YES" | is.na(consent) == T)
  }
  
  # Create master levels data frame for somatic 
  master_df <- rbind(cna_df, fus_df)
  master_df <- rbind(master_df, mut_df)
  master_df <- master_df %>%
    dplyr::select(SAMPLE_ID, Highest_level) %>%
    filter(Highest_level != "") %>%
    mutate_if(is.factor, as.character) %>%
    group_by(SAMPLE_ID) %>%
    dplyr::arrange(Highest_level) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  # Collapse oncogenic alterations
  cna_df <- collapse_oncogenic(cna_df, "SAMPLE_ID", "cna")
  fus_df <- collapse_oncogenic(fus_df, "SAMPLE_ID", "fus")
  mut_df <- collapse_oncogenic(mut_df, "SAMPLE_ID", "mut")
  
  # Filter if germline
  if (status == "germline") {
    clin_df <- left_join(clin_df, mut_germ_df[,c("SAMPLE_ID", "Highest_level")])
    clin_df <- clin_df %>% mutate_if(is.factor, as.character)
    clin_df$HIGHEST_LEVEL <- ifelse(clin_df$SAMPLE_ID %in% mut_germ_df$SAMPLE_ID, clin_df$Highest_level, "NO_LEVEL")
    # Get list of sample with oncogenic alteration
    onco_samp_list <- mut_df$SAMPLE_ID
  } else if (status == "somatic") {
    clin_df <- left_join(clin_df, master_df, by = "SAMPLE_ID")
    clin_df <- clin_df %>% mutate_if(is.factor, as.character)
    clin_df$HIGHEST_LEVEL <- ifelse(clin_df$SAMPLE_ID %in% mut_germ_df$SAMPLE_ID,
                                    clin_df$Highest_level, clin_df$HIGHEST_LEVEL)
    # Merge to make master oncogenic list of samples
    all_df <- full_join(cna_df, fus_df, by = "SAMPLE_ID")
    all_df <- full_join(all_df, mut_df, by = "SAMPLE_ID")
    # Get list of sample with oncogenic alteration
    onco_samp_list <- all_df$SAMPLE_ID
  } else {
    # Merge to make master oncogenic list of samples
    all_df <- full_join(cna_df, fus_df, by = "SAMPLE_ID")
    all_df <- full_join(all_df, mut_df, by = "SAMPLE_ID")
    # Get list of sample with oncogenic alteration
    onco_samp_list <- all_df$SAMPLE_ID
  }
  
  # Fill in the highest level blanks:
  clin_df$HIGHEST_LEVEL <- as.character(clin_df$HIGHEST_LEVEL)
  clin_df$HIGHEST_LEVEL[clin_df$HIGHEST_LEVEL == "" | is.na(clin_df$HIGHEST_LEVEL) == T] <- "NO_LEVEL"
  clin_df$HIGHEST_LEVEL[(clin_df$SAMPLE_ID %in% onco_samp_list) & (clin_df$HIGHEST_LEVEL == "NO_LEVEL") ] <- "ONCOGENIC"
  
  # For highest level of actionability, calculate the percentage of each level by subtype
  prop_level_df <- freq_dataframe(clin_df, group_col, "HIGHEST_LEVEL")
  
  # Set level order
  level_order <- c("LEVEL_1_MSI-H_TMB-H","LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B", "LEVEL_4", "ONCOGENIC", "NO_LEVEL")
  prop_level_df$HIGHEST_LEVEL <- factor(prop_level_df$HIGHEST_LEVEL, levels = level_order)
  
  # Add counts for labels
  # Check the number of oncotree codes and their frequency
  data_freeze <- data_freeze[data_freeze$SAMPLE_ID %in% clin_df$SAMPLE_ID,]
  clin_oncotree_freq <- as.data.frame(table(data_freeze[,group_col]))
  clin_oncotree_freq <- clin_oncotree_freq[order(clin_oncotree_freq$Freq, decreasing = T),]
  colnames(clin_oncotree_freq)[1] <- group_col
  prop_level_df <- left_join(prop_level_df, clin_oncotree_freq, by = group_col)
  
  if (status == "both") {
    data_freeze_2 <- data_freeze[data_freeze$SAMPLE_ID %in% clin_germ_df$SAMPLE_ID,]
    clin_oncotree_freq_germ <- as.data.frame(table(data_freeze_2[,group_col]))
    colnames(clin_oncotree_freq_germ)[1] <- group_col
    prop_level_df <- left_join(prop_level_df, clin_oncotree_freq_germ, by = group_col)
    prop_level_df$label <- apply(prop_level_df[ ,c(group_col, "Freq.x")], 1, paste0, collapse = " n=" )
    prop_level_df$label <- apply(prop_level_df[ ,c("label", "Freq.y")], 1, paste0, collapse = ":" )
  } else {
    prop_level_df$label <- apply(prop_level_df[ ,c(group_col, "Freq")], 1, paste0, collapse = " n=" )
  }
  
  # # Arrange by frequency of actionability
  # prop_level_df <- prop_level_df %>%
  #   arrange(HIGHEST_LEVEL, desc(freq))
  
  # Arrange by frequency of combined top 4 levels of actionability
  prop_level_df_order <- prop_level_df %>%
    filter(HIGHEST_LEVEL %in% c("LEVEL_1_MSI-H_TMB-H", "LEVEL_1", "LEVEL_2", "LEVEL_3A")) %>%
    group_by(CANCER_TYPE) %>%
    dplyr::mutate(sum_freq = sum(freq)) %>%
    right_join(prop_level_df) %>%
    dplyr::arrange(desc(sum_freq), HIGHEST_LEVEL, desc(freq)) %>%
    dplyr::rename(total_count = Freq) %>%
    mutate(CANCER_TYPE = factor(CANCER_TYPE, levels = unique(CANCER_TYPE)))
  
  # Save
  write.table(prop_level_df_order, "./actionability_levels_barplot_table.txt", sep = "\t", row.names = F, quote = F)
  
  # Set orders
  cancer_order <- unique(prop_level_df_order$label)
  
  # Plot breakdown of levels of evidence as a percentage by sarcoma subtype
  percent_bar_plot <- ggplot(prop_level_df_order, aes(y = freq, x = label, fill = HIGHEST_LEVEL)) +
    geom_col(position = position_stack(reverse = TRUE)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 6),
          axis.text.y = element_text(size = 6),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm"),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0.1, 0.05), "cm"),
          legend.justification="left",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-5)) +
    scale_y_continuous(limits=c(0, 1.00), expand = c(0, 0)) +
    scale_fill_manual(values = c("#88E281","#33A02C", "#1F78B4", "#984EA3", "#BE98CE", "#a8a8a8", "#ffdab9", "gray90"),
                      limits = c("LEVEL_1_MSI-H_TMB-H","LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B", "LEVEL_4", "ONCOGENIC", "NO_LEVEL"),
                      labels = c("LEVEL 1 MSI/TMB-H","LEVEL 1", "LEVEL 2", "LEVEL 3A", "LEVEL 3B", "LEVEL 4", "ONCOGENIC", "NO LEVEL")) +
    scale_x_discrete(position = "top",
                     limits = cancer_order) +
    ylab("Frequency") +
    labs(fill = "Highest Level of Evidence")
  
  return(percent_bar_plot)
  
}

# Create actionability alteration barplot
action_alterations_barplot_fun <- function(cna_df, mut_df, fus_df, clin_df, data_freeze,
                                           status = c("somatic", "germline", "both"),
                                           group_col, consent_col,
                                           prop_level_df = "./actionability_levels_barplot_table.txt",
                                           only_highest_level = F,
                                           msi_tmb_status,
                                           msi_tmb_df){
  # Read in data
  cna_df <- read.delim(cna_df)
  fus_df <- read.delim(fus_df)
  mut_df <- read.delim(mut_df)
  clin_df <- read.delim(clin_df)
  data_freeze <- read.delim(data_freeze)
  prop_level_df <- read.delim(prop_level_df)
  
  # Set order
  cancer_order_other <- as.character(unique(prop_level_df[,c(group_col)]))
  
  
  ######
  
  # Optional MSI/TMB addition
  if (msi_tmb_status == TRUE){
    msi_tmb_df <- read.delim(msi_tmb_path)
    msi_tmb_df <- msi_tmb_df %>%
      dplyr::select(SAMPLE_ID) %>% 
      mutate_if(is.factor, as.character) %>%
      mutate(Highest_level = "LEVEL_1_MSI-H_TMB-H") %>%
      distinct()
    data_freeze <- filter(data_freeze, !SAMPLE_ID %in% msi_tmb_df$SAMPLE_ID)
  } 
  
  #####
  
  # Clean & filter clinical data
  # Add group column
  clin_df <- clin_df %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID) %>%
    left_join(data_freeze[,c("SAMPLE_ID", group_col, consent_col)], by = c("SAMPLE_ID"))
  group_col_dup <- paste0(group_col, ".y")
  colnames(clin_df)[which(names(clin_df) == group_col_dup)] <- group_col
  colnames(clin_df)[which(names(clin_df) == consent_col)] <- "consent"
  
  # Clean, filter, rename genomic data 
  # Fix column names if upper
  fus_df <- fus_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID) %>%
    mutate(Fusion = gsub(" fusion", "", Fusion)) %>%
    mutate(Fusion = gsub(" - Archer", "", Fusion)) %>%
    dplyr::select(SAMPLE_ID, oncogenic, Highest_level, Fusion) %>%
    rowwise() %>%
    mutate(Fusion = ifelse(grepl("intragenic", Fusion), Fusion, 
                           paste(sort(unlist(strsplit(Fusion, "-", fixed = TRUE))), collapse = "-"))) %>%
    ungroup() %>%
    distinct() %>%
    dplyr::select(-Fusion) %>%
    mutate(ALTERATION = "Fusion") %>%
    filter(grepl("Oncogenic", oncogenic) == T, is.na(Highest_level) == F & Highest_level != "") %>%
    dplyr::select(SAMPLE_ID, ALTERATION, oncogenic, Highest_level)
  
  # CNA
  cna_df <- cna_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID)  %>%
    dplyr::select(SAMPLE_ID, ALTERATION, oncogenic, Highest_level) %>%
    filter(grepl("Oncogenic", oncogenic) == T, is.na(Highest_level) == F & Highest_level != "") %>%
    dplyr::select(SAMPLE_ID, ALTERATION, oncogenic, Highest_level)
  
  # Mutations
  # Filter for status
  if (status == "somatic") {
    mut_df <- filter(mut_df, Mutation_Status != "GERMLINE" | is.na(Mutation_Status) == T)
  } else if (status == "germline") {
    clin_df <- filter(clin_df, consent == "YES")
    mut_df <- filter(mut_df, Mutation_Status == "GERMLINE")
    mut_df <- mut_df[mut_df$SAMPLE_ID %in% clin_df$SAMPLE_ID,]
  }
  # Clean & Filter
  mut_df <- mut_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID) %>%
    dplyr::select(SAMPLE_ID, oncogenic, Highest_level) %>%
    filter(grepl("Oncogenic", oncogenic) == T, is.na(Highest_level) == F & Highest_level != "") %>%
    mutate(ALTERATION = "Mutation") %>%
    dplyr::select(SAMPLE_ID, ALTERATION, oncogenic, Highest_level)
  
  
  # rbind to create master alterations data frame
  # Filter for status for mutation data frame
  if (status == "somatic" | status == "both") {
    alt_final <- rbind(cna_df, fus_df)
    alt_final <- rbind(alt_final, mut_df)
  } else if (status == "germline") {
    alt_final <- mut_df
  }
  alt_final <- left_join(alt_final, data_freeze[,c("SAMPLE_ID", group_col)], by = "SAMPLE_ID")
  
  # Save
  write.table(alt_final, "actionability_master_alterations_table.txt", sep = "\t", row.names = F, quote = F)
  
  ########## optional select only the highest level ##########
  
  if (only_highest_level == T){
    alt_final <- alt_final %>% 
      left_join(dplyr::select(clin_df, SAMPLE_ID, HIGHEST_LEVEL), by = "SAMPLE_ID") %>%
      mutate_if(is.factor, as.character) %>%
      filter(HIGHEST_LEVEL == Highest_level)
  }
  
  ###########
  
  # Save
  write.table(alt_final, "actionability_master_alterations_highest_level_table.txt", sep = "\t", row.names = F, quote = F)
  
  # Calculate the percentage of each alteration by subtype
  prop_alteration_df <- as.data.frame(freq_dataframe(alt_final, group_col, "ALTERATION"))
  prop_alteration_df$freq[is.na(prop_alteration_df$freq)] <- 0
  prop_alteration_df$ALTERATION <- factor(prop_alteration_df$ALTERATION,
                                          levels = c("Amplification", "Deletion", "Fusion", "Mutation"))
  prop_alteration_df$group <- factor(prop_alteration_df[,group_col],
                                     levels = cancer_order_other)
  
  # Save
  write.table(prop_alteration_df, "actionability_alterations_barplot_table.txt", sep = "\t", row.names = F, quote = F)
  
  # Plot for ACTIONABLE ALTERATIONS
  alt_freq_bar_plot <- ggplot(prop_alteration_df, aes(y = freq, x = group, fill = ALTERATION)) +
    geom_col(position = position_stack(reverse = FALSE)) +
    ylab("Frequency") +
    labs(fill = "Actionable Alteration") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm"),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
          legend.justification="left",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-5)) +
    scale_y_continuous(limits=c(0, 1.00), expand = c(0, 0)) +
    ylab("Frequency") +
    scale_fill_manual(values = c("#A11111",  "#02488E", "#660066", "#037903"),
                      limits = c("Amplification", "Deletion", "Fusion", "Mutation"),
                      labels = c("Amplification", "Deletion", "Fusion", "Mutation")) +
    scale_x_discrete(limits = cancer_order_other)
  alt_freq_bar_plot
  
  return(alt_freq_bar_plot)
  
}

# Create actionability count barplot
action_count_barplot_fun <- function(clin_df, data_freeze, group_col,
                                     prop_level_df = "./actionability_levels_barplot_table.txt",
                                     status = c("somatic", "germline", "both"),
                                     consent_col,
                                     alt_final_df = "./actionability_master_alterations_table.txt",
                                     msi_tmb_status,
                                     msi_tmb_df){
  
  # Read in files
  prop_level_df <- read.delim(prop_level_df)
  alt_final <- read.delim(alt_final_df)
  clin_df <- read.delim(clin_df)
  data_freeze <- read.delim(data_freeze)
  
  # Filter for samples in data freeze and clean consent column
  data_freeze$SAMPLE_ID <- as.character(data_freeze$SAMPLE_ID)
  
  ######
  
  # Optional MSI/TMB addition
  if (msi_tmb_status == TRUE){
    msi_tmb_df <- read.delim(msi_tmb_path)
    msi_tmb_df <- msi_tmb_df %>%
      dplyr::select(SAMPLE_ID) %>% 
      mutate_if(is.factor, as.character) %>%
      mutate(Highest_level = "LEVEL_1_MSI-H_TMB-H") %>%
      distinct()
    data_freeze <- filter(data_freeze, !SAMPLE_ID %in% msi_tmb_df$SAMPLE_ID)
  } 
  
  #####
  
  # Clean
  clin_df <- clin_df[as.character(clin_df$SAMPLE_ID) %in% data_freeze$SAMPLE_ID,]
  colnames(data_freeze)[which(names(data_freeze) == consent_col)] <- "consent"
  
  # Set order
  cancer_order_other <- as.character(unique(prop_level_df[,c(group_col)]))
  
  # Create data frame that counts the number of actionable oncogenic alterations
  alt_final$alt_count <- 1
  alt_final <- dplyr::select(alt_final, SAMPLE_ID, alt_count)
  
  # Filter for status
  if (status == "germline") {
    clin_df <- clin_df[clin_df$SAMPLE_ID %in% as.character(filter(data_freeze, consent == "YES")$SAMPLE_ID),]
  }
  
  # Add in samples that don't have an actionable alteration
  alt_final_none <- as.data.frame(clin_df[,c("SAMPLE_ID")])
  colnames(alt_final_none)[1] <- "SAMPLE_ID"
  alt_final_none$alt_count <- 0
  alt_final <- rbind(alt_final, alt_final_none)
  alt_final <- aggregate(alt_count ~ SAMPLE_ID, alt_final, sum)
  
  # Add cancer subtypes to clinical data frame and create labels
  alt_final <- left_join(alt_final, data_freeze[,c("SAMPLE_ID", group_col)], by = "SAMPLE_ID")
  alt_final$label <- ifelse(alt_final$alt_count >= 3, "3+", alt_final$alt_count)
  
  # Calculate the percentage of each count by subtype
  prop_alt_count_df <- as.data.frame(freq_dataframe(alt_final, group_col, "label"))
  prop_alt_count_df$freq[is.na(prop_alt_count_df$freq)] <- 0
  
  # Set order
  prop_alt_count_df$label <- factor(prop_alt_count_df$label,  levels = c("0", "1", "2", "3+"))
  prop_alt_count_df$group <- factor(prop_alt_count_df[,group_col],
                                    levels = cancer_order_other)
  
  # Save
  write.table(prop_alt_count_df, "actionability_count_table.txt", sep = "\t", row.names = F, quote = F)
  
  # Number of alterations plot
  alt_per_num_prop_plot <- ggplot(prop_alt_count_df, aes(y = freq, x = group, fill = label)) +
    geom_col(position = position_stack(reverse = FALSE)) +
    ylab("Frequency") +
    labs(fill = "# of Actionable Alterations") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm"),
          axis.title.x = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
          legend.justification="left",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-5)) +
    scale_y_continuous(limits=c(0, 1.00), expand = c(0, 0))  +
    scale_fill_manual(values = c("#F7E690", "#F7AA14", "#E17202" ,"#701C5A"),
                      limits = c("0", "1", "2", "3+"),
                      labels = c("0", "1", "2", "3+")) +
    scale_x_discrete(limits = cancer_order_other)
  alt_per_num_prop_plot
  
  return(alt_per_num_prop_plot)
  
}

# Create actionability alterations main plot
action_main_fun <- function(cna_df, mut_df, fus_df, clin_df, data_freeze,
                            path_df,
                            tsg_list, fusion_list,
                            prop_level_df = "./actionability_levels_barplot_table.txt",
                            group_col,
                            consent_col,
                            alt_min = 1,
                            status = c("somatic", "germline", "both"),
                            gene_order,
                            only_highest_level = F,
                            msi_tmb_status,
                            msi_tmb_df,
                            include_oncogenic = F){
  
  # Read in data
  cna_df <- read.delim(cna_df)
  fus_df <- read.delim(fus_df)
  mut_df <- read.delim(mut_df)
  clin_df <- read.delim(clin_df)
  data_freeze <- read.delim(data_freeze)
  tsg_df <- read.delim(tsg_list, header = F)
  prop_level_df <- read.delim(prop_level_df)
  
  # Set order
  cancer_order_other <- as.character(unique(prop_level_df[,c(group_col)]))
  
  # Clean and filter data
  # Data freeze
  colnames(data_freeze)[which(names(data_freeze) == group_col)] <- "cancer_type"
  colnames(data_freeze)[which(names(data_freeze) == consent_col)] <- "consent"
  data_freeze <- data_freeze %>%
    mutate_if(is.factor, as.character)
  
  # Optional MSI/TMB addition
  if (msi_tmb_status == TRUE){
    msi_tmb_df <- read.delim(msi_tmb_path)
    msi_tmb_df <- msi_tmb_df %>%
      dplyr::select(SAMPLE_ID) %>% 
      mutate_if(is.factor, as.character) %>%
      mutate(Highest_level = "LEVEL_1_MSI-H_TMB-H") %>%
      distinct()
    data_freeze <- filter(data_freeze, !SAMPLE_ID %in% msi_tmb_df$SAMPLE_ID)
  } 
  
  # Clinical
  clin_df <- clin_df %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID)
  
  # CNA
  cna_df <- cna_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID),
           Highest_level = as.character(Highest_level)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID,
           grepl("Oncogenic", oncogenic)) 
  
  # Fusions
  fus_df <- fus_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID),
           Highest_level = as.character(Highest_level)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID,
           grepl("Oncogenic", oncogenic))
  
  # Mutations
  mut_df <- mut_df %>%
    dplyr::rename_all(recode, 
                      Tumor_Sample_Barcode = "SAMPLE_ID", 
                      HIGHEST_LEVEL = "Highest_level",
                      ONCOGENIC = "oncogenic") %>%
    mutate(SAMPLE_ID = as.character(SAMPLE_ID),
           Highest_level = as.character(Highest_level)) %>%
    filter(SAMPLE_ID %in% data_freeze$SAMPLE_ID,
           grepl("Oncogenic", oncogenic))
  
  # Set tumor suppresor list
  tumor_suppressor_list <- as.character(tsg_df$V1)
  
  # Make count data frame - consider somatic/germline/both
  if (status == "germline") {
    data_freeze <- filter(data_freeze, consent == "YES")
    clin_oncotree_freq <- as.data.frame(table(data_freeze$cancer_type))
    colnames(clin_oncotree_freq)[] <- c("cancer_type", "total_count")
    mut_df <- mut_df[mut_df$SAMPLE_ID %in% data_freeze$SAMPLE_ID,]
  } else if (status == "somatic") {
    clin_oncotree_freq <- as.data.frame(table(data_freeze$cancer_type))
    colnames(clin_oncotree_freq)[] <- c("cancer_type", "total_count")
  } else {
    clin_oncotree_freq <- as.data.frame(table(data_freeze$cancer_type))
    data_freeze_1 <- filter(data_freeze, consent == "YES")
    data_freeze_1$SAMPLE_ID <- as.character(data_freeze_1$SAMPLE_ID)
    clin_oncotree_freq_1 <- as.data.frame(table(data_freeze_1$cancer_type))
    clin_oncotree_freq <- left_join(clin_oncotree_freq, clin_oncotree_freq_1, by = "Var1")
    colnames(clin_oncotree_freq)[] <- c("cancer_type", "total_count", "germ_count")
    # Remove samples that have germline alterations but ARE NOT Part C consented
    remove_list <- intersect(filter(mut_df, Mutation_Status == "GERMLINE")$SAMPLE_ID,
                             filter(data_freeze, consent == "NO")$SAMPLE_ID)
    mut_df <- mut_df[!(mut_df$SAMPLE_ID %in% remove_list),]
  }
  
  
  # Create CNA data frame, combine with pathways and tumor suppresor list
  cna_df <- cna_df %>%
    inner_join(dplyr::select(data_freeze, SAMPLE_ID, cancer_type), by = "SAMPLE_ID") %>%
    dplyr::select(SAMPLE_ID, HUGO_SYMBOL, ALTERATION, LEVEL_1, LEVEL_2, LEVEL_3A,
                  LEVEL_3B, LEVEL_4, Highest_level, oncogenic, cancer_type) %>%
    distinct() %>%
    filter(is.na(Highest_level) == F) %>%
    mutate(ALTERATION = substring(ALTERATION, 1, 3)) %>%
    dplyr::select(SAMPLE_ID, HUGO_SYMBOL, ALTERATION, Highest_level, oncogenic, cancer_type) %>%
    dplyr::rename(sample_id = SAMPLE_ID, 
                  gene_symbol = HUGO_SYMBOL, 
                  alteration = ALTERATION, 
                  highest_level = Highest_level) %>%
    mutate(onco_type = ifelse(gene_symbol %in% tumor_suppressor_list, "tumor_suppresor", NA))
  
  # Create fusion data frame
  # Combine fusions where the hugo gene symbol is counted twice (impact and archer)
  fus_df <- fus_df %>%
    inner_join(dplyr::select(data_freeze, SAMPLE_ID, cancer_type), by = "SAMPLE_ID") %>%
    dplyr::select(SAMPLE_ID, Hugo_Symbol, Fusion, LEVEL_1, LEVEL_2, LEVEL_3A,
                  LEVEL_3B, LEVEL_4, Highest_level, oncogenic, cancer_type) %>%
    mutate_if(is.factor, as.character) %>%
    mutate(Fusion = gsub(" fusion", "", Fusion)) %>%
    mutate(Fusion = gsub(" - Archer", "", Fusion)) %>%
    rowwise() %>%
    mutate(Fusion = ifelse(grepl("intragenic", Fusion), Fusion,
                           paste(sort(unlist(strsplit(Fusion, "-", fixed = TRUE))), collapse = "-"))) %>%
    ungroup() %>%
    distinct()
  
  # If fusion list is provided, select the gene partner of interest based on the list
  if (missing(fusion_list) == FALSE) {
    # Read in fusion list
    fusion_list <- read.delim(fusion_list, header = F)
    fusion_list <- as.character(fusion_list$V1)
    fusion_list_collapse <- paste0("\\b", paste(fusion_list , collapse="\\b|\\b"), "\\b")
    # Filter for fusion list or full fusion name
    fus_df <- fus_df %>% 
      mutate_if(is.factor, as.character) %>%
      mutate(Fusion = ifelse(Hugo_Symbol %in% fusion_list, Hugo_Symbol,
                             ifelse(grepl(fusion_list_collapse, Fusion) == F, Fusion, "REMOVE"))) %>%
      filter(Fusion != "REMOVE") %>%
      mutate(Fusion = gsub("-intragenic", "", Fusion))
  }
  
  # Clean, add tumor suppresor columns
  fus_df <- fus_df %>%
    filter(Highest_level != "") %>%
    mutate(Alteration = "Fus") %>%
    dplyr::select(SAMPLE_ID, Fusion, Alteration, Highest_level, oncogenic, cancer_type) %>%
    dplyr::rename(sample_id = SAMPLE_ID,
                  gene_symbol = Fusion,
                  alteration = Alteration,
                  highest_level = Highest_level) %>%
    mutate(onco_type = ifelse(gene_symbol %in% tumor_suppressor_list, "tumor_suppresor", NA)) %>%
    distinct()
  
  # Collapse NTRK fusions
  # Other fusions can be added to this list moving forward
  fus_df <- fus_df %>%
    mutate_if(is.factor, as.character) %>%
    mutate(gene_symbol = ifelse(gene_symbol %in% c("NTRK1", "NTRK2", "NTRK3"), "NTRK1/2/3", gene_symbol)) %>%
    distinct()
  
  # Filter for mutation status
  if (status == "somatic") {
    mut_df <- filter(mut_df, Mutation_Status != "GERMLINE" | is.na(Mutation_Status) == T)
  } else if (status == "germline") {
    mut_df <- filter(mut_df, Mutation_Status == "GERMLINE")
  }
  
  # Mutation
  mut_df <- mut_df %>%
    inner_join(dplyr::select(data_freeze, SAMPLE_ID, cancer_type), by = "SAMPLE_ID") %>%
    dplyr::select(SAMPLE_ID, Hugo_Symbol, Variant_Type, LEVEL_1, LEVEL_2, LEVEL_3A,
                  LEVEL_3B, LEVEL_4, Highest_level, oncogenic, cancer_type, HGVSp_Short, Mutation_Status)
  # Add in oncogenic here if included
  if (include_oncogenic == T) {
    mut_df <- mut_df %>%
      mutate(ONCOGENIC = "ONCOGENIC")
  } 
  mut_df <- melt(mut_df, id.vars = c("SAMPLE_ID", "Hugo_Symbol", "Variant_Type", "Highest_level",
                                     "oncogenic", "cancer_type", "HGVSp_Short", "Mutation_Status"))
  
  # Aggregate by everything but strip for the highest level
  # This is just in case there is a gene alteration that has more than one level
  # Add pathways and tumor suppressor column
  # Remove duplicates if they are in the same pathway (use order of input df)
  mut_df <- mut_df %>%
    mutate_if(is.factor, as.character) %>%
    filter(value != "") %>%
    dplyr::select(-value, -Highest_level) %>%
    dplyr::rename(highest_level = variable) %>%
    filter(is.na(highest_level) == F) %>%
    mutate(highest_level == as.character(highest_level),
           Mutation_Status = ifelse(is.na(Mutation_Status) == T, "", Mutation_Status)) %>%
    group_by(SAMPLE_ID, Hugo_Symbol, Variant_Type, oncogenic, cancer_type, HGVSp_Short, Mutation_Status) %>%
    dplyr::summarise(highest_level = toString(highest_level)) %>%
    ungroup() %>%
    mutate(highest_level = gsub(",.*", "", highest_level),
           alteration = "Mut") %>%
    mutate(onco_type = ifelse(Hugo_Symbol %in% tumor_suppressor_list, "tumor_suppresor", NA)) %>%
    ###
    ### work in progress
    mutate(Hugo_Symbol = ifelse(Hugo_Symbol == "BRAF" & HGVSp_Short == "p.V600E", "BRAF_V600E", Hugo_Symbol)) %>%
    mutate(Hugo_Symbol = ifelse(Hugo_Symbol == "BRAF" & HGVSp_Short != "p.V600E", "BRAF_Other", Hugo_Symbol)) %>%
    ###
    ###
    distinct() %>%
    rename(gene_symbol = Hugo_Symbol) %>%
    mutate(gene_symbol = as.character(gene_symbol)) %>%
    dplyr::select(SAMPLE_ID, gene_symbol, alteration, highest_level, oncogenic, cancer_type, Mutation_Status, onco_type) %>%
    dplyr::rename(sample_id = SAMPLE_ID) %>%
    mutate(onco_type = ifelse(Mutation_Status == "GERMLINE", "germline", onco_type),
           gene_symbol = ifelse(Mutation_Status == "GERMLINE", paste0(gene_symbol, "*"), gene_symbol)) %>%
    dplyr::select(-Mutation_Status) %>%
    group_by(sample_id, gene_symbol, alteration, highest_level, oncogenic, cancer_type, onco_type) %>%
    dplyr::slice(1) %>%
    ungroup()
  
  # Filter for status
  # Combine CNA, FUS, and MUT - create final df
  if (status == "somatic" | status == "both") {
    gene_final_df <- rbind(cna_df, fus_df)
    gene_final_df <- rbind(gene_final_df, mut_df)
  } else if (status == "germline") {
    gene_final_df <- mut_df
  }
  
  # Optional include oncogenic alterations in plot
  if (include_oncogenic == T){
    gene_final_df <- gene_final_df %>%
      mutate(highest_level = ifelse((is.na(highest_level) == T | highest_level == "") &
                                      grepl("Oncogenic", oncogenic) == T, "ONCOGENIC", highest_level))
  } 
  
  # Combine all tumor suppressor alterations (del, mut, fus)
  # If the alteration is on a tumor suppresor, ignore alteration label
  # Clean up gene symbol, remove everything after the comma
  # Remove mutation label to clean up y axis
  gene_final_df <- gene_final_df %>%
    filter(is.na(highest_level) == F & highest_level != "") %>%
    mutate_if(is.factor, as.character) %>%
    mutate(onco_type = ifelse(is.na(onco_type) == T, "oncogene", onco_type)) %>%
    group_by(sample_id, gene_symbol, highest_level, cancer_type, onco_type) %>%
    dplyr::summarise(alteration = toString(alteration)) %>%
    ungroup() %>%
    mutate(alteration = as.character(alteration)) %>%
    mutate(alteration = ifelse(onco_type == "tumor_suppresor", "Del", alteration)) %>%
    mutate(alteration = gsub(",.*", "", alteration)) %>%
    mutate(gene_symbol_label = gsub(" Mut", "", paste0(gene_symbol, " ", alteration)))

  # Optional select only the highest level
  if (only_highest_level == T){
    colnames(clin_df)[which(names(clin_df) == "SAMPLE_ID")] <- "sample_id"
    gene_final_df <- gene_final_df %>%
      left_join(dplyr::select(clin_df, sample_id, HIGHEST_LEVEL), by = "sample_id") %>%
      mutate_if(is.factor, as.character) %>%
      filter(HIGHEST_LEVEL == highest_level)
  }
  
  # Manual alterations
  ###
  ### work in progress
  gene_final_df <- gene_final_df %>% 
    mutate_if(is.factor, as.character) %>%
    mutate(gene_symbol = ifelse(gene_symbol %in% c("BRCA1", "BRCA2"), "BRCA1/2", gene_symbol),
           gene_symbol_label = ifelse(gene_symbol == "BRCA1/2", "BRCA1/2 Del", gene_symbol_label)) %>%
    mutate(gene_symbol = ifelse(gene_symbol %in% c("CHEK1", "CHEK2"), "CHEK1/2", gene_symbol),
           gene_symbol_label = ifelse(gene_symbol == "CHEK1/2", "CHEK1/2 Del", gene_symbol_label)) %>%
    mutate(gene_symbol = ifelse(gene_symbol %in% c("TSC1", "TSC2"), "TSC1/2", gene_symbol),
           gene_symbol_label = ifelse(gene_symbol == "TSC1/2", "TSC1/2 Del", gene_symbol_label)) %>% 
    distinct()
  ###
  ###
  
  # Calculate the percentage of each count by subtype
  # Only select the highest level
  prop_main_plot_df <- gene_final_df %>%
    group_by(cancer_type, gene_symbol_label, highest_level) %>%
    dplyr::summarise(n = n()) %>%
    ungroup() %>%
    left_join(clin_oncotree_freq, by = "cancer_type") %>%
    dplyr::mutate(freq = n /total_count) %>%
    group_by(cancer_type, gene_symbol_label) %>%
    arrange(highest_level) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::mutate(percentage = 100*freq,
                  label_text = round(percentage, 0),
                  label_text = ifelse(percentage > 0 & percentage < 1, " ", label_text))
  
  # Optional add pathway if provided, if not use it to set gene list
  if (missing(path_df) == T) {
    path_df <- gene_final_df  %>%
      left_join(prop_main_plot_df) %>%
      dplyr::select(gene_symbol_label, highest_level, cancer_type, percentage) %>%
      distinct() %>%
      group_by(gene_symbol_label, highest_level) %>%
      mutate(count = n()) %>%
      ungroup() %>%
      arrange(highest_level, desc(count), desc(percentage), gene_symbol_label) %>%
      group_by(gene_symbol_label) %>%
      dplyr::slice(1) %>%
      ungroup() %>%
      arrange(highest_level, desc(count), desc(percentage), gene_symbol_label) %>%
      mutate(pathway = row_number()) %>%
      dplyr::select(gene_symbol_label, pathway)
    gene_final_df <- gene_final_df %>% left_join(path_df)
  } else {
    path_df <- read.delim(path_df)
    colnames(path_df)[] <- c("gene_symbol", "pathway")
    gene_final_df <- gene_final_df %>% left_join(path_df)
  }
  
  # Add germline label if figure includes both somatic and germline
  if (status == "both") {
    prop_main_plot_df$freq <- ifelse(grepl("\\*",prop_main_plot_df$gene_symbol_label) == TRUE,
                                     prop_main_plot_df$n/prop_main_plot_df$germ_count,
                                     prop_main_plot_df$freq)
  }
  
  # Add pathways
  prop_main_plot_df <- prop_main_plot_df %>%
    left_join(dplyr::select(gene_final_df, gene_symbol_label, gene_symbol), by = "gene_symbol_label") %>%
    left_join(dplyr::select(gene_final_df, gene_symbol, cancer_type, highest_level, pathway, onco_type),
              by = c("gene_symbol", "cancer_type", "highest_level")) %>%
    group_by(gene_symbol, gene_symbol_label, cancer_type, pathway, onco_type, n, total_count, percentage, freq, label_text) %>%
    dplyr::summarise(highest_level = toString(highest_level)) %>%
    ungroup() %>%
    mutate(highest_level_label = gsub(",.*", "", highest_level)) %>%
    dplyr::select(-highest_level) %>%
    dplyr::arrange(pathway, gene_symbol, highest_level_label, cancer_type)
  
  # Only keep rows where at least one subtype meets the percetage threshold (alt_min)
  prop_main_plot_df_filter <- prop_main_plot_df %>%
    dplyr::select(gene_symbol_label, percentage) %>%
    group_by(gene_symbol_label) %>%
    filter(percentage == max(percentage)) %>%
    filter(percentage < alt_min)
  prop_main_plot_df <- prop_main_plot_df %>%
    filter(!gene_symbol_label %in%prop_main_plot_df_filter$gene_symbol_label) %>%
    mutate(cancer_type = factor(cancer_type, levels = cancer_order_other))
  
  # Set gene order manually
  if (missing(gene_order) == F) {
    gene_order <- read.delim(gene_order, header = F)
    gene_order <- as.data.frame(gene_order[rep(seq_len(nrow(gene_order)), each = 2), ])
    colnames(gene_order)[] <- c("gene_symbol")
    gene_order <- gene_order %>%
      mutate_if(is.factor, as.character) %>%
      mutate(order = seq(1:nrow(gene_order)),
             gene_symbol = ifelse(order %% 2 == 0, paste0(gene_symbol, "*"), gene_symbol))
    prop_main_plot_df <- prop_main_plot_df %>%
      left_join(gene_order, by = "gene_symbol") %>%
      dplyr::arrange(order)
  }
  
  # Get text color order
  text_tsg_col <- prop_main_plot_df %>%
    dplyr::select(gene_symbol_label, onco_type) %>%
    distinct() %>%
    dplyr::arrange(onco_type) %>%
    group_by(gene_symbol_label) %>%
    dplyr::summarise(onco_type = toString(onco_type)) %>%
    ungroup() %>%
    mutate(col = ifelse(onco_type != "tumor_suppresor", ifelse(onco_type == "oncogene", "#7E1116", "#4F0043"), "#191A57")) %>%
    group_by(gene_symbol_label, col) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(gene_symbol_label = factor(gene_symbol_label, levels = unique(prop_main_plot_df$gene_symbol_label))) %>%
    dplyr::arrange(gene_symbol_label)
  
  # Write out data frame
  write.table(prop_main_plot_df, "actionability_main_plot_data.txt", sep = "\t", quote = F, row.names = F)
  
  # Create main plot
  action_tile_plot_all <- ggplot(data = prop_main_plot_df, aes(x = cancer_type, y = gene_symbol_label)) +
    geom_tile(aes(fill = highest_level_label)) +
    geom_text(aes(label = label_text), colour = "white", size = 2) +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6), # colour = rev(text_tsg_col$col)),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 8),
          plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.justification="left",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,0,-10,-5)) +
    geom_vline(xintercept=seq(1.5, length(levels(prop_main_plot_df$cancer_type))-0.5, 1),
               lwd=0.25, colour="gray80") +
    geom_hline(yintercept=seq(1.5, length(unique(prop_main_plot_df$gene_symbol_label))-0.5, 1),
               lwd=0.25, colour="gray80") +
    scale_fill_manual(values = c("#88E281","#33A02C", "#1F78B4", "#984EA3", "#BE98CE", "#a8a8a8", "#ffdab9", "gray90"),
                      limits = c("LEVEL_1_MSI-H_TMB-H","LEVEL_1", "LEVEL_2", "LEVEL_3A", "LEVEL_3B", "LEVEL_4", "ONCOGENIC", "NO_LEVEL"),
                      labels = c("LEVEL 1 MSI/TMB-H","LEVEL 1", "LEVEL 2", "LEVEL 3A", "LEVEL 3B", "LEVEL 4", "ONCOGENIC", "NO LEVEL")) +
    scale_y_discrete(limits = rev(unique(prop_main_plot_df$gene_symbol_label)),
                     labels = gsub("_", " ", rev(unique(prop_main_plot_df$gene_symbol_label))),
                     expand = c(0,0)) +
    scale_x_discrete(limits = cancer_order_other) +
    labs(fill = "Highest Level\nof Evidence") +
    guides(fill = guide_legend(override.aes = list(size = 1)))
  
  return(action_tile_plot_all)
}


# Create actionability TMB-H & MSI-H main plot add-on
action_main_msi_tmb_fun <- function(clin_df, 
                                    data_freeze,
                                    group_col,
                                    prop_level_df = "./actionability_levels_barplot_table.txt",
                                    msi_tmb_df){
  # Read in data
  data_freeze <- read.delim(data_freeze)
  clin_df <- read.delim(clin_df)
  msi_tmb_df <- read.delim(msi_tmb_df)
  prop_level_df <-  read.delim(prop_level_df)
  
  # Set order
  cancer_order_other <- as.character(unique(prop_level_df[,c(group_col)]))
  
  # Get MSI/TMB frequency
  aty_alt_df <- msi_tmb_df %>% 
    mutate_if(is.factor, as.character) %>%
    dplyr::select(SAMPLE_ID, ALTERATION, HIGHEST_LEVEL) %>%
    right_join(data_freeze) %>%
    dplyr::select(CANCER_TYPE, ALTERATION) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(ALTERATION = ifelse(is.na(ALTERATION) == T, "NONE", ALTERATION))
  aty_alt_df <- aty_alt_df %>%
    left_join(dplyr::count(dplyr::select(group_by(aty_alt_df, CANCER_TYPE), CANCER_TYPE))) %>%
    dplyr::rename(total_count = n) %>%
    mutate_if(is.character, as.factor) %>%
    group_by(ALTERATION, CANCER_TYPE, total_count) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(freq = n/total_count,
                  percentage = 100*freq,
                  label_text = as.character(round(percentage, 0)),
                  label_text_final = ifelse(percentage > 0 & percentage < 1, "", label_text)) %>%
    filter(ALTERATION != "NONE") %>%
    dplyr::mutate(ALTERATION = factor(ALTERATION, levels = c("TMB-H", "MSI-H")),
                  label = "MSI-H & TMB-H") 
  
  # Plot
  aty_alt_tile_plot_all <- ggplot(data = aty_alt_df, aes(x = CANCER_TYPE, y = ALTERATION)) +
    geom_tile(aes(fill = label)) +
    geom_text(aes(label = label_text_final), colour = "black", size = 2) +
    theme(panel.grid.major = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),
          panel.background = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.05, 0.05, 0.1, 0.05), "cm"),
          legend.justification="left",
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(-10,0,-10,-5),
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.key.size = unit(0.4, "cm")) +
    geom_vline(xintercept=seq(1.5, length(levels(aty_alt_df$CANCER_TYPE))-0.5, 1),
               lwd=0.5, colour="white") +
    geom_hline(yintercept=seq(1.5, length(unique(aty_alt_df$ALTERATION))-0.5, 1),
               lwd=0.25, colour="white") +
    scale_fill_manual(values = c("#88E281"),
                      limits = c("MSI-H & TMB-H"),
                      labels = c("MSI-H & TMB-H")) +
    scale_y_discrete(limits = levels(aty_alt_df$ALTERATION),
                     expand = c(0,0)) +
    scale_x_discrete(limits = cancer_order_other)
  
  return(aty_alt_tile_plot_all)
  
}



#--------------