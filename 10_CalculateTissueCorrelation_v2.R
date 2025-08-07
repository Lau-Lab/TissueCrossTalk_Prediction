library(tidyverse)
library(hpar)
library(DESeq2)
library(Hmisc)
library(reshape2)
library(beepr)

time_0 <- Sys.time()
#Requirements
#library(devtools)
#install_github("wgmao/PLIER")


main_path <- '.'
plier_path <- file.path(main_path, '03_out_old')
gtex_path <- file.path(main_path, '01_out_old')
save_path <- file.path(main_path, '10_out')
if(!dir.exists(save_path)){
  dir.create(save_path)
}

source(file.path(main_path, '..', 'utils', 'plier_util.R'))


file <- c("Correlation_data.txt", "Correlation_results.csv", "Ssec_results.csv")

#Read files
print("Reading in PLIER model ...")
PLIER_model <- readRDS(file.path(plier_path, 'GTExv8VSTcomBATAgeCorrNoSecretome_canonicalPathwaysBloodCellSVM_multiplier_v2.Rds'))
print("Reading in GTEx data ...")
GTEx_data <- readRDS(file.path(gtex_path, 'gtexv8_vst_combat4_dds.Rds'))


gsub_pattern <- 'GTEX-([A-Z0-9]*)-.*$'
gsub_replacement <- "\\1"

S_cap <- 1e-6

rowNormZ <- function(x){
  standard_deviation <- apply(x,1,sd)
  mean <- apply(x,1,mean);
  x=sweep(x,1,mean)
  x=sweep(x,1,standard_deviation,"/")
  return(x)
}


tissue_variable <- 'SMTSD'
# SMTS - General tissue types
# SMTSD - More focused tissue types
# To do - Find a way to check what group type was used, insert it as a global variable and have 01_filter data read it, currently it is fixed to tissue_variable

data_annotation <- colData(GTEx_data) |>as.data.frame()
data_annotation$SUBJECT <- gsub(gsub_pattern, gsub_replacement, data_annotation$SAMPID)
all_tissues <- unique(data_annotation[[tissue_variable]])

secreted_to_blood_annotation <- hpaSecretome() |>
  dplyr::filter(Secretome.location == "Secreted to blood") |>
  dplyr::select('Gene.name', 'Secretome.location') 


counter <- 0
#Loop through tissue pairs
print("Beginning loop ...")
for(source_tissue in all_tissues){
  for(target_tissue in all_tissues){
    start_time <- Sys.time()
    #Avoid redundancies
    if(source_tissue != target_tissue){
      counter <- counter + 1
      print(paste0("counter: ", counter, " of ", length(all_tissues)*(length(all_tissues)-1)))
      
      #Save folder items
      save_folder <- file.path(save_path, paste0(stringr::str_replace_all(source_tissue, " ", ""),
                                                 "_to_",
                                                 stringr::str_replace_all(target_tissue, " ", "")))
      #Create folder if it doesn't exist, otherwise move on
      print(paste("Checking for folder:", basename(save_folder)))
      if(!dir.exists(save_folder)){
        print(paste(basename(save_folder), "folder does not exist, creating folder"))
        dir.create(save_folder)
      }else{
        print(paste(basename(save_folder), "folder exists, checking for files ..."))
      }
      #If folder already exists, check to see if all files are present
      check_file <- all(file.exists(file.path(save_folder, file)))
      skip_iteration <- FALSE
      if(check_file){
        print("All files exists, moving on ...")
        skip_iteration <<- TRUE
      }else{
        print("Files are missing, creating files ...")
      }
      
      #Skip iteration if folder and file are present
      if(skip_iteration){
        print("Skipping")
        next
      }
      
      print("1) Creating annotations")
      select_annotations <- data_annotation %>% 
        dplyr::select(sample_id = SAMPID, subject_id = SUBJECT, tissue = !!sym(tissue_variable), SMAFRZE) %>%
        dplyr::filter(SMAFRZE == "RNASEQ", tissue %in% c(source_tissue, target_tissue)) %>%
        dplyr::distinct(tissue, subject_id, .keep_all = T)
      
      both_tissues_annotations <- select_annotations %>%
        dplyr::group_by(subject_id) %>%
        dplyr::filter(n_distinct(tissue) == 2) %>%
        dplyr::ungroup()
      
      target_annotations <- both_tissues_annotations %>%
        dplyr::filter(tissue == target_tissue)
      source_annotations <- both_tissues_annotations %>%
        dplyr::filter(tissue == source_tissue)
      
      print(paste("2) Processing target tissue:", target_tissue))

      target.data <- assay(GTEx_data)[ , target_annotations$sample_id]
      print(paste("Starting number of target genes:", nrow(target.data)))
      print(paste("Starting number of target samples:", ncol(target.data)))
      
      if(length(target_annotations$sample_id) < 5){
        print("Insufficient number of samples to calculate a correlation matrix, moving on ...")
        next
      } else if(ncol(target.data) < 5){
        print("Insufficient number of samples to calculate a correlation matrix, moving on ...")
        next
      }
      
      print("Removing zero variance genes")
      target.data <- target.data[-caret::nearZeroVar(t(target.data)), ]
      print(paste("Remaining number of target genes:", nrow(target.data)))
      print(paste("Remaining number of target samples:", ncol(target.data)))
      print("Filtering for genes retained within tissue")
      target.data <- target.data[!rownames(target.data) %in% secreted_to_blood_annotation$Gene.name, ]
      print(paste("Generating B matrix"))
      target.data <- GetNewDataB(target.data, PLIER_model) #Function loaded in from plier_util.R
      colnames(target.data) <- gsub(gsub_pattern, gsub_replacement, colnames(target.data))
      rownames(target.data) <- paste0("Target_", rownames(target.data))
      ################################################################################
      print(paste("3) Processing source tissue:", source_tissue))
      source.data <- assay(GTEx_data)[ ,source_annotations$sample_id]
      print(paste("Starting number of source genes:", nrow(source.data)))
  
      #print("Removing zero variance genes")
      #Removing zero variance genes is important for machine learning purposes,
      #which makes sense when we are creating a new B matrix for the target genes
      #but since we are selecting genes based on whether it is secreted from the source tissue or not, 
      #I felt that the remove zero var component for the source was not necessary
#      source.data <- source.data[-caret::nearZeroVar(t(source.data)), ]
      print("Filtering for genes secreted to blood")
      source.data <- source.data[rownames(source.data) %in% secreted_to_blood_annotation$Gene.name, ]
      print(paste("Remaining number of source genes:", nrow(source.data)))
      colnames(source.data) <- gsub(gsub_pattern, gsub_replacement, colnames(source.data))
      rownames(source.data) <- paste0("Source_", rownames(source.data))
      
      combined.data <- rbind(source.data, target.data)
      
      print("4) Creating correlation matrices")
      
      tryCatch({
        correlation_matrices <- Hmisc::rcorr(t(source.data), t(target.data))
      }, error = function(e){
        message("Error encountered: ", e$message)
        skip_iteration <<- TRUE
      })
      
      if(skip_iteration){
        print("Moving on ...")
        next
      }
      
      print("Summarizing correlations")
      correlations_p <- melt(correlation_matrices$P)
      correlations_r <- melt(correlation_matrices$r)
      
      colnames(correlations_p) <- c("Source_gene", "Target_gene", "p_value")
      colnames(correlations_r) <- c("Source_gene", "Target_gene", "r_correlation")
      
      correlations_summary <- merge(correlations_r, correlations_p, by = c("Source_gene", "Target_gene")) %>%
        filter(str_detect(Source_gene, "Source_"),
               str_detect(Target_gene, "Target_")) %>%
        mutate(p_adjusted = p.adjust(p_value, method = 'bonferroni'))
      
      print("5) Scoring correlation matrices")
      
      unique_sources <- unique(correlations_summary$Source_gene)
      unique_targets <- unique(correlations_summary$Target_gene)
      
      correlation_matrix_r <- correlation_matrices$r[(rownames(correlation_matrices$r) %in% unique_sources),
                                                   (colnames(correlation_matrices$r) %in% unique_targets)]
      correlation_matrix_p <- correlation_matrices$P[(rownames(correlation_matrices$P) %in% unique_sources),
                                                   (colnames(correlation_matrices$P) %in% unique_targets)]
      
      ################################################################################
      ############################ Z score function ##################################
      ################################################################################
      
      z_score <- function(x){
        return((x-mean(x))/sd(x))
      }
      
      ################################################################################
      ################################ Scoring by Ssec ###############################
      ################################################################################
      
      n_target_genes <- ncol(correlation_matrix_p)
      
      p_scores_list <- rowSums(-log10(pmax(correlation_matrix_p, S_cap)))
      r_scores_list <- rowSums(abs(correlation_matrix_r))
      
      Ssec_p <- data.frame(Source_gene = names(p_scores_list), p_score = p_scores_list) %>%
        dplyr::mutate(p_Ssec = p_score/n_target_genes,
                      p_Ssec_z = z_score(p_Ssec)) 
      
      Ssec_r <- data.frame(Source_gene = names(r_scores_list), r_score = r_scores_list) %>%
        dplyr::mutate(r_Ssec = r_score/n_target_genes,
                      r_Ssec_z = z_score(r_Ssec)) 
      
      Ssec <- merge(Ssec_p, Ssec_r, by = "Source_gene")
      
      
      
      
      
      ################################################################################
      ################################Write files#####################################
      ################################################################################
      print("6) Writing files")
      
      write.table(combined.data, file.path(save_folder, file[1]), sep = "\t", row.names = TRUE)
      write.csv(correlations_summary, file.path(save_folder, file[2]), row.names = FALSE)
      write.csv(Ssec, file.path(save_folder, file[3]), row.names = FALSE)
      
      end_time <- Sys.time()
      
#      elasped_time <- round((end_time - start_time), 2)
      elasped_time <- difftime(start_time, end_time, unit = "mins")
      print(paste0("Elasped time ", round(elasped_time, 1), " minutes"))
      
    }
  }
}


beep(sound = 3)
time_fin <- Sys.time()
run_time <- difftime(time_0, time_fin, unit = "hours")
print(paste0("Script run time ", round(run_time), 1), " hours")

