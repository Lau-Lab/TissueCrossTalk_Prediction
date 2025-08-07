# 12: Read Interorgan Crosstalk
# This script will read the inter-organ correlation results for all the significant co-expression between 
# source secretomes and tissue latent variables in each of 54 tissuesand generate tables for all correlations



library(dplyr)
library(readr)
library(stringr)


source("../utils/plier_util.R")

in_dir <- "10_out"

out_dir <- "12_out"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

S_cap <- 1e-6


all_gtex_tissues <- readr::read_tsv("../data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
all_gtex_tissues <- unique(all_gtex_tissues$SMTSD)


# Read in all the correlation results
all_correlations <- data.frame()

for (source_todo in all_gtex_tissues){
  for(target_todo in all_gtex_tissues){
    if (source_todo != target_todo){
      print(paste(source_todo, target_todo))
      
      matrix_write_path = file.path(in_dir, 
                                    paste0(stringr::str_replace_all(source_todo, " ", ""), 
                                           "_to_", 
                                           stringr::str_replace_all(target_todo, " ", "")),
                                    "Correlation_results.csv"
      )
      

      # If the file exists, read it and add to all_correlations
      if(file.exists(matrix_write_path)){
        
        corr_result <- readr::read_csv(matrix_write_path, show_col_types = F)
        
        # Filter by p value, the padj is already bonferroni adjusted by the number of genes/LV comparison. Here we should divide it further by the number 
        # of tissue pairs (54 x 54), so the final threshold is 0.01/(54*54) = 3.4e-6
        corr_result <- corr_result |> dplyr::filter(p_adjusted < 0.01) #/(54*54))
        
        corr_result$source <- source_todo
        corr_result$target <- target_todo
        all_correlations <- rbind(all_correlations, corr_result)
        
      }
    }
  }
}

readr::write_tsv(all_correlations, file.path(out_dir, "all_correlations.tsv"))

# Read back
all_correlations <- readr::read_tsv(file.path(out_dir, "all_correlations.tsv"))
all_correlations <- all_correlations |> dplyr::filter(p_adjusted < 0.005)

all_correlations_filtered <- all_correlations |> 
  dplyr::filter(!is.na(Source_gene), !is.na(Target_gene)) |> 
  dplyr::select(source, target, Source_gene, Target_gene, p_adjusted, r_correlation)

# Remove reciprocal correlations where a source secreted protein (Source_gene) and target LV (Target_gene) are significantly correlated in both directions
# This is to avoid double counting
all_correlations_filtered$forward_direction <- paste(all_correlations_filtered$Source_gene, all_correlations_filtered$Target_gene, all_correlations_filtered$source, all_correlations_filtered$target, sep = "_")
all_correlations_filtered$reverse_direction <- paste(all_correlations_filtered$Source_gene, all_correlations_filtered$Target_gene, all_correlations_filtered$target, all_correlations_filtered$source, sep = "_")
all_correlations_filtered_reciprocal_removed <- all_correlations_filtered |> dplyr::filter(!(forward_direction %in% reverse_direction))

# Get the number of source-LV pairs that recur in multiple tissue comparisons
source_lv_frequency <- all_correlations_filtered_reciprocal_removed |> 
  dplyr::group_by(Source_gene, Target_gene) |> dplyr::summarize(n = n_distinct(source, target)) |> arrange(desc(n))

# hist(source_lv_frequency$n, breaks=100)

# Remove source_LV pairs that appear in more than 100 source-target tissue pairs
source_lv_frequency_filtered <- source_lv_frequency|> dplyr::filter(n <= 100)

# Remove these from the all_correlations_filtered_reciprocal_removed
all_correlations_filtered_reciprocal_removed_filtered <- all_correlations_filtered_reciprocal_removed |>
  dplyr::filter((paste(Source_gene, Target_gene) %in% paste(source_lv_frequency_filtered$Source_gene, source_lv_frequency_filtered$Target_gene)))

# Write this out to a file
readr::write_tsv(all_correlations_filtered_reciprocal_removed_filtered, file.path(out_dir, "all_correlations_filtered.tsv"))


