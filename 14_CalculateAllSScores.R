# 14: Calculate All S Scores
# This script will calculate all organokine S score within a source-target tissue pairs
# It will then output a list of prioritized organokines for each tissue pair, and their associated LVs


library(readr)
library(dplyr)
library(stringr)


source("../utils/plier_util.R")

in_dir <- "12_out"

out_dir <- "14_out"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

S_cap <- 1e-6

all_gtex_tissues <- readr::read_tsv("../data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
all_gtex_tissues <- unique(all_gtex_tissues$SMTSD)

all_correlations_filtered_reciprocal_removed_filtered <- readr::read_tsv(file.path(in_dir, "all_correlations_filtered.tsv"))


all_S <- data.frame()

for (source_todo in all_gtex_tissues){
  for(target_todo in all_gtex_tissues){
    if (source_todo != target_todo){
      
      print(paste(source_todo, target_todo))
      
      # Subset the correlations for this source-target tissue pair
      subset_correlations <- all_correlations_filtered_reciprocal_removed_filtered |>
        dplyr::filter(source == source_todo, target == target_todo)
      
      # If there are no correlations, skip this pair
      if (nrow(subset_correlations) == 0) {
        next
      }
      
      # Calculate S score for each Source_gene (source secreted protein)
      S_scores <- subset_correlations %>%
        dplyr::group_by(Source_gene) %>%
        dplyr::summarize(S = sum(-log10(pmax(p_adjusted, S_cap)))) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(desc(S))
      
      # Calculate Z score of S
      S_scores$Z <- (S_scores$S - mean(S_scores$S)) / sd(S_scores$S)
      
      # Turn the Z scores into p values
      S_scores$p_S <- pnorm(S_scores$Z, lower.tail = FALSE)
      
      # Add source and target tissue information
      S_scores$source <- source_todo
      S_scores$target <- target_todo
      
      # Combine with all_S data frame
      all_S <- rbind(all_S, S_scores)
      
    }
  }
}

readr::write_tsv(all_S, file.path(out_dir, "all_S_scores.tsv"))

# Read back
# all_S <- readr::read_tsv(file.path(out_dir, "all_S_scores.tsv"))

# Now filter all_correlations_filtered_reciprocal_removed_filtered by prioritized organokines (Z > 2)
all_S_filtered <- all_S %>%
  dplyr::filter(Z >= 2) %>%
  dplyr::select(Source_gene, source, target, S, Z, p_S)

all_correlations_filtered_reciprocal_removed_filtered_Z <- all_S_filtered |>
  dplyr::left_join(all_correlations_filtered_reciprocal_removed_filtered, 
                        by = c("Source_gene", "source", "target"))

readr::write_tsv(all_correlations_filtered_reciprocal_removed_filtered_Z, 
                 file.path(out_dir, "all_correlations_filtered_Z_prioritized.tsv"))
