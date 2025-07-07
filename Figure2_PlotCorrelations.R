library(tidyverse)

#correlation_matrix <- read.delim('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Adipose - Subcutaneous_to_Heart - Left Ventricle/Correlation_data.txt', sep = "\t")
#correlation_file <- read.csv('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Adipose - Subcutaneous_to_Heart - Left Ventricle/Correlation_results.csv')

correlation_matrix <- read.delim('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Liver_to_Heart - Left Ventricle/Correlation_data.txt', sep = "\t")
correlation_file <- read.csv('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Liver_to_Heart - Left Ventricle/Correlation_results.csv')

#gene <- "ADIPOQ"
gene <- "AHSG"

select_correlation <- correlation_file %>%
  mutate(Target_gene = gsub("Target", "LV", Target_gene),
         Source_gene = gsub("Source_","",Source_gene)) %>%
  filter(Source_gene == gene,
         p_adjusted <= 0.05,
         grepl(",\\s*.+$", Target_gene)) %>%
  #filter(grepl("KEGG_CYSTEINE_AND_METHIONINE_METABOLISM", Target_gene, ignore.case = TRUE)) %>%
  arrange(desc(abs(r_correlation)))

rownames(correlation_matrix) <- gsub("Source_", "", rownames(correlation_matrix))
rownames(correlation_matrix) <- gsub("Target", "LV", rownames(correlation_matrix))

#group <- c(21, 85, 136)
group <- c(9, 14, 83)
for(i in group){  
  correlation_matrix_select <- t(correlation_matrix[rownames(correlation_matrix) %in% c(select_correlation$Source_gene[i], 
                                                                                        select_correlation$Target_gene[i]), ])
  
  g <- ggplot(data = correlation_matrix_select, 
              mapping = aes(x = correlation_matrix_select[,1], 
                            y = correlation_matrix_select[,2])) +
    geom_point(size = 1) +
    geom_smooth(method = "lm") +
    labs(title = paste0(colnames(correlation_matrix_select)[1], " to ", colnames(correlation_matrix_select)[2], " correlation"),
         x = colnames(correlation_matrix_select)[1],
         y = colnames(correlation_matrix_select)[2]) +
    theme(aspect.ratio = 0.619,
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 8, color = 'black'),
          axis.text.y = element_text(size = 8, color = 'black'),
          axis.title.x = element_text(size = 8, color = 'black'),
          axis.title.y = element_text(size = 8, color = 'black'), 
          panel.grid = element_blank(),
          axis.line = element_line(color = "black"),
          panel.background = element_blank()) +
    annotate("Text", x = min(correlation_matrix_select[,1]), 
             y = (max(correlation_matrix_select[,2]) - 0.05*(max(correlation_matrix_select[,2]))), 
             label = paste("R = ", round(select_correlation$r_correlation[i], 3), 
                           "\np(adjusted) = ", 
                           sprintf("%.2e", select_correlation$p_adjusted[i])),
                           size = 3, 
                           hjust = 0)
  print(g)
  
#  ggsave(filename = paste0(colnames(correlation_matrix_select)[1], "_to_", colnames(correlation_matrix_select)[2], ".png"), g, width = 6, height = 6, dpi = 300)
}
