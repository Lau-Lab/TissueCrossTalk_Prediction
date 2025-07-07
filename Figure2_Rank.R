library(tidyverse)
library(ggrepel)

#correlation_file <- read.csv('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Adipose - Subcutaneous_to_Heart - Left Ventricle/Correlation_results.csv')
correlation_file <- read.csv('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Liver_to_Heart - Left Ventricle/Correlation_results.csv')
#correlation_file <- read.csv('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Heart - Left Ventricle_to_Adipose - Subcutaneous/Correlation_results.csv')
#correlation_file <-read.csv('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Adipose - Subcutaneous_to_Muscle - Skeletal/Correlation_results.csv')
#correlation_file <-read.csv('/Volumes/Calvin/TissueCrosstalk/Tissue_Pairs/Adipose - Subcutaneous_to_Liver/Correlation_results.csv')

z_score <- function(x){
  return((x-mean(x))/sd(x))
}

min_max_norm <- function(x){
  return((x-min(x))/(max(x) - min(x)))
}

select_correlation <- correlation_file %>%
#  filter(Source_gene == "Source_ADIPOQ",
  filter(Source_gene == "Source_AHSG",
         grepl(",\\s*.+$", Target_gene)) %>%
  mutate(
    Source_gene = gsub("Source_","", Source_gene),
    Target_gene = gsub("Target", "LV", Target_gene),
    min_max_r = min_max_norm(abs(r_correlation)),
    p_score = -log(p_adjusted + 1e-10),
    r_score = -log(abs(r_correlation))
  ) %>%
  arrange(desc(min_max_r)) %>%
  mutate(
    rank = row_number() 
  )

specific_correlations <- select_correlation %>%
  filter(
    p_adjusted <= 0.05) %>%
#  filter(grepl("LIPID|INSULIN|DNA", Target_gene, ignore.case = TRUE)) %>%
#  filter(grepl("APOPTOSIS|CALCIUM|SPLICEOSOME", Target_gene, ignore.case = TRUE)) %>%
  arrange(p_adjusted)

#ADIPOQ Adipose - Sucutaneous ~ Heart - Left Ventricle:
##Chose LV_214,REACTOME_IMMUNE_SYSTEM, LV_383,KEGG_CALCIUM_SIGNALING_PATHWAY, LV_212,SVM Macrophages M2
#select_correlation[c(21, 85, 136), ]
#AHSG Liver - Heart - Left Ventricle:


#group <- select_correlation[c(21, 85, 136), ]
group <- select_correlation[c(9, 14 ,83), ]
g <- ggplot(data = select_correlation,
            mapping = aes(x = rank, y = p_score)) +
  labs(#title = paste0(select_correlation$Source_gene[1], ":\nAdipose(Subcutaneous) ~ Heart(Left Ventricle)"),
       title = paste0(select_correlation$Source_gene[1], ":\nLiver ~ Heart(Left Ventricle)"),
       x = "rank",
       y = "p score") +
  geom_hline(yintercept = -log(0.05), 
             color = 'black', 
             alpha = 0.25,
             linewidth = 1,
             linetype = "dashed") +
  geom_point(size = 1, 
             alpha = 0.2) +
  geom_point(data = group,
             #data = select_correlation[c(21, 85, 136), ],
             aes(color = Target_gene),
             size = 2,
             alpha = 0.7,
             color = c("blue", "violet", "red")) +
  geom_text_repel(data = group,
                  #data = select_correlation[c(21, 85, 136), ],
                  aes(label = Target_gene),
                  size = 2.5,
                  box.padding = 0.5,
                  point.padding = 0,
                  nudge_x = 150,
                  nudge_y = 1,
                  segment.color = "black",
                  fontface = 'bold') +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank()
  )

print(g)

#write.csv(x = specific_correlations, "ADIPOQ_significant_correlations.csv", row.names = F)
#ggsave(paste0("2A_002_rank_", select_correlation$Source_gene[1], ".png"), plot = g, width = 6, height = 6, dpi = 300)
#ggsave(paste0("2A_002_rank_", select_correlation$Source_gene[1], ".pdf"), plot = g, width = 6, height = 6, dpi = 300)

