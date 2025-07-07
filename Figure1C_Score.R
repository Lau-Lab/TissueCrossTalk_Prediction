library(tidyverse)
library(ggrepel)

main_directory <- '/Volumes/Calvin/TissueCrosstalk' 

#Key tissue types of interest
conditions_key <- c("Adipose - Subcutaneous",
                    "Heart - Left Ventricle",
                    "Liver",
                    "Adipose - Visceral (Omentum)",
                    "Heart - Atrial Appendage",
                    "Muscle - Skeletal")

#Adipose - Subcutaneous
#Adipose - Visceral (Omentum)
#Heart - Left Ventricle
#Liver
#Muscle - Skeletal


#Adipose to heart and vice versa
#Adipose - Subcutaneous_to_Heart - Left Ventricle
#condition_directory <- file.path(main_directory, 'Tissue_Pairs', 'Adipose - Subcutaneous_to_Heart - Left Ventricle')
#Adipose - Visceral (Omentum)_to_Heart - Left Ventricle
#condition_directory <- file.path(main_directory, 'Tissue_Pairs', 'Adipose - Visceral (Omentum)_to_Heart - Left Ventricle')
#Heart - Left Ventricle_to_Adipose - Subcutaneous
#condition_directory <- file.path(main_directory, 'Tissue_Pairs', 'Heart - Left Ventricle_to_Adipose - Subcutaneous') 
#Heart - Left Ventricle_to_Adipose - Visceral (Omentum)
#condition_directory <- file.path(main_directory, 'Tissue_Pairs', 'Heart - Left Ventricle_to_Adipose - Visceral (Omentum)')

#Modulable
condition_directory <- file.path(main_directory, 
                                 'Tissue_Pairs', 
                                 paste0(conditions_key[2], 
                                        "_to_", 
                                        conditions_key[1]))

#Adipose to heart
highlighted_genes <- c("ADIPOQ","THPO","LGALS3")
gene_colors <- setNames(c("red", "gold", "brown", "brown2"), highlighted_genes)

#Liver to adipose
#highlighted_genes <- c("AHSG", "SHBG", "F11", "INHBE") 
#gene_colors <- setNames(c("peru", "brown", "maroon", "orchid"), highlighted_genes)

#Muscle to liver
#highlighted_genes <- c("ERFE")
#gene_colors <- setNames(c("peru"), highlighted_genes)

##Liver to heart
#highlighted_genes <- c(#"ENHO", 
##                       "FGF21", 
##                       "SELENOP",
#                       "AHSG",
##                       "FETUB",
##                       "AMBP",
##                       "F12",
#                       "F11",
#                       "INHBE"
#)
#
#
##highlighted_genes <- c("LGFBP7", 
##                       "LIPC", 
##                       "EMILIN1", 
##                       "LGALS9", 
##                       "ST6GA1", 
##                       "GHR", 
##                       "CRLF2", 
##                       "LCAT")
##
#gene_colors <- setNames(c("peru", 
#                          "brown", 
#                          "maroon", 
#                          "darkorchid", 
#                          "darksalmon", 
#                          "darkseagreen2", 
#                          "cadetblue4", 
#                          "rosybrown2"),
#                        highlighted_genes)
#
#condition_directory <- file.path(main_directory, 'Tissue_Pairs', 'Heart - Left Ventricle_to_Adipose - Visceral (Omentum)')
#condition_directory <- file.path(main_directory, 'Tissue_Pairs', 'Heart - Left Ventricle_to_Adipose - Visceral (Omentum)')
#condition_directory <- file.path(main_directory, 'Tissue_Pairs', 'Heart - Left Ventricle_to_Adipose - Visceral (Omentum)')
#


Ssec_results <- read.csv(file.path(condition_directory, 'Ssec_results.csv')) %>%
  mutate(Source_gene = gsub("Source_", "", Source_gene))

condition <- trimws(strsplit(basename(condition_directory), "_to_")[[1]])
source_tissue <- condition[1]
target_tissue <- condition[2]


highlighted_df <- Ssec_results %>%
  filter(Source_gene %in% highlighted_genes)

g <- ggplot(Ssec_results, aes(x = p_Ssec_z, y = -r_Ssec_z)) +
  geom_point(size = 1, 
             aes(color = Source_gene),
             alpha = 0.5) +
  labs(title = paste(source_tissue, "to", target_tissue),
       x = "p score",
       y = "r score",
       color = "Source Gene") + 
  xlim(-2, 4) +
  theme(plot.title = element_text(size = 14, h = 0.5, face = 'bold'),
        axis.text.x = element_text(size = 12, color = 'black'),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 12, face = 'bold'),
        axis.line = element_line(color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
#        legend.title = element_text(size = 20, face = 'bold'),
#        legend.text = element_text(size = 24),
        legend.position = "none") +
  geom_point(data = highlighted_df, 
             aes(color = Source_gene), 
             size = 3,
             alpha = 0.7) + 
  scale_color_manual(values = gene_colors) +
  geom_text_repel(data = highlighted_df, 
                  aes(label = Source_gene),
                  size = 5,
                  box.padding = 2.5,
                  point.padding = 0.8,
                  fontface = "bold",
                  segment.color = "black",
                  nudge_x = 0,
#                  nudge_x = 0.5,
                  nudge_y = 0.5)
  
print(g)
#ggsave(file.path(main_directory, paste0("PNG_Figure_", source_tissue,"_to_",target_tissue,"_8.png")), width = 6.5, height = 6.5, g, dpi = 300)
#ggsave(file.path(main_directory, paste0("PDF_Figure_", source_tissue,"_to_",target_tissue,"_8.pdf")), width = 6.5, height = 6.5, g, dpi = 300)

