library(tidyverse)
library(reshape2)

z_score_normalization <- function(x){
  mu_x <- mean(x)
  sigma_x <- sd(x)
  return((x - mu_x)/sigma_x)
}

summary_file <- read.csv('/Volumes/Calvin/TissueCrosstalk/Figure1 items/correlation_count_summary.csv') %>%
  mutate(Target = gsub("_secretomics", "", Target),
         Z_score = z_score_normalization(Sig_correlations),
         log_10 = log10(Sig_correlations + 1e-6),
         log_10z = z_score_normalization(log_10)) %>%
  rename(Significant_Correlations = Sig_correlations) %>%
#  arrange(Significant_Correlations)
  arrange(desc(Significant_Correlations))

unique_ID <- summary_file %>%
  distinct(Source)

unique_identifiers <- rev(unique_ID$Source)

n <- length(unique_identifiers)
square_matrix <- matrix(min(summary_file$log_10), nrow = n, ncol = n)
rownames(square_matrix) <- unique_identifiers
colnames(square_matrix) <- unique_identifiers

rownames_equal_colnames <- rownames(square_matrix) %in% colnames(square_matrix)

# Loop through the row names and set the diagonal values
for (i in seq_along(rownames(square_matrix))) {
  if (rownames(square_matrix)[i] %in% colnames(square_matrix)) {
    square_matrix[rownames(square_matrix)[i], rownames(square_matrix)[i]] <- -Inf
  }
}
#

for(i in 1:nrow(summary_file)){
  row_name <- summary_file$Source[i]
  col_name <- summary_file$Target[i]
  count <- summary_file$log_10[i]
  
  if(row_name %in% rownames(square_matrix) && col_name %in% colnames(square_matrix)) {
    print(paste(row_name, col_name, count))
    square_matrix[row_name, col_name] <- count  # Assign the count to the matrix
  } else {
    warning(paste("Row or column name not found:", row_name, col_name))
  }
}

for_heatmap <- melt(square_matrix, value.name = "log10")
hm <- ggplot(for_heatmap, aes(Var1, Var2, fill = log10)) +
  geom_tile() +
  scale_fill_gradient2(low = "white", 
                       mid = "pink", 
                       high = "purple4", 
                       midpoint = 2.5,
                       limits = c(0,5),
                       oob = scales::squish,
                       na.value = rgb(0,0,0, alpha = 0.3)) +
  labs(title = "Heatmap of Significant Correlations", 
      x = "Target Tissues",
      y = "Source Tissues") +
  theme_minimal() + 
  scale_y_discrete(position="right") +
  theme(plot.title = element_text(size = 12, hjust = 0.5, face = 'bold'),
        axis.title.x = element_text(size = 10, face = 'bold', color = 'black'),
        axis.title.y = element_text(size = 10, face = 'bold', color = 'black'),
        axis.text.x = element_text(size = 9, angle = 90, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 9, color = 'black'),
        axis.line = element_line(color = "black",),
        legend.position = "left"
       )

print(hm)
#ggsave('/Volumes/Calvin/TissueCrosstalk/F1B_Heatmap_002.png', plot = hm, width = 10, height = 10, dpi = 300)
#ggsave('/Volumes/Calvin/TissueCrosstalk/F1B_Heatmap_002.pdf', plot = hm, width = 10 , height = 10, dpi = 300)

