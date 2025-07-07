library(tidyverse)

summary_file <- read.csv('/Volumes/Calvin/TissueCrosstalk/Figure1 items/correlation_count_summary.csv')

organized_file <- summary_file %>%
  group_by(Source) %>%
  mutate(IQR = IQR(Sig_correlations)) %>%
  ungroup() %>%
  mutate(Source = reorder(Source, -IQR))

main_plot <- ggplot(organized_file, aes(x=Source, y=Sig_correlations)) + 
  geom_boxplot(fill = "lightgreen") +
  theme(
    plot.title = element_text(size = 12, hjust = 0.5, face = 'bold'),
    axis.title.x = element_text(size = 10, face = 'bold'),
    axis.title.y = element_text(size = 10, face = 'bold'),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = 'black'),
    axis.text.y = element_text(size = 9, angle = 45, hjust = 1, color = 'black'),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")) +
  labs(title="Number of Signficant Correlations per Tissue", y = "Significant Correlations", x = "Source Tissues")

inset_data <- organized_file %>%
  filter(Source %in% levels(organized_file$Source)[50:52])
y_min <- min(inset_data$Sig_correlations, na.rm = TRUE)
y_max <- max(inset_data$Sig_correlations, na.rm = TRUE)

inset_plot <- ggplot(inset_data,
                     aes(x=Source, y=Sig_correlations)) +
  geom_boxplot(fill = "lightblue") +
  coord_cartesian(ylim = c(0, NA)) +
  theme_minimal() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = 'black'),
        axis.text.y = element_text(size = 9, color = 'black'),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black")
#        axis.text.x = element_blank())
)

combined_plot <- main_plot +
  annotation_custom(ggplotGrob(inset_plot),
                    xmin = 35, xmax = 45, ymin = 25000, ymax = 90000)

print(combined_plot)

#ggsave('/Volumes/Calvin/TissueCrosstalk/F1B_BoxPlot_001.pdf', combined_plot, width = 12, height = 8, dpi = 300)
#ggsave('/Volumes/Calvin/TissueCrosstalk/F1B_BoxPlot_001.png', combined_plot, width = 12, height = 8, dpi = 300)


