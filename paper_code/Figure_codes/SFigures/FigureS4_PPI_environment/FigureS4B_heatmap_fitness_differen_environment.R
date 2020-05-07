#### Make a heatmap to show different groups of PPIs
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
PPI_heatmap = dataFrameReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
heatmap_matrix = PPI_heatmap[,4:12]
fitness_all = unique(as.vector(as.matrix(heatmap_matrix)))
bk2 = seq(0, 1, by = 0.01)
bk3 = seq(1.05, 1.6, by = 0.05)
col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])

my_palette = c(colorRampPalette(col_chosen)(length(bk2)),
               rep(apple_colors[7], length(bk3)))
fit_heatmap = pheatmap(heatmap_matrix, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames=FALSE,
                       show_colnames=T, col = my_palette, 
                       labels_col = c("SD",expression('H'[2]* 'O'[2]), "Hydroxyurea", "Doxorubicin",
                                      "Forskolin", "Raffinose", "NaCl", "16 \u00B0C", "FK506"),
                       breaks = c(bk2, bk3), treeheight_row = 0, angle_col = 45)

save_pheatmap_pdf <- function(x, filename, width=4.5, height=5) {
        pdf(filename, width = width, height = height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
}

save_pheatmap_pdf(fit_heatmap, "Figures/SFigures/SFigure4/FigureS4B_fitness_environment_heatmap.pdf")
