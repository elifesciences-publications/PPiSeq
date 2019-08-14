# Source some basic functions froma function.R in Github repository
source_https <- function(u, unlink.tmp.certs = FALSE) {
        # load package
        require(RCurl)
        # read script lines from website using a security certificate
        if(!file.exists("cacert.pem")) download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.pem")
        script <- getURL(u, followlocation = TRUE, cainfo = "cacert.pem")
        if(unlink.tmp.certs) unlink("cacert.pem")
        
        # parase lines and evealuate in the global environement
        eval(parse(text = script), envir= .GlobalEnv)
}
source_https("https://raw.githubusercontent.com/sashaflevy/PPiSeq/master/working_code/function.R", unlink.tmp.certs = TRUE)

#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
#### Make a heatmap to show different groups of PPIs
setwd("~/Dropbox/PPiSeq_02/")
PPI_heatmap = dataFrameReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
### Combine the fitness from the two
PPI_count = dataFrameReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
PPI_heatmap$Environment_number = PPI_count[match(as.character(PPI_heatmap$PPI), as.character(PPI_count[,1])),2]
PPI_heatmap_order = PPI_heatmap[order(PPI_heatmap$Environment_number, decreasing=F),]
heatmap_matrix = PPI_heatmap_order[,4:13]
#colnames(heatmap_matrix) = c("SD","SD2", "H2O2", "Hydroxyurea", "Doxorubicin",
                             #"Forskolin", "Raffinose", "NaCl", "16 \u00B0C", "FK506")

#row_ann = data.frame(Environment = as.character(PPI_heatmap_order$Environment_number))
#row.names(row_ann) = rownames(PPI_heatmap_order)
#pdf("~/Desktop/PPI_fitness_across_environment.pdf", height = 5, width =5)
#color_scale=c("#FFFFFF",colorRampPalette((RColorBrewer::brewer.pal(n=7,name="YlGnBu")))(99))
col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
#color_scale = colorRampPalette(col_chosen)(n=100)
fitness_all = unique(as.vector(as.matrix(heatmap_matrix)))
min(fitness_all) # -0.6331257
max(fitness_all) # 1.555071
bk1 = seq(-0.7, -0.05, by = 0.05)
bk2 = seq(0, 1, by = 0.01)
bk3 = seq(1.05, 1.6, by = 0.05)
my_palette = c(rep(apple_colors[5], length(bk1)), colorRampPalette(col_chosen)(length(bk2)),
               rep(apple_colors[7], length(bk3)))
fit_heatmap = pheatmap(heatmap_matrix, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames=FALSE,
                       show_colnames=T, col = my_palette, 
                       labels_col = c("SD", "SD2",expression('H'[2]* 'O'[2]), "Hydroxyurea", "Doxorubicin",
                                      "Forskolin", "Raffinose", "NaCl", "16 \u00B0C", "FK506"),
                       breaks = c(bk1, bk2, bk3), treeheight_row = 0)

save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
        pdf(filename, width = width, height = height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
}

save_pheatmap_pdf(fit_heatmap, "Working_figure/SFigures/Figure2_related/Fitness_environment_correlated_Figure2D.pdf")
