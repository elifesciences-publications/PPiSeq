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
PPI_heatmap = dataFrameReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment.csv")
install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
PPI_heatmap_order = PPI_heatmap[order(PPI_heatmap$Environment_number, decreasing=F),]
heatmap_matrix = PPI_heatmap_order[,4:12]
colnames(heatmap_matrix) = c("SD", "H2O2", "Hydroxyurea", "Doxorubicin",
                                "Forskolin", "Raffinose", "NaCl", "16 \u00B0C", "FK506")

row_ann = data.frame(Environment = as.character(PPI_heatmap_order$Environment_number))
row.names(row_ann) = rownames(PPI_heatmap_order)
#pdf("~/Desktop/PPI_fitness_across_environment.pdf", height = 5, width =5)
#color_scale=c("#FFFFFF",colorRampPalette((RColorBrewer::brewer.pal(n=7,name="YlGnBu")))(99))
color_scale = c("white", colorRampPalette(brewer.pal(9,"YlGnBu"))(n=99))
fit_heatmap = pheatmap(heatmap_matrix, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames=FALSE,
         annotation_row = row_ann, show_colnames=T, col = color_scale)

save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
        pdf(filename, width = width, height = height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
}

save_pheatmap_pdf(fit_heatmap, "Working_figure/Figure6/Consider_neg_PPI_zero/Figure6B_fitness_environment.pdf")
