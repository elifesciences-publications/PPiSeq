########################### THis code try to calculate the coefficient variation score for each GO_GO pair
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
### Check the biological process of carbonhydrate transport
setwd("~/Dropbox/PPiSeq_02/working_data/Positive_PPI_environment/PPI_pair_GO/environment")
DMSO_matrix = csvReader_T("DMSO/Network_density_p_value_matrix_BP.csv")
Forskolin_matrix = csvReader_T("Forskolin/Network_density_p_value_matrix_BP.csv")
FK506_matrix = csvReader_T("FK506/Network_density_p_value_matrix_BP.csv")
NaCl_matrix = csvReader_T("NaCl/Network_density_p_value_matrix_BP.csv")
Raffinose_matrix = csvReader_T("Raffinose/Network_density_p_value_matrix_BP.csv")
HU_matrix = csvReader_T("HU/Network_density_p_value_matrix_BP.csv")
H2O2_matrix = csvReader_T("H2O2/Network_density_p_value_matrix_BP.csv")
Dox_matrix = csvReader_T("Dox/Network_density_p_value_matrix_BP.csv")
cold_matrix = csvReader_T("16C/Network_density_p_value_matrix_BP.csv")

variation_matrix = matrix(0, nrow(DMSO_matrix), 3)
variation_matrix[,1:2] = DMSO_matrix[,1:2]
for(i in 1:nrow(DMSO_matrix)){
        density_all = as.numeric(c(DMSO_matrix[i,3],
                             Forskolin_matrix[i,3],
                             FK506_matrix[i,3],
                             NaCl_matrix[i,3],
                             Raffinose_matrix[i,3],
                             HU_matrix[i,3],
                             H2O2_matrix[i,3],
                             Dox_matrix[i,3],
                             cold_matrix[i,3]))
        variation_matrix[i,3] = sd(density_all)/mean(density_all)
}

variation_matrix[which(variation_matrix[,3] == "NaN"),3] = 0
GO_order = read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/environment/GO_BP_order.txt",
                      header = T, sep = "\t")

variation_score = as.numeric(variation_matrix[,3]) # 3
max(variation_score[which(variation_score !=3)]) # 2.679
hist(variation_score, breaks = seq(0,3, by = 0.1))
dataf = data.frame(variation_matrix[,1],variation_matrix[,2], variation_score)
colnames(dataf) = c("rowv", "columnv", "variation")
#library(grid)
#text_high <- textGrob("Network density", gp=gpar(fontsize=8, fontface="bold"))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =variation, color = variation), dataf)  + 
        scale_color_gradientn(name = "P value", colors = apple_colors[c(10,5,7)], 
                              limits = c(0, 2.7),
                              oob = scales::squish) +
        #scale_color_gradientn(name = "P value",colors = apple_colors[c(5,3,7)])+  
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Network density variation", breaks = c(0, 0.5, 1, 0.5, 2, 2.5, 3), 
                   label =as.character(c(0, 0.5, 1, 0.5, 2, 2.5, 3)), 
                   range = c(0, 2.5))+ ## to tune the size of circles
        theme(legend.position= "right",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"),
              axis.text.y.left = element_text(size = 9, color = "black"), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_BP_variation_network.pdf", width = 16, height = 12)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure3/Figure3C_heatmap_density_BP_variation_network.pdf", 
       width = 16, height = 12)



