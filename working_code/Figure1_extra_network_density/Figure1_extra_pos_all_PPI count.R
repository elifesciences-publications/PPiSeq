###########################
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

##### Plot the number of protein pairs across each cellular compartment or biological processes
setwd("~/Dropbox/PPiSeq_02")
# Cellular compartment all protein protein pair count
transform_list_matrix = function(all_count_name, output_all_count_name){
        all_count= as.matrix(read.table(all_count_name,sep = "\t", header = T))
        GO_unique = unique(c(all_count[,1], all_count[,2]))
        matrix_all_count = matrix(0, length(GO_unique), length(GO_unique))
        for(i in 1:length(GO_unique)){
                a = GO_unique[i]
                for(j in 1:length(GO_unique)){
                        b = GO_unique[j]
                        if(length(which(all_count[,1] == a & all_count[,2] == b)) != 0){
                                matrix_all_count[i,j] = as.numeric(all_count[which(all_count[,1] == a & all_count[,2] == b),3])
                                
                        }
                        else{
                                matrix_all_count[i,j] = 0
                        }
                        
                }
        }
        colnames(matrix_all_count) = GO_unique
        write.table(matrix_all_count, output_all_count_name, sep = "\t", row.names = F, quote=F)
}

#Transform a list to a matrix to make a heatmap
#Cellular compartment
all_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_CC_new.txt"
output_all_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_CC_heatmap.txt"
transform_list_matrix(all_count_name, output_all_count_name)
#Biological process
all_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/BP_primary/Network_all_count_PPI_BP_new.txt"
output_all_count_name = "Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_BP_heatmap.txt"
transform_list_matrix(all_count_name, output_all_count_name)

### Cellular Compartment all protein protein pair count
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_CC_heatmap.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column

GO_CC = colnames(network_density)
rowv = rep(GO_CC, each = 22)
columnv = rep(GO_CC, 22)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network


GO_order = read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/environment/GO_CC_order.txt",
                      header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) # 883490
dataf = data.frame(rowv,columnv,Network_density)
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density), color = apple_colors[5], dataf)  + 
        
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Count", breaks = c(20,200, 20000, 200000, 800000), 
                   label = c("20","200", "20000", "200000", "800000"), range = c(0,4))+ ## to tune the size of circles
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("Working_figure/Figure1_extra/Supplementaryheatmap_count_CC_all_PPI_network.pdf", width = 8, height =4.5)

##### Cellular compartment positive PPI count

network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/environment/DMSO/Network_pos_count_PPI_CC_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column

GO_CC = colnames(network_density)
rowv = rep(GO_CC, each = 22)
columnv = rep(GO_CC, 22)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### This name will be used to match the summary or random network


GO_order = read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/environment/GO_CC_order.txt",
                      header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) # 4112
dataf = data.frame(rowv,columnv,Network_density)
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density), color = apple_colors[5], dataf)  + 
        
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Count", breaks = c(0,500, 1000, 2000, 4000), 
                   label = c("0","500", "1000", "2000", "4000"), range = c(0,4))+ ## to tune the size of circles
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("Working_figure/Figure1_extra/Supplementaryheatmap_count_CC_pos_PPI_network.pdf", width = 8, height =4.5)


### Biological process positive PPI count
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/environment/DMSO/Network_pos_count_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
network_density_vector = as.vector(network_density) # by column

chosen_BP = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/chosen_BP_GO_59.txt",
                                 header = F, sep = "\t"))
all_BP = colnames(network_density)
row_column_good = which(all_BP %in% chosen_BP[,1])
network_density = network_density[row_column_good,row_column_good] # 59

network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 59)
columnv = rep(GO_BP, 59)
GO_GO_name_order = paste(rowv, columnv, sep= "_")### T

GO_order = read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/environment/GO_BP_order_chosen.txt",
                      header = T, sep = "\t")

Network_density = as.numeric(network_density_vector) #216
dataf = data.frame(rowv,columnv,Network_density)
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density), color = apple_colors[5], dataf)  + 
        
        scale_x_discrete(limits = GO_order$x) + 
        scale_y_discrete(limits = GO_order$x) +## color of the corresponding aes
        scale_size(name = "Count", breaks = c(0, 50, 100, 150, 200), 
                   label = c("0","50", "100", "150","200"), range = c(0,4))+ ## to tune the size of circles
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("Working_figure/Figure1_extra/Supplementaryheatmap_count_BP_pos_PPI_network.pdf", width = 12, height =7.5)

