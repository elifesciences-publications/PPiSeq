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

## Make a heatmap to show the network density for each pair of GO terms
setwd("~/Dropbox/PPiSeq_02")
## Cellular Compartment
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_CC.csv",
                             sep = "\t", header = T))
rownames(network_density) = colnames(network_density)
# dendrogram
library(ggdendro)

hc = hclust(dist(network_density), method="ward.D")

dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 3) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
                  size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom")
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Dendrogram_CC_GO_pair_cluster.pdf", width = 5, height = 5)
label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),]

## Get network density of 100 random networks and take the mean value
## Color the dots with different colors: high than the mean, lower or equal than the mean
setwd("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/")
## Cellular Compartment
network_density = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_CC.csv",
                                       sep = "\t", header = T))
GO_file_name_front = "random_network_CC_density_"
quantile = c(0.025, 0.5, 0.975)
random_network_quantile = function(network_density, GO_file_name_front, quantile){
        density_CC_random = rep(0, nrow(network_density) * ncol(network_density))
        for (i in 1:100){
                density_CC_name = paste(GO_file_name_front, as.character(i), ".txt", sep = "")
                CC_matrix = as.matrix(read.table(density_CC_name, sep = "\t", header = T))
                vector_density = as.vector(CC_matrix)
                density_CC_random = cbind(density_CC_random, vector_density)
        }
        density_CC_random = density_CC_random[,2:ncol(density_CC_random)]
        #mean_density_CC_random = rowMeans(density_CC_random)
        density_CC_quantile = matrix(0, nrow(density_CC_random), 3)
        for(i in 1:nrow(density_CC_quantile)){
                density_CC_quantile[i,] = quantile(density_CC_random[i,], probs = quantile)
        }  
        return(density_CC_quantile)
}
density_CC_quantile = random_network_quantile(network_density, GO_file_name_front, quantile)

GO = colnames(network_density)
network_density_vector = as.vector(network_density)
random_comp = rep("No", length(network_density_vector))
random_comp[which(as.numeric(network_density_vector) >= as.numeric(density_CC_quantile[,3]))] = "High"
#random_comp[which(network_density_vector <= density_CC_quantile[,1])] = "Low"

Network_density = as.numeric(network_density_vector)
dataf = data.frame(rowv = rep(GO, each = 22),
                   columnv = rep(GO, 22),
                   Network_density = as.numeric(network_density_vector),
                   label = random_comp)
#library(grid)
#text_high <- textGrob("Network density", gp=gpar(fontsize=8, fontface="bold"))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        scale_color_manual(name = "", values = apple_colors[c(7, 10)], breaks = c("NA", "NA", "NA"))+  
        scale_x_discrete(limits = GO_order$label) + 
        scale_y_discrete(limits = GO_order$label) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = seq(0.01, 0.04, by = 0.01), 
                   label =as.character(seq(0.01,0.04, by = 0.01)), 
                   range = c(0, 3))+ ## to tune the size of circles
        theme(legend.position= "top",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = apple_colors[11], angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
#ggsave("~/Desktop/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_CC_PPI_network.pdf", width = 6, height = 5)

###### Make a heatmap for all count
network_all_count = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_CC.csv",
                                       sep = "\t", header = T))

GO = colnames(network_all_count)
network_count_vector = as.vector(network_all_count)

dataC = data.frame(rowv = rep(GO, each = 22),
                   columnv = rep(GO, 22),
                   Network_count = as.numeric(network_count_vector)
)
Network_count = as.numeric(network_count_vector)
count_label = c(expression('3 x 10'^5), expression('8 x 10'^5), expression('1.2 X 10'^6))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_count), col= apple_colors[5], dataC)  + 
        scale_x_discrete(limits = GO_order$label) + 
        scale_y_discrete(limits = GO_order$label) +## color of the corresponding aes
        scale_size(name = "Protein protein pair count", breaks = c(3e5, 8e5, 12e5), 
                   label = count_label, 
                   range = c(0, 4))+
        theme(legend.position ="bottom", legend.key = element_blank(), 
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_all_count_CC_PPI_network.pdf", width = 6, height = 5)

########################################
## Biological process
setwd("~/Dropbox/PPiSeq_02")
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_BP.csv",
                                       sep = "\t", header = T))
rownames(network_density) = colnames(network_density)
GO_BP_index = cbind(1:nrow(network_density), rownames(network_density))
write.csv(GO_BP_index, "~/Desktop/BG_GO_index.csv", quote =F, row.names = F)
# dendrogram
library(ggdendro)
#rownames(network_density) = as.character(1:98)
hc = hclust(dist(network_density), method="ward.D")

dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 10) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
                  size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom")
#ggsave("~/Desktop/Dendrogram_BP_GO_pair_cluster.pdf", width = 6, height = 10)
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Dendrogram_BP_GO_pair_cluster.pdf", width = 6, height = 10)

label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),] # order the BP GO term by the cluster analysis
setwd("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/")
network_density = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_BP.csv",
                                       sep = "\t", header = T))
GO_file_name_front = "random_network_BP_density_"
quantile = c(0.025, 0.5, 0.975)
random_network_quantile = function(network_density, GO_file_name_front, quantile){
        density_CC_random = rep(0, nrow(network_density) * ncol(network_density))
        for (i in 1:100){
                density_CC_name = paste(GO_file_name_front, as.character(i), ".txt", sep = "")
                CC_matrix = as.matrix(read.table(density_CC_name, sep = "\t", header = T))
                vector_density = as.vector(CC_matrix)
                density_CC_random = cbind(density_CC_random, vector_density)
        }
        density_CC_random = density_CC_random[,2:ncol(density_CC_random)]
        #mean_density_CC_random = rowMeans(density_CC_random)
        density_CC_quantile = matrix(0, nrow(density_CC_random), 3)
        for(i in 1:nrow(density_CC_quantile)){
                density_CC_quantile[i,] = quantile(density_CC_random[i,], probs = quantile)
        }  
        return(density_CC_quantile)
}
density_CC_quantile = random_network_quantile(network_density, GO_file_name_front, quantile)

GO = colnames(network_density)
network_density_vector = as.vector(network_density)
random_comp = rep("No", length(network_density_vector))
random_comp[which(network_density_vector >= density_CC_quantile[,3])] = "High"

Network_density = as.numeric(network_density_vector)
dataf = data.frame(rowv = rep(GO, each = 98),
                   columnv = rep(GO, 98),
                   Network_density = as.numeric(network_density_vector),
                   label = random_comp)
#library(grid)
#text_high <- textGrob("Network density", gp=gpar(fontsize=8, fontface="bold"))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        #scale_color_gradient2(low = apple_colors[5], high = apple_colors[7])+
        scale_color_manual(name = "", values = apple_colors[c(7,10)], breaks = c("NA", "NA"))+  
        scale_x_discrete(limits = GO_order$label) + 
        scale_y_discrete(limits = GO_order$label) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = seq(0, 0.16, by = 0.04), 
                   label =as.character(seq(0,0.16, by = 0.04)), 
                   range = c(0, 1.5))+ ## to tune the size of circles
        theme(legend.position= "top",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 6, color = apple_colors[11], angle = 90, hjust =1),
              axis.text.y.left = element_text(size = 6, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_BP_PPI_network.pdf", width = 12, height = 12)

###### Make a heatmap for all count
network_all_count = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_BP.csv",
                                         sep = "\t", header = T))

GO = colnames(network_all_count)
network_count_vector = as.vector(network_all_count)

dataC = data.frame(rowv = rep(GO, each = 98),
                   columnv = rep(GO, 98),
                   Network_count = as.numeric(network_count_vector)
)
Network_count = as.numeric(network_count_vector)
count_label = c(expression('100'), expression('10'^3), expression('10'^4))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_count), col= apple_colors[5], dataC)  + 
        scale_x_discrete(limits = GO_order$label) + 
        scale_y_discrete(limits = GO_order$label) +## color of the corresponding aes
        scale_size(name = "Protein protein pair count", breaks = c(100, 1e3, 1e4), 
                   label = count_label, 
                   range = c(0, 4))+
        theme(legend.position ="bottom", legend.key = element_blank(), 
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 5, color = "black", angle = 90, hjust =1),
              axis.text.y.left = element_text(size = 5, color = "black"), axis.title = element_blank())

ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_all_count_BP_PPI_network.pdf", width = 12, height = 12)

#### Molecular function
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_MF.txt",
                                       sep = "\t", header = T))
dim(network_density) # 41
GO = colnames(network_density)
network_density_vector = as.vector(network_density)
network_density_vector[which(network_density_vector == "2")] = 0
dataf = data.frame(rowv = rep(GO, each = 41),
                   columnv = rep(GO, 41),
                   network_density, 
                   Network_density = as.numeric(network_density_vector))
GO_matrix = cbind(1:41, GO)
colnames(GO_matrix) = c("Index", "GO_term")
csvWriter(GO, "Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/MF_GO_term_index.csv")
library(ggplot2)
ggplot(dataf, aes(y = rowv,
                  x = columnv)) +        ## global aes
        #geom_tile(aes(fill = circlefill)) +         ## to get the rect filled
        geom_point(aes(size =Network_density, color = Network_density), show.legend = FALSE)  +    ## geom_point for circle illusion
        scale_color_gradient(low = "yellow",  
                             high = apple_colors[7])+       ## color of the corresponding aes
        scale_size(range = c(0.1,1.5))+             ## to tune the size of circles
        theme(legend.position ="right", legend.key=element_blank(), legend.text=element_blank()) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        #scale_y_continuous(name = "",
        # limits=c(1,98),
        #breaks=seq(1,98, by =1),
        #labels = seq(1,98, by= 1)) +
        #scale_x_continuous(name = "",
        #limits=c(1,98),
        #breaks=seq(1,98, by =1),
        #labels = seq(1,98, by= 1))+
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 5, color = "black", angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 5, color = "black"), axis.title = element_blank())
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_MF_PPI_network.pdf", width =6, height = 5)

######################################################################
# Check random network-1 functional enrichment
# Cellular compartment
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/random_network_CC_density_1.txt",
                                       sep = "\t", header = T))
dim(network_density)
GO = colnames(network_density)
network_density_vector = as.vector(network_density)
dataf = data.frame(rowv = rep(GO, each = 22),
                   columnv = rep(GO, 22),
                   network_density, 
                   Network_density = as.numeric(network_density_vector))
library(ggplot2)
ggplot(dataf, aes(y = rowv,
                  x = columnv)) +        ## global aes
        #geom_tile(aes(fill = circlefill)) +         ## to get the rect filled
        geom_point(aes(size =Network_density, color = Network_density), show.legend = FALSE)  +    ## geom_point for circle illusion
        scale_color_gradient(low = "yellow",  
                             high = apple_colors[7])+       ## color of the corresponding aes
        scale_size(range = c(0.5,3))+             ## to tune the size of circles
        theme(legend.position ="right", legend.key=element_blank(), legend.text=element_blank()) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust =1),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_CC_random_network_1.pdf", width =6, height = 5)

# Biological process
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/random_network_BP_density_1.txt",
                                       sep = "\t", header = T))
dim(network_density)
GO = colnames(network_density)
network_density_vector = as.vector(network_density)
dataf = data.frame(rowv = rep(1:98, each = 98),
                   columnv = rep(1:98, 98),
                   network_density, 
                   Network_density = as.numeric(network_density_vector))
GO[93] # "transcription.from.RNA.polymerase.I.promoter"
GO_matrix = cbind(1:98, GO)
colnames(GO_matrix) = c("Index", "GO_term")
csvWriter(GO, "Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Random_BP_GO_term_index.csv")
library(ggplot2)
ggplot(dataf, aes(y = rowv,
                  x = columnv)) +        ## global aes
        #geom_tile(aes(fill = circlefill)) +         ## to get the rect filled
        geom_point(aes(size =Network_density, color = Network_density), show.legend = FALSE)  +    ## geom_point for circle illusion
        scale_color_gradient(low = "yellow",  
                             high = apple_colors[7])+       ## color of the corresponding aes
        scale_size(range = c(0.1,0.8))+             ## to tune the size of circles
        theme(legend.position ="right", legend.key=element_blank(), legend.text=element_blank()) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        #scale_y_continuous(name = "",
        # limits=c(1,98),
        #breaks=seq(1,98, by =1),
        #labels = seq(1,98, by= 1)) +
        #scale_x_continuous(name = "",
        #limits=c(1,98),
        #breaks=seq(1,98, by =1),
        #labels = seq(1,98, by= 1))+
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 5, color = "black"),
              axis.text.y.left = element_text(size = 5, color = "black"), axis.title = element_blank())
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_BP_random_network_1.pdf", width =6, height = 5)


########################################
## Molecular function
setwd("~/Dropbox/PPiSeq_02")
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_MF.txt",
                                       sep = "\t", header = T))
network_density[which(network_density == 2)] = 0 # No any possible PPIs will be 2. Here I set it up to be zero
rownames(network_density) = colnames(network_density)
GO_BP_index = cbind(1:nrow(network_density), rownames(network_density))
write.csv(GO_BP_index, "~/Desktop/MF_GO_index.csv", quote =F, row.names = F)
# dendrogram
library(ggdendro)
#rownames(network_density) = as.character(1:98)
hc = hclust(dist(network_density), method="ward.D")

dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 5) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster), 
                  size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom")
#ggsave("~/Desktop/Dendrogram_MF_GO_pair_cluster.pdf", width = 6, height = 6)
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/Dendrogram_MF_GO_pair_cluster.pdf", width = 5, height = 5)

label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),] # order the BP GO term by the cluster analysis
setwd("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/")
network_density = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_MF.txt",
                                       sep = "\t", header = T))
network_density[which(network_density == 2)] = 0 # No any possible PPIs will be 2. Here I set it up to be zero
GO_file_name_front = "random_network_MF_density_"
quantile = c(0.025, 0.5, 0.975)
random_network_quantile = function(network_density, GO_file_name_front, quantile){
        density_CC_random = rep(0, nrow(network_density) * ncol(network_density))
        for (i in 1:100){
                density_CC_name = paste(GO_file_name_front, as.character(i), ".txt", sep = "")
                CC_matrix = as.matrix(read.table(density_CC_name, sep = "\t", header = T))
                vector_density = as.vector(CC_matrix)
                density_CC_random = cbind(density_CC_random, vector_density)
        }
        density_CC_random = density_CC_random[,2:ncol(density_CC_random)]
        #mean_density_CC_random = rowMeans(density_CC_random)
        density_CC_quantile = matrix(0, nrow(density_CC_random), 3)
        for(i in 1:nrow(density_CC_quantile)){
                density_CC_quantile[i,] = quantile(density_CC_random[i,], probs = quantile)
        }  
        return(density_CC_quantile)
}
density_CC_quantile = random_network_quantile(network_density, GO_file_name_front, quantile)

GO = colnames(network_density)
network_density_vector = as.vector(network_density)
random_comp = rep("No", length(network_density_vector))
random_comp[which(network_density_vector >= density_CC_quantile[,3])] = "High"

Network_density = as.numeric(network_density_vector)
dataf = data.frame(rowv = rep(GO, each = 41),
                   columnv = rep(GO, 41),
                   Network_density = as.numeric(network_density_vector),
                   label = random_comp)
#library(grid)
#text_high <- textGrob("Network density", gp=gpar(fontsize=8, fontface="bold"))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_density, color = label), dataf)  + 
        #scale_color_gradient2(low = apple_colors[5], high = apple_colors[7])+
        scale_color_manual(name = "", values = apple_colors[c(7,10)], breaks = c("NA", "NA"))+  
        scale_x_discrete(limits = GO_order$label) + 
        scale_y_discrete(limits = GO_order$label) +## color of the corresponding aes
        scale_size(name = "Network density", breaks = seq(0, 0.12, by = 0.04), 
                   label =as.character(seq(0,0.12, by = 0.04)), 
                   range = c(0, 2))+ ## to tune the size of circles
        theme(legend.position= "top",legend.key = element_blank(),
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        #guides(fill=guide_legend(ncol = 2))+
        #annotation_custom(grob = text_high,xmin= 10,xmax=13,ymin=23.5,ymax=24.5)  +
        #coord_cartesian(clip = 'off')+# This focuses the x-axis on the range of interest
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 6, color = apple_colors[11], angle = 90, hjust =1),
              axis.text.y.left = element_text(size = 6, color = "black"), axis.title = element_blank())
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_MF_PPI_network.pdf", width = 6, height = 6)

###### Make a heatmap for all count
network_all_count = as.matrix(read.table("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/Network_all_count_PPI_MF.txt",
                                         sep = "\t", header = T))

GO = colnames(network_all_count)
network_count_vector = as.vector(network_all_count)

dataC = data.frame(rowv = rep(GO, each = 41),
                   columnv = rep(GO, 41),
                   Network_count = as.numeric(network_count_vector)
)
Network_count = as.numeric(network_count_vector)
count_label = c(expression('10'), expression('10'^3), expression('10'^5))
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_count), col= apple_colors[5], dataC)  + 
        scale_x_discrete(limits = GO_order$label) + 
        scale_y_discrete(limits = GO_order$label) +## color of the corresponding aes
        scale_size(name = "Protein protein pair count", breaks = c(10, 1e3, 1e5), 
                   label = count_label, 
                   range = c(0, 3))+
        theme(legend.position ="bottom", legend.key = element_blank(), 
              legend.text = element_text(size = 8, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_text(size = 5, color = "black", angle = 90, hjust =1),
              axis.text.y.left = element_text(size = 5, color = "black"), axis.title = element_blank())

ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_all_count_MF_PPI_network.pdf", width = 6, height = 6)
