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
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_CC_PPI_network.pdf", width =6, height = 5)

## Biological Process
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/Network_density_PPI_BP.csv",
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
csvWriter(GO, "Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/BP_GO_term_index.csv")
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
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_BP_PPI_network.pdf", width =6, height = 5)

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

## Molecular function
network_density = as.matrix(read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/random_network_MF_density_1.txt",
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
csvWriter(GO, "Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/random_MF_GO_term_index.csv")
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
ggsave("Working_figure/Figure2/Figure2G_PPI_pair_GO_enrichment/heatmap_density_MF_random_network_1.pdf", width =6, height = 5)
