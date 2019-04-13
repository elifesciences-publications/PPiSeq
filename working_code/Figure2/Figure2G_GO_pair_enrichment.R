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
network_density = read.table("Working_data/Positive_PPI_environment/PPI_pair_GO/random_network/random_network_density/random_network_density_1.txt",
                             sep = "\t", header = T)
dim(network_density)
network_density= network_density[,1:98]
network_density = as.vector(network_density)
dataf = data.frame(rowv = rep(1:98, each = 98),
                   columnv = rep(1:98, 98),
                   network_density, 
                   circlesize = as.numeric(network_density),
                    circlefill = as.numeric(network_density))
library(ggplot2)
ggplot(dataf, aes(y = factor(rowv),
                  x = factor(columnv))) +        ## global aes
        #geom_tile(aes(fill = circlefill)) +         ## to get the rect filled
        geom_point(aes(size =circlesize, color = circlefill), )  +    ## geom_point for circle illusion
        scale_color_gradient(low = "yellow",  
                            high = "red")+       ## color of the corresponding aes
        scale_size(range = c(0.1,0.6))+             ## to tune the size of circles
        theme_bw()
ggsave("~/Desktop/heatmap_BP_PPI_network.pdf", width =5, height = 5)