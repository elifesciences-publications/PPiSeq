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

#### Find neighbour proteins and check their variation scores at each stability bin
setwd("~/Dropbox/PPiSeq_02/")
vScore_PPI = csvReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")

##### Calculate the vScore for each protein
PPI_pair = split_string_vector(vScore_PPI[,1])
protein_unique = unique(as.vector(PPI_pair))
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
        index_protein = unique(c(which(PPI_pair[,1] == protein_unique[i]),
        which(PPI_pair[,2] == protein_unique[i])))
        vScore_protein[i,2] = mean(as.numeric(vScore_PPI[index_protein,3]))
}

## Input a network and find its neighbours for a specific node
PPI_network = data.frame(protein_A = PPI_pair[,1], protein_B = PPI_pair[,2], 
                         Env = as.numeric(vScore_PPI[,3]), vScore = as.numeric(vScore_PPI[,4]))
library(igraph)
g = graph_from_data_frame(PPI_network, directed = FALSE)

neighbour_vScore_bin = function(vScore_PPI, bin, g, vScore_protein){
        stability_1_PPI = vScore_PPI[which(as.numeric(vScore_PPI[,2]) == bin), 1]
        PPI_pair = split_string_vector(stability_1_PPI)
        protein_unique = unique(as.vector(PPI_pair))
        adj_protein = "0"
        for (i in 1: length(protein_unique)){
                temp = V(g)$name[neighbors(g, protein_unique[i])]
                adj_protein = c(adj_protein, temp)
        }
        adj_protein = adj_protein[2: length(adj_protein)]
        adj_protein_unique = unique(adj_protein)
        #neighbours = adjacent_vertices(g, protein_unique)
        
        all_degree = as.numeric(vScore_protein[match(adj_protein_unique, as.character(vScore_protein[,1])), 2])
        return(all_degree)
}

stability_1 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 1, g, vScore_protein))
stability_2 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 2, g, vScore_protein))
stability_3 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 3, g, vScore_protein))
stability_4 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 4, g, vScore_protein))
stability_5 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 5, g, vScore_protein))
stability_6 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 6, g, vScore_protein))
stability_7 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 7, g, vScore_protein))
stability_8 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 8, g, vScore_protein))
stability_9 = data.frame(degree = neighbour_vScore_bin(vScore_PPI, 9, g, vScore_protein))

stability_1$label = "1"
stability_2$label = "2"
stability_3$label = "3"
stability_4$label = "4"
stability_5$label = "5"
stability_6$label = "6"
stability_7$label = "7"
stability_8$label = "8"
stability_9$label = "9"

stability_bin_degree = rbind(stability_1, stability_2, stability_3,
                             stability_4, stability_5, stability_6,
                             stability_7, stability_8, stability_9)
col_chosen = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
library(ggplot2)
ggplot(stability_bin_degree, aes(x = degree, fill = label, col = label))+
        geom_density(alpha = 0.3)+
        scale_color_manual(name = "Positive environment number", values = col_chosen)+
        scale_fill_manual(name = "Positive environment number", values = col_chosen)+
        scale_x_continuous(name = "Neighbours' stability score", 
                           limits=c(0, 3),
                           breaks = seq(0,3, by =0.5),
                           labels = seq(0,3, by= 0.5)) +
        ylab("Density") +
        guides(fill=guide_legend(ncol=3), col = guide_legend(ncol= 3))+
        theme(legend.key = element_blank(), legend.position = c(0.75,0.9),
              legend.text=element_text(size=10),legend.title=element_text(size=10),
              legend.key.size = unit(0.4, "cm"))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/Stability_score_environment_bin.pdf", width =4, height =4)




