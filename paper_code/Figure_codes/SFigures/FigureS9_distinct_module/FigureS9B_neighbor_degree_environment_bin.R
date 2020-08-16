#### Find neighbor proteins and check their degree at each stability bin
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
vScore_PPI = csvReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
PPI_pair = split_string_vector(vScore_PPI[,1])

## Calculate the PPI degree for neighbor proteins after removing the target (center) protein
protein_neighbor_degree = function(PPI_pair, target_protein, neighbor_proteins){
  index_not_target = which(PPI_pair[,1] != target_protein & PPI_pair[,2] != target_protein)
  PPI_pair_select = PPI_pair[index_not_target,]
  protein_degree = as.data.frame(table(as.character(c(PPI_pair_select[,1], PPI_pair_select[,2]))))
  protein_degree_select = protein_degree[which(protein_degree[,1] %in% neighbor_proteins),2]
  return(protein_degree_select)
}

### Find neighbor proteins based on igraph and calculate their stability scores after removing the target proteins
neighbor_protein_degree_bin = function(vScore_PPI, PPI_pair, bin, g){
  stability_1_PPI = vScore_PPI[which(as.numeric(vScore_PPI[,2]) == bin), 1]
  PPI_pair_bin = split_string_vector(stability_1_PPI)
  protein_unique = unique(as.vector(PPI_pair_bin))
  neighbor_degree = 0
  for (i in 1: length(protein_unique)){
    temp = V(g)$name[neighbors(g, protein_unique[i])]
    temp_degree = protein_neighbor_degree(PPI_pair, protein_unique[i], temp)
    neighbor_degree = c(neighbor_degree, temp_degree)
  }
  neighbor_degree = neighbor_degree[2: length(neighbor_degree)]
  return(neighbor_degree)
}
## Input a network and find its neighbours for a specific node
PPI_network = data.frame(protein_A = PPI_pair[,1], protein_B = PPI_pair[,2], 
                         Env = as.numeric(vScore_PPI[,3]), vScore = as.numeric(vScore_PPI[,4]))
library(igraph)
g = graph_from_data_frame(PPI_network, directed = FALSE)
stability_1 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 1, g))
stability_2 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 2, g))
stability_3 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 3, g))
stability_4 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 4, g))
stability_5 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 5, g))
stability_6 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 6, g))
stability_7 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 7, g))
stability_8 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 8, g))
stability_9 = data.frame(degree = neighbor_protein_degree_bin(vScore_PPI, PPI_pair, 9, g))

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
        geom_density(alpha = 0.03)+
        scale_color_manual(name = "Number of positive environments", values = col_chosen)+
        scale_fill_manual(name = "Number of positive environments", values = col_chosen)+
        scale_x_continuous(name = "Neighbor's degree", 
                           limits=c(0, 40),
                           breaks = seq(0,40, by =5),
                           labels = seq(0,40, by= 5)) +
        ylab("Density") +
        guides(fill=guide_legend(ncol=3), col = guide_legend(ncol= 3))+
        theme(legend.key = element_blank(), legend.position = c(0.65,0.92),
              legend.text=element_text(size=10),legend.title=element_text(size=10),
              legend.key.size = unit(0.4, "cm"))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"),
              axis.title.y=element_text(size=10)) + 
        theme(text = element_text(size=10))

#ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/neighbor_degree_environment_bin.pdf", width =4, height =4)
ggsave("Figures/SFigures/SFigure9/FigureS9B_Neighbor_degree_environment_bin.pdf", width =4, height =4)


