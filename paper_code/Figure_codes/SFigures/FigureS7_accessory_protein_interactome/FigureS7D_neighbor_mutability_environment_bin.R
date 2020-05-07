#### Find neighbour proteins and check their stability score at each stability bin
#### stability_score = Sum(number_of_positve_environments_PPI)/number_of_PPI
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
vScore_PPI = csvReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
PPI_pair = split_string_vector(vScore_PPI[,1])

## Calculate the PPI stability score for neighbor proteins after removing the target (center) protein
protein_neighbor_stability = function(vScore_PPI, PPI_pair, target_protein, neighbor_proteins){
  index_not_target = which(PPI_pair[,1] != target_protein & PPI_pair[,2] != target_protein)
  vScore_select = vScore_PPI[index_not_target,]
  PPI_pair_select = PPI_pair[index_not_target,]
  stability_score = rep(0, length(neighbor_proteins))
  for (i in 1:length(neighbor_proteins)){
    index_protein = unique(c(which(PPI_pair_select[,1] == neighbor_proteins[i]),
                             which(PPI_pair_select[,2] == neighbor_proteins[i])))
    stability_score[i]= mean(as.numeric(vScore_select[index_protein,3]))
  }
  return(stability_score)
}

### Find neighbor proteins based on igraph and calculate their stability scores after removing the target proteins
neighbour_vScore_bin = function(vScore_PPI, PPI_pair, bin, g){
  stability_1_PPI = vScore_PPI[which(as.numeric(vScore_PPI[,2]) == bin), 1]
  PPI_pair_bin = split_string_vector(stability_1_PPI)
  protein_unique = unique(as.vector(PPI_pair_bin))
  neighbor_stability_score = 0
  for (i in 1: length(protein_unique)){
    temp = V(g)$name[neighbors(g, protein_unique[i])]
    temp_score = protein_neighbor_stability(vScore_PPI, PPI_pair, protein_unique[i], temp)
    neighbor_stability_score = c(neighbor_stability_score, temp_score)
  }
  neighbor_sCore = neighbor_stability_score[2: length(neighbor_stability_score)]
  return(neighbor_sCore)
}
## Input a network and find its neighbours for a specific node
PPI_network = data.frame(protein_A = PPI_pair[,1], protein_B = PPI_pair[,2], 
                         Env = as.numeric(vScore_PPI[,3]), vScore = as.numeric(vScore_PPI[,4]))
library(igraph)
g = graph_from_data_frame(PPI_network, directed = FALSE)

stability_1 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 1, g))
stability_2 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 2, g))
stability_3 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 3, g))
stability_4 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 4, g))
stability_5 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 5, g))
stability_6 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 6, g))
stability_7 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 7, g))
stability_8 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 8, g))
stability_9 = data.frame(sCore = neighbour_vScore_bin(vScore_PPI, PPI_pair, 9, g))

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
ggplot(stability_bin_degree, aes(x = sCore, fill = label, col = label))+
        geom_density(alpha = 0.03)+
        scale_color_manual(name = "Number of positive environments", values = col_chosen)+
        scale_fill_manual(name = "Number of positive environments", values = col_chosen)+
        scale_x_continuous(name = "Neighbor's variability score", 
                           limits=c(0, 3),
                           breaks = seq(0,3, by =0.5),
                           labels = seq(0,3, by= 0.5)) +
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

ggsave("Figures/SFigures/SFigure7/FigureS7D_neighbor_variability_score_environment_bin.pdf", width =4, height =4)


