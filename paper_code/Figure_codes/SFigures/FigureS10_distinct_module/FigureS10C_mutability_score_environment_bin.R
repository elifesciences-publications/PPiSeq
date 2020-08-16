setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
vScore_PPI = csvReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
PPI_pair = split_string_vector(vScore_PPI[,1])
protein_unique = unique(as.vector(PPI_pair)) # 2082
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
  index_protein = unique(c(which(PPI_pair[,1] == protein_unique[i]),
                           which(PPI_pair[,2] == protein_unique[i])))
  vScore_protein[i,2] = mean(as.numeric(vScore_PPI[index_protein,3]))
}
#### I put all the protein degree data into each bin

vScore_bin = function(vScore_PPI, bin, vScore_protein){
  stability_1_PPI = vScore_PPI[which(as.numeric(vScore_PPI[,2]) == bin), 1]
  PPI_pair = split_string_vector(stability_1_PPI)
  protein_unique = unique(as.vector(PPI_pair))
  all_vScore = as.numeric(vScore_protein[match(protein_unique, as.character(vScore_protein[,1])), 2])
  return(all_vScore)
}

stability_1 = data.frame(vScore = vScore_bin(vScore_PPI, 1, vScore_protein))
stability_2 = data.frame(vScore = vScore_bin(vScore_PPI, 2, vScore_protein))
stability_3 = data.frame(vScore = vScore_bin(vScore_PPI, 3, vScore_protein))
stability_4 = data.frame(vScore = vScore_bin(vScore_PPI, 4, vScore_protein))
stability_5 = data.frame(vScore = vScore_bin(vScore_PPI, 5, vScore_protein))
stability_6 = data.frame(vScore = vScore_bin(vScore_PPI, 6, vScore_protein))
stability_7 = data.frame(vScore = vScore_bin(vScore_PPI, 7, vScore_protein))
stability_8 = data.frame(vScore = vScore_bin(vScore_PPI, 8, vScore_protein))
stability_9 = data.frame(vScore = vScore_bin(vScore_PPI, 9, vScore_protein))

stability_1$label = "1"
stability_2$label = "2"
stability_3$label = "3"
stability_4$label = "4"
stability_5$label = "5"
stability_6$label = "6"
stability_7$label = "7"
stability_8$label = "8"
stability_9$label = "9"

stability_bin_vScore = rbind(stability_1, stability_2, stability_3,
                             stability_4, stability_5, stability_6,
                             stability_7, stability_8, stability_9)

min(stability_bin_vScore$vScore)
col_chosen = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
library(ggplot2)
ggplot(stability_bin_vScore, aes(x = vScore, fill = label, col = label))+
  geom_density(alpha = 0.05)+
  scale_color_manual(name = "Number of positive environments", values = col_chosen)+
  scale_fill_manual(name = "Number of positive environments", values = col_chosen)+
  scale_x_continuous(name = "Variability score", 
                     limits=c(0, 3),
                     breaks = seq(0,3, by =0.5),
                     labels = seq(0,3, by= 0.5)) +
  ylab("Density") +
  guides(fill=guide_legend(ncol=3), col = guide_legend(ncol= 3))+
  theme(legend.key = element_blank(), legend.position = c(0.6,0.8),
        legend.text=element_text(size=10),legend.title=element_text(size=10),
        legend.key.size = unit(0.5, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"),
        axis.title.y=element_text(size=10)) + 
  theme(text = element_text(size=10))

ggsave("Figures/SFigures/SFigure10/FigureS10C_mutability_score_environment_bin.pdf", width =4, height =4)


