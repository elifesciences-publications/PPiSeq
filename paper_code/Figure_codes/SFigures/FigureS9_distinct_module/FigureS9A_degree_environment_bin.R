setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

vScore_PPI = csvReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")

protein_degree_count = function(PPI){
        all_PPI_gene = split_string_vector(PPI)
        protein_degree = as.data.frame(table(as.character(c(all_PPI_gene[,1], all_PPI_gene[,2]))))
        protein_degree_order= protein_degree[order(protein_degree[,2], decreasing = T),]
        return(protein_degree_order)
}

degree_protein = protein_degree_count(vScore_PPI[,1])

#### I put all the protein degree data into each bin

degree_bin = function(vScore_PPI, bin, degree_protein){
        stability_1_PPI = vScore_PPI[which(as.numeric(vScore_PPI[,2]) == bin), 1]
        PPI_pair = split_string_vector(stability_1_PPI)
        protein_unique = unique(as.vector(PPI_pair))
        all_degree = degree_protein[match(protein_unique, as.character(degree_protein[,1])), 2]
        return(all_degree)
}

stability_1 = data.frame(degree = degree_bin(vScore_PPI, 1, degree_protein))
stability_2 = data.frame(degree = degree_bin(vScore_PPI, 2, degree_protein))
stability_3 = data.frame(degree = degree_bin(vScore_PPI, 3, degree_protein))
stability_4 = data.frame(degree = degree_bin(vScore_PPI, 4, degree_protein))
stability_5 = data.frame(degree = degree_bin(vScore_PPI, 5, degree_protein))
stability_6 = data.frame(degree = degree_bin(vScore_PPI, 6, degree_protein))
stability_7 = data.frame(degree = degree_bin(vScore_PPI, 7, degree_protein))
stability_8 = data.frame(degree = degree_bin(vScore_PPI, 8, degree_protein))
stability_9 = data.frame(degree = degree_bin(vScore_PPI, 9, degree_protein))

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
        geom_density(alpha = 0.05)+
        scale_color_manual(name = "Number of positive environments", values = col_chosen)+
        scale_fill_manual(name = "Number of positive environments", values = col_chosen)+
        scale_x_continuous(name = "Degree", 
                           limits=c(0, 120),
                           breaks = seq(0,120, by =20),
                           labels = seq(0,120, by= 20)) +
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

#ggsave("Working_figure/Figure3_accessory_PPIs/accessory_PPI/degree_environment_bin.pdf", width =4, height =4)

ggsave("Figures/SFigures/SFigure9/FigureS9A_Degree_environment_bin.pdf", width =4, height =4)




