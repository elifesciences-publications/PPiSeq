# Input the normalized fitness values for all PPIs in each environment
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
PPI_fit = read.csv("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv", header = T)
PPI_fit$PPI = as.character(PPI_fit$PPI)
GO_slim = as.matrix(read.table("Outsourced_datasets/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
Gene_carbon = unique(GO_slim[which(GO_slim[,6] == "GO:0008643"), 1])

# A function to extract specific PPIs
check_specific_protein = function(PPI, Gene_Carbon){
        PPI_chosen = "0"
        protein_pair = split_string_vector(PPI[,1])
        if(length(Gene_Carbon) > 1){
                for(i in 1:nrow(protein_pair)){
                        if(protein_pair[i,1] %in% Gene_Carbon | protein_pair[i,2] %in% Gene_Carbon){
                                PPI_chosen = c(PPI_chosen, PPI[i,1])
                        }
                }  
        }else {
                for(i in 1:nrow(protein_pair)){
                        if(protein_pair[i,1] == Gene_Carbon | protein_pair[i,2] == Gene_Carbon){
                                PPI_chosen = c(PPI_chosen, PPI[i,1])
                        }
                }  
        }
        
        PPI_chosen = PPI_chosen[2:length(PPI_chosen)]
        return(PPI_chosen)
}

PPI_carbon = check_specific_protein(PPI_fit, Gene_carbon) # 277

PPI_carbon_fit = PPI_fit[which(as.character(PPI_fit[,1]) %in% PPI_carbon),] #277
PPI_carbon_order = sort(PPI_carbon_fit[,1]) #277
PPI_carbon_fit_order = PPI_carbon_fit[match(PPI_carbon_order, PPI_carbon_fit[,1]),]
PPI_HXT1 = PPI_carbon_fit[grep("YHR094C", PPI_carbon_fit[,1]),] #71
Group = rep("HXT1", nrow(PPI_HXT1))
PPI_HXT1 = data.frame(PPI_HXT1,Group )
PPI_HXT2 = PPI_carbon_fit[grep("YMR011W", PPI_carbon_fit[,1]),] #5
Group = rep("HXT2", nrow(PPI_HXT2))
PPI_HXT2 = data.frame(PPI_HXT2,Group )
PPI_HXT3 = PPI_carbon_fit[grep("YDR345C", PPI_carbon_fit[,1]),] #68
Group = rep("HXT3", nrow(PPI_HXT3))
PPI_HXT3 = data.frame(PPI_HXT3, Group)
PPI_HXT5 = PPI_carbon_fit[grep("YHR096C", PPI_carbon_fit[,1]),] 
Group = rep("HXT5", nrow(PPI_HXT5))
PPI_HXT5 = data.frame(PPI_HXT5, Group)
PPI_HXT7 = PPI_carbon_fit[grep("YDR342C", PPI_carbon_fit[,1]),] 
Group = rep("HXT7", nrow(PPI_HXT7))
PPI_HXT7 = data.frame(PPI_HXT7, Group)
PPI_HXT = rbind(PPI_HXT1, PPI_HXT7, PPI_HXT3, PPI_HXT5, PPI_HXT2) 
PPI_HXT= PPI_HXT[!duplicated(PPI_HXT[,1]),] 

library(tidyverse)
transform_matrix = function(PPI_HXT1){
      PPI_HXT1 = PPI_HXT[,c(1,4,9,10,13)]
      a =PPI_HXT1 %>% gather("SD", "Raffinose", "NaCl", key = "Environment", value = "Fitness")
      return (a)
}
HXT_bar = transform_matrix(PPI_HXT)


HXT_bar$Environment = factor(HXT_bar$Environment, levels = c("SD","Raffinose", "NaCl"))
library(ggplot2)
ggplot()+
        #geom_violin(aes(x = Group, y = Fitness, col = Environment), HXT_bar) +
        geom_boxplot(aes(x = Group, y = Fitness, col = Environment), HXT_bar) +
        scale_color_manual(name = "Environment", breaks = c("SD", 'Raffinose', "NaCl"),
                           values  = apple_colors[c(1,2,4)]) +
        #theme(legend.key = element_blank(), legend.position =c(0.9,0.8))+
        theme(legend.key = element_blank())+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.title.x = element_blank(),axis.text.y.left = element_text(size = 10, color = "black")) + 
        theme(text = element_text(size=10))
ggsave("Figures/Figure5/Figure5B_barplot_HXTs_three_conditions.pdf", width= 4, height = 3.5)
