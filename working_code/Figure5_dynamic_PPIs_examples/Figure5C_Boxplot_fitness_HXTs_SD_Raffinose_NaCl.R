########################### Make a heatmap to show the normalized fitness values for Carbohydrate transport PPIs
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

# Input the normalized fitness values for all PPIs in each environment
setwd("~/Dropbox/PPiSeq_02/")
PPI_fit = csvReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
GO_slim = as.matrix(read.table("Paper_data/Outside_datasets/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
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

# Make a heatmap for each groups

PPI_fit = dataFrameReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
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
PPI_HXT5 = PPI_carbon_fit[grep("YHR096C", PPI_carbon_fit[,1]),] #45
Group = rep("HXT5", nrow(PPI_HXT5))
PPI_HXT5 = data.frame(PPI_HXT5, Group)
PPI_HXT7 = PPI_carbon_fit[grep("YDR342C", PPI_carbon_fit[,1]),] #74
Group = rep("HXT7", nrow(PPI_HXT7))
PPI_HXT7 = data.frame(PPI_HXT7, Group)
PPI_HXT = rbind(PPI_HXT1, PPI_HXT7, PPI_HXT3, PPI_HXT5, PPI_HXT2) 
PPI_HXT= PPI_HXT[!duplicated(PPI_HXT[,1]),] #255
csvWriter(PPI_HXT, "~/Dropbox/PPiSeq_02/Working_data_2/Glucose_transporter/PPI_HXT_environments.csv")
## Only get SD, Raffinose, NaCl

library(tidyverse)
transform_matrix = function(PPI_HXT1){
      PPI_HXT1 = PPI_HXT[,c(1,4,9,10, 13)]
      a =PPI_HXT1 %>% gather("SD", "Raffinose", "NaCl", key = "Environment", value = "Fitness")
      return (a)
}
HXT_bar = transform_matrix(PPI_HXT)

transform_matrix_all = function(PPI_HXT1){
        a =PPI_HXT1 %>% gather("SD","H2O2", "HU","Dox", "Forskolin",
                               "Raffinose", "NaCl", "X16C", "FK506",
                               key = "Environment", value = "Fitness")
        return (a)
}
HXT_bar_all = transform_matrix_all(PPI_HXT)
#HXT_bar$Environment = factor(HXT_bar$Environment, levels = c("SD","Raffinose", "NaCl"))
HXT_bar = transform_matrix(PPI_HXT)
#violin plot to shown the fitness of each group in SD, Raffinose, NaCl
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
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure5_environment_dynamics/Figure5C_barplot_HXTs_three.pdf", width= 4, height = 3.5)


ggplot()+
        #geom_violin(aes(x = Group, y = Fitness, col = Environment), HXT_bar) +
        geom_boxplot(aes(x = Group, y = Fitness, col = Environment), HXT_bar_all) +
        scale_color_manual(values  = apple_colors[c(1:8,11)]) +
        theme(legend.key = element_blank())+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 20, color = "black"),
              axis.title.x = element_blank(),axis.text.y.left = element_text(size = 20, color = "black")) + 
        theme(text = element_text(size=20))
ggsave("~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/FigureS9_carbohydrate_transport/FigureS9F_barplot_HXTs_all_supplementary.pdf", width= 18, height = 5)


