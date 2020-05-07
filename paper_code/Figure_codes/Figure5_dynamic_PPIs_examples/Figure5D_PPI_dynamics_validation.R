############### Make a scatter plot to show the comparison between Tecan and PPiSeq
##### Separate each group by HXT1, HXT3, HXT5, and HXT7
##### Add error bars onto each spot (x-axis, y-axis)
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
SD_Raff_final = dataFrameReader_T("Figure5_related_data/SD_Raff_Tecan_PPiSeq_all_comparison.csv")
name_exchange = csvReader_T("Outsourced_datasets/Systematic_standard_protein.csv")
PPI_split = split_string_vector(SD_Raff_final[,1])
protein_1 = name_exchange[match(PPI_split[,1], name_exchange[,1]),2]
protein_2 = name_exchange[match(PPI_split[,2], name_exchange[,1]),2]
label = rep("0", nrow(SD_Raff_final))
for(i in 1:length(protein_1)){
        if(protein_1[i] == "HXT1" | protein_2[i] == "HXT1"){
                label[i] = "HXT1"
        }
        else if(protein_1[i] == "HXT5" | protein_2[i] == "HXT5"){
                label[i] = "HXT5"
        }
        else if(protein_1[i] == "HXT3" | protein_2[i] == "HXT3"){
                label[i] = "HXT3"
        }
        else if(protein_1[i] == "HXT7" | protein_2[i] == "HXT7"){
                label[i] = "HXT7"
        }else if(protein_1[i] == "HXT2" | protein_2[i] == "HXT2"){
                label[i] = "HXT2"
        }else{
                label[i] = "Other"
        }
}
SD_Raff_final = data.frame(protein_1, protein_2, SD_Raff_final, label)
SD_Raff_final$p.value.1 = P_value
PPiseq_min = SD_Raff_final$Mean_diff - SD_Raff_final$SD_diff
PPiseq_max = SD_Raff_final$Mean_diff + SD_Raff_final$SD_diff
Tecan_min = SD_Raff_final$Mean_dynamics - SD_Raff_final$SD_dynamics
Tecan_max = SD_Raff_final$Mean_dynamics + SD_Raff_final$SD_dynamics
SD_Raff_final = data.frame(SD_Raff_final, PPiseq_min, PPiseq_max, Tecan_min, Tecan_max)
cor(SD_Raff_final$Mean_dynamics, SD_Raff_final$Mean_diff, method = "spearman") # 0.61
SD_Raff_final$label = factor(SD_Raff_final$label, levels = c("HXT1", "HXT7", "HXT3", "HXT5", "HXT2", "Other"))
SD_Raff_final$p.value.1 = factor(SD_Raff_final$p.value.1, levels = c("Significant", "Non-significant", "NA"))
library(ggplot2)
#"HXT1" = "#7b3294", "HXT3" = "#c2a5cf", "HXT5" = "#d01c8b", "HXT7" = "#a6dba0", "Others" ="#008837"
ggplot(data = SD_Raff_final, aes(x = Mean_diff, y = Mean_dynamics))+ 
        geom_point(aes(col = label), size =3, alpha = 0.8)+
        #geom_errorbarh(aes(xmin = PPiseq_min, xmax = PPiseq_max), col = apple_colors[8], size = 0.2 )+
        #geom_errorbar(aes(ymin = Tecan_min, ymax = Tecan_max), col = apple_colors[8], size = 0.2)+
        geom_vline(xintercept = 0, col = apple_colors[11], linetype = 2,size = 0.2)+
        geom_hline(yintercept = 0, col = apple_colors[11], linetype = 2, size = 0.2)+
        annotate("text", x = -0.26, y = 0.6, label = expression(paste("Spearman's ", italic(r), " = 0.61")),  
                 parse = TRUE, col = apple_colors[11], size = 4) +
        scale_color_manual(name = "", values = c("#1b9e77","#e7298a", "#d95f02", "#7570b3", "#1f78b4",  "#CECED2"))+
        scale_y_continuous(name = "Fitness change in Raffinose by OD595",
                           limits=c(-0.5, 0.7),
                           breaks=seq(-0.5,0.7, by =0.2),
                           labels = seq(-0.5,0.7, by= 0.2)) +
        scale_x_continuous(name = "Fitness change in Rraffinose by PPiSeq", 
                           limits=c(-0.5, 0.3),
                           breaks=seq(-0.5,0.3, by =0.1),
                           labels = seq(-0.5,0.3, by= 0.1))+
        theme(legend.key=element_blank(), legend.text=element_text(size=10)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("Figures/Figure5/Figure5D_Raffinose_SD_Tecan_PPiSeq_comparison.pdf", width =4, height =3 )       

################ NaCl Environment
SD_NaCl_final = dataFrameReader_T("Figure5_related_data/SD_NaCl_Tecan_PPiSeq_all_comparison.csv")
PPI_split = split_string_vector(SD_NaCl_final[,1])
protein_1 = name_exchange[match(PPI_split[,1], name_exchange[,1]),2]
protein_2 = name_exchange[match(PPI_split[,2], name_exchange[,1]),2]
label = rep("0", nrow(SD_Raff_final))
for(i in 1:length(protein_1)){
        if(protein_1[i] == "HXT1" | protein_2[i] == "HXT1"){
                label[i] = "HXT1"
        }
        else if(protein_1[i] == "HXT5" | protein_2[i] == "HXT5"){
                label[i] = "HXT5"
        }
        else if(protein_1[i] == "HXT3" | protein_2[i] == "HXT3"){
                label[i] = "HXT3"
        }
        else if(protein_1[i] == "HXT7" | protein_2[i] == "HXT7"){
                label[i] = "HXT7"
        }else if(protein_1[i] == "HXT2" | protein_2[i] == "HXT2"){
                label[i] = "HXT2"
        }else{
                label[i] = "Other"
        }
}
SD_NaCl_final = data.frame(protein_1, protein_2, SD_NaCl_final, label)

PPiseq_min = SD_NaCl_final$Mean_diff - SD_NaCl_final$SD_diff
PPiseq_max = SD_NaCl_final$Mean_diff + SD_NaCl_final$SD_diff
Tecan_min = SD_NaCl_final$Mean_dynamics - SD_NaCl_final$SD_dynamics
Tecan_max = SD_NaCl_final$Mean_dynamics + SD_NaCl_final$SD_dynamics
SD_NaCl_final = data.frame(SD_NaCl_final, PPiseq_min, PPiseq_max, Tecan_min, Tecan_max)
cor(SD_NaCl_final$Mean_dynamics, SD_NaCl_final$Mean_diff, method = "spearman") # 0.21
SD_NaCl_final$label = factor(SD_NaCl_final$label, levels = c("HXT1", "HXT7", "HXT3", "HXT5", "HXT2",  "Other"))
library(ggplot2)
#"HXT1" = "#7b3294", "HXT3" = "#c2a5cf", "HXT5" = "#d01c8b", "HXT7" = "#a6dba0", "Others" ="#008837"
ggplot(data = SD_NaCl_final, aes(x = Mean_diff, y = Mean_dynamics))+ 
        geom_point(aes(col = label), size = 3, alpha = 0.8)+
        #geom_errorbarh(aes(xmin = PPiseq_min, xmax = PPiseq_max), col = apple_colors[8], size = 0.2 )+
        #geom_errorbar(aes(ymin = Tecan_min, ymax = Tecan_max), col = apple_colors[8], size = 0.2)+
        geom_vline(xintercept = 0, col = apple_colors[11], linetype = 2,size = 0.2)+
        geom_hline(yintercept = 0, col = apple_colors[11], linetype = 2, size = 0.2)+
        annotate("text", x = -0.24, y = 0.2, label = expression(paste("Spearman's ", italic(r), " = 0.21")),  
                 parse = TRUE, col = apple_colors[11], size = 3) +
        scale_color_manual(name = "", values = c("#1b9e77","#e7298a", "#d95f02", "#7570b3", "#1f78b4", "#CECED2"))+
        scale_y_continuous(name = "Fitness change in NaCl by OD595",
                           limits=c(-0.4, 0.2),
                           breaks=seq(-0.4,0.2, by =0.1),
                           labels = seq(-0.4,0.2, by= 0.1)) +
        scale_x_continuous(name = "Fitness change in NaCl by PPiSeq", 
                           limits=c(-0.4, 0.3),
                           breaks=seq(-0.4,0.3, by =0.1),
                           labels = seq(-0.4,0.3, by= 0.1))+
        theme(legend.key=element_blank(), legend.text=element_text(size=10)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("Figures/Figure5/Figure5D_NaCl_SD_Tecan_PPiSeq_comparison.pdf", width =4, height =3 )       

