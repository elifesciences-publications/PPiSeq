### Combine all the data and re-calculate the an adjusted p-value
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
# Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
SD = csvReader_T("Growth_curve_validation_data/False_negative_rate/AUC_SD_MTX_T7.csv")
H2O2 = csvReader_T("Growth_curve_validation_data/False_negative_rate/AUC_H2O2_MTX_T7.csv")
HU = csvReader_T("Growth_curve_validation_data/False_negative_rate/AUC_HU_MTX_T5.csv")
Forskolin = csvReader_T("Growth_curve_validation_data/False_negative_rate/AUC_Forskolin_MTX_T5.csv")
FK506 = csvReader_T("Growth_curve_validation_data/False_negative_rate/AUC_FK506_MTX_T7.csv")
NaCl = csvReader_T("Growth_curve_validation_data/False_negative_rate/AUC_NaCl_MTX_T7.csv")
Raffinose = csvReader_T("Growth_curve_validation_data/False_negative_rate/AUC_Raffinose_MTX_T7.csv")

combine = rbind(SD, H2O2, HU, Forskolin, FK506, NaCl, Raffinose)
p_value = as.numeric(combine[,4])
Q_value = p.adjust(p_value, method = "BH")
PPI_calling = as.numeric(Q_value <= 0.05)

combine_matrix = data.frame(SD = PPI_calling[1:30], H2O2 = PPI_calling[31:60], HU = PPI_calling[61:90],
                            Forskolin = PPI_calling[91:120], FK506 = PPI_calling[121:150], 
                            NaCl = PPI_calling[151:180], Raffinose = PPI_calling[181:210])
rownames(combine_matrix) = SD[,1]
env_count = rowSums(combine_matrix)
combine_matrix_final = data.frame(PPI= SD[,1], env_count = env_count, combine_matrix)
csvWriter(combine_matrix_final, "Growth_curve_validation_data/False_negative_rate/Summary_false_negative_validation_environment_BH.csv")
SD_matrix_final = combine_matrix_final[which(combine_matrix_final$SD != 0),]

average = (sum(SD_matrix_final$env_count) - nrow(SD_matrix_final))/nrow(SD_matrix_final) # 2.5
average # 2.5

library(pheatmap)
combine_matrix_SD = combine_matrix[which(combine_matrix$SD != 0),]
heatmap = pheatmap(combine_matrix_SD, cluster_rows = T, cluster_cols = TRUE, show_rownames=T,
                   labels_col  = c("SD", expression('H'[2]* 'O'[2]), "Hydroxyurea",
                                   "Forskolin", "FK506", "NaCl", "Raffinose"),show_colnames=T)
pdf("Figures/SFigures/SFigure5/FigureS5C_PPI_calling_SD_only_PPIs.pdf", width = 7, height = 6)
grid::grid.newpage()
grid::grid.draw(heatmap$gtable)
dev.off()

## Plot the relative dynamics (FigureS5D)
SD_dynamics = (as.numeric(SD[,2]) - as.numeric(SD[,3]))/as.numeric(SD[,3])
H2O2_dynamics = (as.numeric(H2O2[,2]) - as.numeric(H2O2[,3]))/as.numeric(H2O2[,3])
HU_dynamics = (as.numeric(HU[,2]) - as.numeric(HU[,3]))/as.numeric(HU[,3])
Forskolin_dynamics = (as.numeric(Forskolin[,2]) - as.numeric(Forskolin[,3]))/as.numeric(Forskolin[,3])
FK506_dynamics = (as.numeric(FK506[,2]) - as.numeric(FK506[,3]))/as.numeric(FK506[,3])
NaCl_dynamics = (as.numeric(NaCl[,2]) - as.numeric(NaCl[,3]))/as.numeric(NaCl[,3])
Raffinose_dynamics = (as.numeric(Raffinose[,2]) - as.numeric(Raffinose[,3]))/as.numeric(Raffinose[,3])

### Make a density plot
SD_frame = data.frame(label = "SD", dynamics = SD_dynamics[which(combine_matrix$SD != 0)])
H2O2_frame = data.frame(label = "H2O2", dynamics = H2O2_dynamics[which(combine_matrix$SD != 0)])
HU_frame = data.frame(label = "HU", dynamics = HU_dynamics[which(combine_matrix$SD != 0)])
Forskolin_frame = data.frame(label = "Forskolin", dynamics = Forskolin_dynamics[which(combine_matrix$SD != 0)])
FK506_frame = data.frame(label = "FK506", dynamics = FK506_dynamics[which(combine_matrix$SD != 0)])
NaCl_frame = data.frame(label = "NaCl", dynamics = NaCl_dynamics[which(combine_matrix$SD != 0)])
Raffinose_frame = data.frame(label = "Raffinose", dynamics = Raffinose_dynamics[which(combine_matrix$SD != 0)])

density_matrix = rbind(SD_frame, H2O2_frame, HU_frame, Forskolin_frame, FK506_frame,
                       NaCl_frame, Raffinose_frame)
#density_matrix = rbind(SD_frame, H2O2_frame,  FK506_frame,Raffinose_frame)
library(tidyverse)
means <- group_by(density_matrix, label) %>% 
  summarise(mean = mean(dynamics))

col_chosen = apple_colors[c(7,6,8,11,5,10,3)]
library(ggplot2)
ggplot(density_matrix, aes(x = dynamics, fill = label, col = label))+
  geom_density(alpha = 0.03)+
  geom_vline(aes(xintercept = mean, color = label), data = means, linetype = 2) +
  
  scale_color_manual(name = "Environment", values = col_chosen)+
  scale_fill_manual(name = "Environment", values = col_chosen)+
  scale_x_continuous(name = "Relative increased AUC", 
                     limits=c(-0.1, 0.5),
                     breaks = seq(-0.1,0.5, by =0.1),
                     labels = seq(-0.1,0.5, by= 0.1)) +
  ylab("Density") +
  guides(fill=guide_legend(ncol=2), col = guide_legend(ncol= 2))+
  theme(legend.key = element_blank(), legend.position = c(0.76,0.8),
        legend.text=element_text(size=10),legend.title=element_text(size=10),
        legend.key.size = unit(0.4, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"),
        axis.title.y=element_text(size=10)) + 
  theme(text = element_text(size=10))

ggsave("Figures/SFigures/SFigure5/FigureS5D_Relative_AUC_SD_only_PPIs.pdf", width =4, height =4)

