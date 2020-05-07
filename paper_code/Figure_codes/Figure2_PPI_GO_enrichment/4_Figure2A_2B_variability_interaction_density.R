##### Check the variation score for each GO-GO pair across different environments
setwd("~/Desktop/PPiSeq_additional_data/Figure2_related_data_generated_by_python_script/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
### Write a function to produce matrix for pos_PPI_count and pos_PPI_density
transform_list_matrix = function(all_count_name, pos_count_name, output_pos_count, output_density){
        all_count= as.matrix(read.table(all_count_name,sep = "\t", header = T))
        pos_count = as.matrix(read.table(pos_count_name,sep = "\t", header = T))
        all_count_name = paste(all_count[,1], all_count[,2], sep = "_")
        pos_count_name = paste(pos_count[,1], pos_count[,2], sep = "_")
        pos_match_all_count = pos_count[match(all_count_name, pos_count_name),3]
        pos_match_all_count[is.na(pos_match_all_count)] = 0
        density_all = as.numeric(pos_match_all_count)/as.numeric(all_count[,3])
        network_density = cbind(all_count, pos_match_all_count, density_all)
        colnames(network_density) = c("GO_1", "GO_2", "all_count", "pos_count", "density")
        GO_unique = unique(network_density[,1])
        
        matrix_pos_count = matrix(0, length(GO_unique), length(GO_unique))
        matrix_density = matrix(0, length(GO_unique), length(GO_unique))
        for(i in 1:length(GO_unique)){
                a = GO_unique[i]
                for(j in 1:length(GO_unique)){
                        b = GO_unique[j]
                        if(length(which(network_density[,1] == a & network_density[,2] == b)) != 0){
                                matrix_pos_count[i,j] = network_density[which(network_density[,1] == a & network_density[,2] == b),4]
                                matrix_density[i,j] = network_density[which(network_density[,1] == a & network_density[,2] == b),5]  
                        }
                        else{
                                matrix_pos_count[i,j] = 0
                                matrix_density[i,j] = 0
                        }
                        
                }
        }
        colnames(matrix_pos_count) = GO_unique
        colnames(matrix_density) = GO_unique
        write.table(matrix_pos_count, output_pos_count, sep = "\t", quote = F)
        write.table(matrix_density, output_density, sep = "\t", quote=F)
}

## SD
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/DMSO/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/DMSO/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/DMSO/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/DMSO/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/DMSO/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/DMSO/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)


## 16C
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/Cold_16C/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/Cold_16C/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/Cold_16C/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/Cold_16C/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/Cold_16C/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/Cold_16C/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

##Dox
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/Doxorubicin/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/Doxorubicin/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/Doxorubicin/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/Doxorubicin/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/Doxorubicin/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/Doxorubicin/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## FK506
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/FK506/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/FK506/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/FK506/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/FK506/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/FK506/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/FK506/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## Forskolin
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/Forskolin/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/Forskolin/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/Forskolin/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/Forskolin/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/Forskolin/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/Forskolin/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## H2O2
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/H2O2/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/H2O2/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/H2O2/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/H2O2/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/H2O2/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/H2O2/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## HU
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/Hydroxyurea/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/Hydroxyurea/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/Hydroxyurea/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/Hydroxyurea/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/Hydroxyurea/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/Hydroxyurea/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## NaCl
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/NaCl/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/NaCl/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/NaCl/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/NaCl/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/NaCl/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/NaCl/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## Raffinose
# Cellular compartment
all_count_name = "Network_all_count_PPI_CC_new.txt"
pos_count_name = "environment/Raffinose/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "environment/Raffinose/Network_pos_count_PPI_CC_matrix.txt"
output_density = "environment/Raffinose/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Network_all_count_PPI_BP_new.txt"
pos_count_name = "environment/Raffinose/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "environment/Raffinose/Network_pos_count_PPI_BP_matrix.txt"
output_density = "environment/Raffinose/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)


##############################################################################
## Calculate the variation score of the network density for each GO_GO pair across different environments

# Cellular compartment
sd_CC_density = as.matrix(read.table("environment/DMSO/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
cold_CC_density = as.matrix(read.table("environment/Cold_16C/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
Dox_CC_density = as.matrix(read.table("environment/Doxorubicin/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
FK506_CC_density = as.matrix(read.table("environment/FK506/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
Forskolin_CC_density = as.matrix(read.table("environment/Forskolin/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
H2O2_CC_density = as.matrix(read.table("environment/H2O2/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
HU_CC_density = as.matrix(read.table("environment/Hydroxyurea/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
NaCl_CC_density = as.matrix(read.table("environment/NaCl/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))
Raffinose_CC_density = as.matrix(read.table("environment/Raffinose/Network_density_PPI_CC_matrix.txt", sep = "\t", header = T))

cc_list = list(sd_CC_density,cold_CC_density,Dox_CC_density,FK506_CC_density, Forskolin_CC_density,
               H2O2_CC_density,HU_CC_density,NaCl_CC_density,Raffinose_CC_density)

cc_list_mean = apply(simplify2array(cc_list), 1:2, mean, na.rm= TRUE)
cc_list_sd = apply(simplify2array(cc_list), 1:2, sd, na.rm = TRUE)

cc_list_CV = cc_list_sd/cc_list_mean

colnames(cc_list_CV) = gsub("\\.", " ", colnames(cc_list_CV))
GO_CC = colnames(cc_list_CV)
rowv = rep(GO_CC, each = 22)
columnv = rep(GO_CC, 22)

network_density_vector = as.vector(cc_list_CV) # by column
CV = as.numeric(network_density_vector)
hist(CV)
max(as.numeric(network_density_vector), na.rm = T) # 3
min(as.numeric(network_density_vector), na.rm = T) # 0.1801661
GO_order = read.table("~/Desktop/PPiSeq_additional_data/Figures/Figure2/GO_CC_order.txt",header = T, sep = "\t")

dataf = data.frame(rowv,columnv,
                   Network_variation = as.numeric(network_density_vector))
csvWriter(dataf, "Variation_CC_all_environments.csv") # will be used to make barplot in next script
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_variation, color = Network_variation), dataf)  + 
        scale_color_gradient(low = apple_colors[10], high = apple_colors[5], 
                             breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3))+  
        scale_size(range = c(0,3), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
        guides(color = guide_legend(), size = guide_legend()) +
        scale_x_discrete(limits = rev(GO_order$x)) + 
        scale_y_discrete(limits = rev(GO_order$x)) +## color of the corresponding aes+ ## to tune the size of circles
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.right = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Desktop/PPiSeq_additional_data/Figures/Figure2/heatmap_density_CC_PPI_network_variation.pdf", width = 8, height =4.5)


# Biological process
sd_BP_density = as.matrix(read.table("environment/DMSO/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
cold_BP_density = as.matrix(read.table("environment/Cold_16C/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
Dox_BP_density = as.matrix(read.table("environment/Doxorubicin/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
FK506_BP_density = as.matrix(read.table("environment/FK506/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
Forskolin_BP_density = as.matrix(read.table("environment/Forskolin/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
H2O2_BP_density = as.matrix(read.table("environment/H2O2/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
HU_BP_density = as.matrix(read.table("environment/Hydroxyurea/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
NaCl_BP_density = as.matrix(read.table("environment/NaCl/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))
Raffinose_BP_density = as.matrix(read.table("environment/Raffinose/Network_density_PPI_BP_matrix.txt", sep = "\t", header = T))

BP_list = list(sd_BP_density,cold_BP_density,Dox_BP_density,FK506_BP_density, Forskolin_BP_density,
               H2O2_BP_density,HU_BP_density,NaCl_BP_density,Raffinose_BP_density)

BP_list_mean = apply(simplify2array(BP_list), 1:2, mean, na.rm = TRUE)
BP_list_sd = apply(simplify2array(BP_list), 1:2, sd, na.rm = TRUE)

BP_list_CV = BP_list_sd/BP_list_mean

colnames(BP_list_CV) = gsub("\\.", " ", colnames(BP_list_CV))
chosen_BP = as.matrix(read.table("chosen_BP_GO_59.txt",header = F, sep = "\t"))
all_BP = colnames(BP_list_CV)
row_column_good = which(all_BP %in% chosen_BP[,1])
network_density = BP_list_CV[row_column_good,row_column_good] # 59

network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 59)
columnv = rep(GO_BP, 59)
GO_GO_name_order = paste(rowv, columnv, sep= "_")

max(as.numeric(network_density_vector), na.rm = T) # 3
min(as.numeric(network_density_vector), na.rm = T) # 0

GO_order = read.table("~/Desktop/PPiSeq_additional_data/Figures/Figure2/GO_BP_order_chosen.txt", header = T, sep = "\t")

dataf = data.frame(rowv,columnv, Network_variation = as.numeric(network_density_vector))
csvWriter(dataf, "Variation_BP_all_environments_chosen.csv")
library(ggplot2)
ggplot() + 
        geom_point(aes(x = rowv, y = columnv, size =Network_variation, color = Network_variation), dataf) + 
        scale_color_gradient(low = apple_colors[10], high = apple_colors[5], 
                             breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3))+  
        scale_size(range = c(0,3), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
        guides(color = guide_legend(), size = guide_legend())+
        scale_x_discrete(limits = rev(GO_order$x)) + 
        scale_y_discrete(limits = rev(GO_order$x)) +## color of the corresponding aes
        theme(legend.justification = "left",
              legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
              legend.key = element_blank(),
              legend.text = element_text(size = 9, color = apple_colors[11])) +
        theme(panel.background = element_blank(), axis.ticks=element_blank(),
              panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
        theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
              axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Desktop/PPiSeq_additional_data/Figures/Figure2/heatmap_density_BP_PPI_network_variation.pdf", width = 12, height =7.5)


##### Get the heatmap of variation for all biological processes

colnames(BP_list_CV) = gsub("\\.", " ", colnames(BP_list_CV))

network_density = BP_list_CV
network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 99)
columnv = rep(GO_BP, 99)
GO_GO_name_order = paste(rowv, columnv, sep= "_")

max(as.numeric(network_density_vector), na.rm = T) # 3
min(as.numeric(network_density_vector), na.rm = T) # 0

GO_order = read.table("~/Desktop/PPiSeq_additional_data/Figures/Figure2/GO_BP_order_primary.txt",
                      header = T, sep = "\t")

dataf = data.frame(rowv,columnv,
                   Network_variation= as.numeric(network_density_vector))

library(ggplot2)
ggplot() + 
  geom_point(aes(x = rowv, y = columnv, size =Network_variation, color = Network_variation), dataf)  + 
  scale_color_gradient(low = apple_colors[10], high = apple_colors[5], 
                       breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3))+  
  scale_size(range = c(0,3), breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)) +
  guides(color = guide_legend(), size = guide_legend())+
  scale_x_discrete(limits = rev(GO_order$x)) + 
  scale_y_discrete(limits = rev(GO_order$x)) +## color of the corresponding aes
  theme(legend.justification = "left",
        legend.position= "right", legend.box = 'vertical',legend.box.just = "left",
        legend.key = element_blank(),
        legend.text = element_text(size = 9, color = apple_colors[11])) +
  theme(panel.background = element_blank(), axis.ticks=element_blank(),
        panel.border = element_rect(colour = apple_colors[10], fill = NA, size = 1))+
  theme(axis.text.x = element_blank(), plot.margin = unit(c(1.5,0.2,0.2,2.2), "cm"),
        axis.text.y.left = element_text(size = 8, color = "black"), axis.title = element_blank())
ggsave("~/Desktop/PPiSeq_additional_data/Figures/Figure2/heatmap_density_BP_PPI_network_variation_all.pdf", 
       width = 16, height =12)

