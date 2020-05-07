###### Make a heatmap for SD environment
## Make a heatmap to show the network density for each pair of GO terms
setwd("~/Desktop/PPiSeq_additional_data/")
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
# Cellular compartment
all_count_name = "Figure2_related_data_generated_by_python_script/Network_all_count_PPI_CC_new.txt"
pos_count_name = "Figure2_related_data_generated_by_python_script/Network_pos_count_PPI_CC_new.txt"
output_pos_count = "Figure2_related_data_generated_by_python_script/Network_pos_count_PPI_CC_matrix.txt"
output_density = "Figure2_related_data_generated_by_python_script/Network_density_PPI_CC_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "Figure2_related_data_generated_by_python_script/Network_all_count_PPI_BP_new.txt"
pos_count_name = "Figure2_related_data_generated_by_python_script/Network_pos_count_PPI_BP_new.txt"
output_pos_count = "Figure2_related_data_generated_by_python_script/Network_pos_count_PPI_BP_matrix.txt"
output_density = "Figure2_related_data_generated_by_python_script/Network_density_PPI_BP_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

# Molecular function
all_count_name = "Figure2_related_data_generated_by_python_script/Network_all_count_PPI_MF_new.txt"
pos_count_name = "Figure2_related_data_generated_by_python_script/Network_pos_count_PPI_MF_new.txt"
output_pos_count = "Figure2_related_data_generated_by_python_script/Network_pos_count_PPI_MF_matrix.txt"
output_density = "Figure2_related_data_generated_by_python_script/Network_density_PPI_MF_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

#####(1) I first hierarchical cluster the GO term based on their network density (here I use all the PPIs)
# Cellular compartment
network_density = as.matrix(read.table("Figure2_related_data_generated_by_python_script/Network_density_PPI_CC_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
# dendrogram
library(ggplot2)
library(ggdendro)
rownames(network_density) = colnames(network_density)
hc = hclust(dist(network_density))
dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 4) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        #geom_text(data=label(dendr), aes(x, y, label=label, hjust=0), col = apple_colors[11],size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              axis.line.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
ggsave("Figures/Figure2/Dendrogram_CC_GO_pair_cluster.pdf", width = 3, height = 4.2)

label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),]
csvWriter(GO_order$label, "Figures/Figure2/GO_CC_order.txt")

#### Biological process 
network_density = as.matrix(read.table("Figure2_related_data_generated_by_python_script/Network_density_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))

# dendrogram
library(ggdendro)
rownames(network_density) = colnames(network_density)
hc = hclust(dist(network_density))

dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 4) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        #geom_text(data=label(dendr), aes(x, y, label=label, hjust=0), col = apple_colors[11], size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              axis.line.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom",
              plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
ggsave("Figures/Figure2/Dendrogram_BP_GO_pair_cluster_primary.pdf", width = 4, height =12)

label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),]
csvWriter(GO_order$label, "Figures/Figure2/GO_BP_order_primary.txt")

#### Biological process (Chosen)

network_density = as.matrix(read.table("Figure2_related_data_generated_by_python_script/Network_density_PPI_BP_matrix.txt", 
                                       sep = "\t", header = T))
colnames(network_density) = gsub("\\.", " ", colnames(network_density))
chosen_BP = as.matrix(read.table("Figure2_related_data_generated_by_python_script/chosen_BP_GO_59.txt",
                                 header = F, sep = "\t"))
all_BP = colnames(network_density)
row_column_good = which(all_BP %in% chosen_BP[,1])
network_density = network_density[row_column_good,row_column_good] # 59

# dendrogram
library(ggdendro)
rownames(network_density) = colnames(network_density)
hc = hclust(dist(network_density))

dendr = dendro_data(hc, type= "rectangle")
clust = cutree(hc, k = 4) 
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
#### Plot the dendrogram
ggplot() + 
        geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
        #geom_text(data=label(dendr), aes(x, y, label=label, hjust=0), col = apple_colors[11], size=2) +
        coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
        theme(axis.line.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              axis.line.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              panel.background=element_rect(fill="white"),
              panel.grid=element_blank(),
              legend.position = "bottom",
              plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

ggsave("Figures/Figure2/Dendrogram_BP_GO_pair_cluster_chosen.pdf", width = 3, height =7.5)

label_GO = label(dendr)
GO_order = label_GO[order(label_GO$x),]
csvWriter(GO_order$label, "Figures/Figure2/GO_BP_order_chosen.txt")





