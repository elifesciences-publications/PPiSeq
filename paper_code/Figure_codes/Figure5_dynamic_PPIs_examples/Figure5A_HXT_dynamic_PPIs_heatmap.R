# Input the normalized fitness values for all PPIs in each environment
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
PPI_fit = dataFrameReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
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

PPI_carbon = check_specific_protein(PPI_fit, Gene_carbon) 

# Make a heatmap for each groups
HXT = c("YHR094C","YMR011W","YDR345C","YHR096C", "YDR342C") # (HXT1, HXT2, HXT3, HXT5, HXT7)
PPI_carbon_fit = PPI_fit[which(as.character(PPI_fit[,1]) %in% PPI_carbon),] #277

PPI_HXT1 = PPI_carbon_fit[grep("YHR094C", PPI_carbon_fit[,1]),] #66
Group = rep("HXT1", nrow(PPI_HXT1))
PPI_HXT1 = data.frame(PPI_HXT1,Group )
PPI_HXT2 = PPI_carbon_fit[grep("YMR011W", PPI_carbon_fit[,1]),] #4
Group = rep("HXT2", nrow(PPI_HXT2))
PPI_HXT2 = data.frame(PPI_HXT2,Group )
PPI_HXT3 = PPI_carbon_fit[grep("YDR345C", PPI_carbon_fit[,1]),] #67
Group = rep("HXT3", nrow(PPI_HXT3))
PPI_HXT3 = data.frame(PPI_HXT3, Group)
PPI_HXT5 = PPI_carbon_fit[grep("YHR096C", PPI_carbon_fit[,1]),] #38
Group = rep("HXT5", nrow(PPI_HXT5))
PPI_HXT5 = data.frame(PPI_HXT5, Group)
PPI_HXT7 = PPI_carbon_fit[grep("YDR342C", PPI_carbon_fit[,1]),] #76
Group = rep("HXT7", nrow(PPI_HXT7))
PPI_HXT7 = data.frame(PPI_HXT7, Group) 
PPI_HXT = rbind(PPI_HXT1, PPI_HXT7, PPI_HXT3, PPI_HXT5, PPI_HXT2) # 251
#PPI_FPS1 = PPI_carbon_fit[grep("YLL043W", PPI_carbon_fit[,1]),] #74
PPI_HXT= PPI_HXT[!duplicated(PPI_HXT[,1]),] #248
PPI_other = PPI_carbon_fit[which(!as.character(PPI_carbon_fit[,1]) %in% as.character(PPI_HXT[,1])),] #29
Group = rep("Others", nrow(PPI_other))
PPI_other = data.frame(PPI_other, Group)
PPI_carbon_final = rbind(PPI_HXT, PPI_other) 
## add another annotations (New versus reported)
PPI_reported = dataFrameReader_T("Outsourced_datasets/BIOGRID/multiple_validated_PPI.csv")
PPI_carbon_final_reported = match_both_direction(PPI_carbon_final, PPI_reported[,1])
Reported = rep("No", nrow(PPI_carbon_final))
Reported[which(PPI_carbon_final[,1] %in% PPI_carbon_final_reported[,1])] = "Yes" # 43
PPI_carbon_final = data.frame(PPI_carbon_final, Reported)

## add another annotations (Environment number)
PPI_count = dataFrameReader_T("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter.csv")
count = PPI_count[match(PPI_carbon_final[,1], PPI_count[,1]),2]
PPI_carbon_final = data.frame(PPI_carbon_final, count)

# Cluster PPIs with each HXT protein, respectively and then get the order for the final matrix
PPI_HXT1_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT1"),] # 66
PPI_HXT3_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT3"),] # 4
PPI_HXT5_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT5"),] # 67
PPI_HXT7_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT7"),] # 38
PPI_HXT2_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT2"),] # 76
PPI_others_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "Others"),] # 29

cluster_order = function(PPI_HXT1_matrix){
        PPI_fitness = PPI_HXT1_matrix[,4:12]
        hc = hclust(dist(PPI_fitness), method = "complete")
        PPI_HXT1_matrix_order = PPI_HXT1_matrix[hc$order,]
        return(PPI_HXT1_matrix_order)
        
}
PPI_HXT1_matrix_order = cluster_order(PPI_HXT1_matrix)
PPI_HXT3_matrix_order = cluster_order(PPI_HXT3_matrix)
PPI_HXT5_matrix_order = cluster_order(PPI_HXT5_matrix)
PPI_HXT7_matrix_order = cluster_order(PPI_HXT7_matrix)
PPI_HXT2_matrix_order = cluster_order(PPI_HXT2_matrix)
PPI_others_matrix_order = cluster_order(PPI_others_matrix)

PPI_carbon_final_order = rbind(PPI_HXT1_matrix_order, PPI_HXT7_matrix_order, PPI_HXT3_matrix_order, 
                               PPI_HXT5_matrix_order, PPI_HXT2_matrix_order,
                               PPI_others_matrix_order)

csvWriter(PPI_carbon_final_order, "Figure5_related_data/PPI_HXT_environments.csv")

# Take the order the same with Figure2A
PPI_carbon_heatmap = PPI_carbon_final_order[,c(4,8,10,12,7,9,5,11,6)]
colnames(PPI_carbon_heatmap) = c("SD", "Forskolin","NaCl", "FK506", "Doxorubicin", 
                                 "Raffinose",  "H2O2",  "16 \u00B0C", "Hydroxyurea")
rownames(PPI_carbon_heatmap) = as.character(PPI_carbon_final[,1])
row_ann = data.frame(Protein = as.character(PPI_carbon_final_order$Group))
my_colour = list(Protein = c("HXT1" = "#1b9e77", "HXT7" = "#e7298a","HXT3" = "#d95f02", 
                             "HXT5" = "#7570b3", "HXT2" = "#1f78b4", "Others" ="#CECED2"))
row.names(row_ann) = row.names(PPI_carbon_heatmap)

fitness = as.vector(as.matrix(PPI_carbon_heatmap))
max(fitness) 
bk2 = seq(0, 1, by = 0.01)
bk3 = seq(1.05, 1.3, by = 0.05)
col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
my_palette = c(colorRampPalette(col_chosen)(length(bk2)),
               rep(apple_colors[7], length(bk3)))
library(pheatmap)
library(grid)
fit_heatmap = pheatmap(PPI_carbon_heatmap, cluster_rows = FALSE, cluster_cols = T, show_rownames=FALSE,
                       annotation_colors = my_colour, angle_col = 45, col = my_palette,
                       breaks = c(bk2, bk3),
                       labels_col  = c("SD", "Forskolin","NaCl", "FK506", "Doxorubicin", 
                                       "Raffinose",  expression('H'[2]* 'O'[2]),  "16 \u00B0C", "Hydroxyurea"),
                       annotation_row = row_ann, show_colnames=T)
pdf("Figures/Figure5/Figure5A_carbonhydrate_transport_fitness_environment_primary_cluster.pdf", 
    width = 7, height = 6)
grid::grid.newpage()
grid::grid.draw(fit_heatmap$gtable)
dev.off()
