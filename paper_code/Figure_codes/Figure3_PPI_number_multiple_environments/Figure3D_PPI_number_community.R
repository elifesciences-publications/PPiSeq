#### Detect communities and check their distribution of mean variation score
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
vScore_PPI = csvReader_T("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
community = csvReader_T("Figure3_related_network_data/community_infomap.csv")
community[,2] = as.numeric(community[,2])

stable = community[which(as.numeric(community[,2]) == 1),] # 312 
middle = community[which(as.numeric(community[,2]) == 3),] # 146
unstable = community[which(as.numeric(community[,2]) == 2),] # 490
combine = rbind(stable, middle, unstable) # 948
others = community[which(!community[,1] %in% combine[,1]),] # 1134
others[,2] = "4"

### Replace each protein with module number
all = rbind(stable, middle, unstable, others)
vScore_protein = split_string_vector(vScore_PPI[,1])
module = matrix(0,nrow(vScore_protein),2)
module[,1] = all[match(vScore_protein[,1], all[,1]),2]
module[,2] = all[match(vScore_protein[,2], all[,1]),2]

### Replace each PPI with module connection
connection = paste(module[,1], module[,2], sep = "_")
### Split the module connection into stable and unstable
connection_unstable = connection[which(as.numeric(vScore_PPI[,2]) <= 4)] # 10299
connection_stable = connection[which(as.numeric(vScore_PPI[,2]) >4)] # 2682

### write a function to calculate number of PPIs within each module or between modules
count_connection = function(connection){
  count_matrix = as.matrix(data.frame(table(connection)))
  unique = mark_duplicates_fast(count_matrix[,1])
  unique_count = rep(0, nrow(unique))
  for(i in 1:nrow(unique)){
    if(unique[i,2] != 0){
      a = as.numeric(count_matrix[which(count_matrix[,1] == unique[i,1]),2])
      b = as.numeric(count_matrix[which(count_matrix[,1] == unique[i,2]),2])
      unique_count[i] = a + b 
    }else{
      a = as.numeric(count_matrix[which(count_matrix[,1] == unique[i,1]),2])
      b = 0
      unique_count[i] = a + b 
    }
    
  }
  connection_count = cbind(unique[,1], unique_count)
  connection_count = connection_count[order(connection_count[,1]),]
  return(connection_count)
}

### count stable and unstable PPIs within each module or between modules
connection_unstable_count = count_connection(connection_unstable)
connection_stable_count = count_connection(connection_stable)

connection_matrix = cbind(connection_unstable_count, as.numeric(connection_stable_count[,2]))
colnames(connection_matrix) = c("Module", "Mutable", "Immutable")
csvWriter(connection_matrix, "Figure3_related_network_data/Figure3D_PPI_number_modules.csv")
### 1: Community_core_low_mutability; 2: Community_accessory_high_mutability; 3: Community_core_intermediate_mutability

### make pie chart for each module
connection_matrix_pie = connection_matrix[c(1,5,8,10),]
unstable_count = as.numeric(connection_matrix_pie[,2])
stable_count = as.numeric(connection_matrix_pie[,3])
pdf("Figures/Figure3/Figure3D_PPI_number_module_core_low_mutability.pdf", 
    width =5, height = 5)
pie(c(unstable_count[1], stable_count[1]), col= apple_colors[c(5,7)], 
    labels = NA, main = NA)
dev.off()
pdf("Figures/Figure3/Figure3D_PPI_number_module_accessory_high_mutability.pdf", 
    width =5, height = 5)
pie(c(unstable_count[2], stable_count[2]), col= apple_colors[c(5,7)], 
    labels= NA, main = NA)
dev.off()
pdf("Figures/Figure3/Figure3D_PPI_number_module_core_intermediate_mutability.pdf", 
    width =5, height = 5)
pie(c(unstable_count[3], stable_count[3]), col= apple_colors[c(5,7)], 
    labels= NA, main = NA)
dev.off()

