###########################
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

##### Check the variation score for each GO-GO pair across different environments
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/dynamic_all_pairwise/")
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

CC_name = rep("0", 36)
BP_name = rep("0", 36)
CC_count_name = rep("0", 36)
CC_density_name = rep("0", 36)
BP_count_name = rep("0", 36)
BP_density_name = rep("0", 36)
for (i in 1:36){
  CC_name[i] = paste("dynamic_PPI_CC", paste(as.character(i), "txt", sep="."), sep = "_")
  BP_name[i] = paste("dynamic_PPI_BP", paste(as.character(i), "txt", sep="."), sep = "_")
  CC_count_name[i] = paste("dynamic_PPI_density/dynamic_PPI_CC_count", paste(as.character(i), "txt", sep="."), sep = "_")
  CC_density_name[i] = paste("dynamic_PPI_density/dynamic_PPI_CC_density", paste(as.character(i), "txt", sep="."), sep = "_")
  BP_count_name[i] = paste("dynamic_PPI_density/dynamic_PPI_BP_count", paste(as.character(i), "txt", sep="."), sep = "_")
  BP_density_name[i] = paste("dynamic_PPI_density/dynamic_PPI_BP_density", paste(as.character(i), "txt", sep="."), sep = "_")
}
# Cellular compartment
all_count_name_CC = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
for(i in 1:36){
  pos_count_name = CC_name[i]
  output_pos_count = CC_count_name[i]
  output_density = CC_density_name[i]
  transform_list_matrix(all_count_name_CC, pos_count_name, output_pos_count, output_density)
}

# Biological process
all_count_name_BP = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
for(i in 1:36){
  pos_count_name = BP_name[i]
  output_pos_count = BP_count_name[i]
  output_density = BP_density_name[i]
  transform_list_matrix(all_count_name_BP, pos_count_name, output_pos_count, output_density)
}

##############################################################################
## Calculate the mean of the changing PPI density for each GO_GO pair across different environments
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/dynamic_all_pairwise/")

# Cellular compartment
cc_list = list(rep("0", 36))
for (i in 1:36){
  cc_list[[i]] = as.matrix(read.table(CC_density_name[i], sep = "\t", header = T))

}
#temp = rep(0, 36)
#for(i in 1:36){
  #temp[i] = cc_list[[i]][1,1]
#}
#mean(temp) # 0.3032314
cc_list_mean = apply(simplify2array(cc_list), 1:2, median, na.rm= TRUE) ## Take the median instead of mean
#cc_list_mean = apply(simplify2array(cc_list), 1:2, mean, na.rm= TRUE)
#cc_list_sd = apply(simplify2array(cc_list), 1:2, sd, na.rm = TRUE)
#cc_list_CV = cc_list_sd/cc_list_mean

colnames(cc_list_mean) = gsub("\\.", " ", colnames(cc_list_mean))
GO_CC = colnames(cc_list_mean)
rowv = rep(GO_CC, each = 22)
columnv = rep(GO_CC, 22)

network_density_vector = as.vector(cc_list_mean) # by column
mean = as.numeric(network_density_vector)
hist(mean)
max(as.numeric(network_density_vector), na.rm = T) # 1
min(as.numeric(network_density_vector), na.rm = T) # 0
GO_order = read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_CC_order.txt",header = T, sep = "\t")

dataf = data.frame(rowv,columnv,
                   Dynamic_PPI_ratio = as.numeric(network_density_vector))
csvWriter(dataf, "Changing_PPI_CC_density_mean_environments.csv")
library(ggplot2)
ggplot() + 
  geom_point(aes(x = rowv, y = columnv, size =Dynamic_PPI_ratio, color = Dynamic_PPI_ratio), dataf)  + 
  scale_color_gradient(low = apple_colors[10], high = apple_colors[5], 
                       #breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
                       )+  
  scale_size(range = c(0,3), 
             #breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
             breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
             ) +
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
ggsave("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_mean_density_changing_PPI_CC_median.pdf", width = 8, height =4.5)


# Biological process
BP_list = list(rep("0", 36))
for (i in 1:36){
  BP_list[[i]] = as.matrix(read.table(BP_density_name[i], sep = "\t", header = T))
  
}
#BP_list_mean = apply(simplify2array(BP_list), 1:2, mean, na.rm = TRUE)
BP_list_mean = apply(simplify2array(BP_list), 1:2, median, na.rm = TRUE)
#BP_list_sd = apply(simplify2array(BP_list), 1:2, sd, na.rm = TRUE)
#BP_list_CV = BP_list_sd/BP_list_mean

colnames(BP_list_mean) = gsub("\\.", " ", colnames(BP_list_mean))
chosen_BP = as.matrix(read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/chosen_BP_GO_59.txt",
                                 header = F, sep = "\t"))
all_BP = colnames(BP_list_mean)
row_column_good = which(all_BP %in% chosen_BP[,1])
network_density = BP_list_mean[row_column_good,row_column_good] # 59


network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 59)
columnv = rep(GO_BP, 59)
GO_GO_name_order = paste(rowv, columnv, sep= "_")

max(as.numeric(network_density_vector), na.rm = T) # 1
min(as.numeric(network_density_vector), na.rm = T) # 0

GO_order = read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_BP_order_chosen.txt",
                      header = T, sep = "\t")

dataf = data.frame(rowv,columnv,
                   Dynamic_PPI_ratio = as.numeric(network_density_vector))
csvWriter(dataf, "Changing_PPI_BP_density_mean_environments_chosen.csv")
library(ggplot2)
ggplot() + 
  geom_point(aes(x = rowv, y = columnv, size =Dynamic_PPI_ratio, color = Dynamic_PPI_ratio), dataf)  + 
  scale_color_gradient(low = apple_colors[10], high = apple_colors[5], 
                       #breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))+  
  scale_size(range = c(0,3), 
             #breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
             breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
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
ggsave("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_mean_density_changing_PPI_BP_chosen_median.pdf", width = 12, height =7.5)

############################ plot all the biological functions
colnames(BP_list_mean) = gsub("\\.", " ", colnames(BP_list_mean))

network_density = BP_list_mean # 99

network_density_vector = as.vector(network_density) # by column
GO_BP = colnames(network_density)
rowv = rep(GO_BP, each = 99)
columnv = rep(GO_BP, 99)
GO_GO_name_order = paste(rowv, columnv, sep= "_")

max(as.numeric(network_density_vector), na.rm = T) # 0.03358209
min(as.numeric(network_density_vector), na.rm = T) # 0

GO_order = read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_BP_order_primary.txt",
                      header = T, sep = "\t")

dataf = data.frame(rowv,columnv,
                   Dynamic_PPI_ratio = as.numeric(network_density_vector))
csvWriter(dataf, "Changing_PPI_BP_density_mean_environments_all.csv")
library(ggplot2)
ggplot() + 
  geom_point(aes(x = rowv, y = columnv, size =Dynamic_PPI_ratio, color = Dynamic_PPI_ratio), dataf)  + 
  scale_color_gradient(low = apple_colors[10], high = apple_colors[5], 
                       #breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03)
                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))+  
  scale_size(range = c(0,3), 
             #breaks = c(0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03)
             breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
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
ggsave("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/all_PPI/heatmap_mean_density_changing_PPI_BP_all_median.pdf", width = 16, height =12)

##### Get barplot of the variation of each compartment
#### Take the mean CV for each GO term across all GO terms
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/dynamic_all_pairwise/")
BP_variation = csvReader_T("Changing_PPI_BP_density_mean_environments_chosen.csv")
CC_variation = csvReader_T("Changing_PPI_CC_density_mean_environments.csv")
GO_BP_order = as.matrix(read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_BP_order_chosen.txt",
                                   header = T, sep = "\t")) 
GO_CC_order = as.matrix(read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_CC_order.txt",
                                   header = T, sep = "\t"))
GO_BP_matrix = cbind(GO_BP_order, rev(as.character(1:59)))
GO_CC_matrix = cbind(GO_CC_order, LETTERS[1:22])
CC_CV = rep(0, nrow(GO_CC_matrix))
BP_CV = rep(0, nrow(GO_BP_matrix))
for(i in 1:length(CC_CV)){
  index = which(CC_variation[,1] == GO_CC_matrix[i,1])
  CC_CV[i] = mean(as.numeric(CC_variation[index, 3]), na.rm = T)
}
for(i in 1:length(BP_CV)){
  index = which(BP_variation[,1] == GO_BP_matrix[i,1])
  BP_CV[i] = mean(as.numeric(BP_variation[index, 3]), na.rm = T)
}

matrix_CC = cbind(GO_CC_matrix, CC_CV)
matrix_BP = cbind(GO_BP_matrix, BP_CV)
colnames(matrix_CC) = c("CC","index","ratio")
colnames(matrix_BP) = c("BP","index","ratio")
#matrix_CC_order = matrix_CC[order(as.numeric(matrix_CC[,3]), decreasing = T),]
#matrix_BP_order = matrix_BP[order(as.numeric(matrix_BP[,3]), decreasing = T),]
#csvWriter(matrix_CC_order, "Changing_PPI_ratio_single_CC_order.csv")
#csvWriter(matrix_BP_order, "Changing_PPI_ratio_single_BP_order.csv")
csvWriter(matrix_CC, "Changing_PPI_ratio_single_CC_primary.csv")
csvWriter(matrix_BP, "Changing_PPI_ratio_single_BP_primary.csv")


### plot the dynamics orderly
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/dynamic_all_pairwise/")
matrix_CC = csvReader_T("Changing_PPI_ratio_single_CC_primary.csv")
matrix_BP = csvReader_T("Changing_PPI_ratio_single_BP_primary.csv")
## CC the same color  
# BP: 1:others, 2: Transcription, 3: RNA: processing, 4:Translation, 5: Ribosome regulation 
#col_chosen = c(apple_colors[1], "#f03b20", "#fd8d3c", "#810f7c", "#8856a7")
pdf("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/all_PPI/Changing_PPI_ratio_CC_primary.pdf", width= 2, height=4.5)
par(mar = c(2,2,0.5,0.5))

barCenter = barplot(as.numeric(matrix_CC[,3]), horiz=T, beside=F, 
                    xlab="Median dynamic PPI ratio",axisnames=F, border=NA, cex.lab = 0.5,
                    col = apple_colors[5], cex.axis = 0.5, ylab = "Cellular compartment")
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(y= barCenter, x = -0.01, labels = matrix_CC[,2], cex = 0.5, xpd = TRUE)

#legend(19,1.2, legend= c("Chromosome", "Nucleolus"), fill = col_chosen[c(2,3)], bty = "n", 
#border = FALSE, xpd = TRUE, cex = 0.6)
dev.off()

pdf("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/all_PPI/Changing_PPI_ratio_BP_primary.pdf", width= 4.8, height=1.8)
par(mar = c(3,2,0.5,0))

barCenter = barplot(as.numeric(matrix_BP[,3]), horiz=F, beside=F,
                    ylab="Coefficient variration",axisnames=F, border=NA, cex.lab = 0.5,
                    col = apple_colors[5], cex.axis = 0.5)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x= barCenter, y = -0.1, labels = matrix_BP[,2], srt = 60, cex = 0.4, xpd = TRUE)
#legend(19,2.8, legend= c("Transcription", "RNA processing", "Translation","Ribosome regulation"), 
       #fill = col_chosen[2:5], bty = "n", border = FALSE, xpd = TRUE, cex = 0.6)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()



############################################################# not use
## 16C
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
pos_count_name = "16C/Changing_PPI_count_CC.txt"
output_pos_count = "16C/Changing_PPI_count_CC_matrix.txt"
output_density = "16C/Changing_PPI_count_CC_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
pos_count_name = "16C/Changing_PPI_count_BP.txt"
output_pos_count = "16C/Changing_PPI_count_BP_matrix.txt"
output_density = "16C/Changing_PPI_count_BP_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

##Dox
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
pos_count_name = "Dox/Changing_PPI_count_CC.txt"
output_pos_count = "Dox/Changing_PPI_count_CC_matrix.txt"
output_density = "Dox/Changing_PPI_count_CC_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
pos_count_name = "Dox/Changing_PPI_count_BP.txt"
output_pos_count = "Dox/Changing_PPI_count_BP_matrix.txt"
output_density = "Dox/Changing_PPI_count_BP_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## FK506
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
pos_count_name = "FK506/Changing_PPI_count_CC.txt"
output_pos_count = "FK506/Changing_PPI_count_CC_matrix.txt"
output_density = "FK506/Changing_PPI_count_CC_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
pos_count_name = "FK506/Changing_PPI_count_BP.txt"
output_pos_count = "FK506/Changing_PPI_count_BP_matrix.txt"
output_density = "FK506/Changing_PPI_count_BP_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## Forskolin
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
pos_count_name = "Forskolin/Changing_PPI_count_CC.txt"
output_pos_count = "Forskolin/Changing_PPI_count_CC_matrix.txt"
output_density = "Forskolin/Changing_PPI_count_CC_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
pos_count_name = "Forskolin/Changing_PPI_count_BP.txt"
output_pos_count = "Forskolin/Changing_PPI_count_BP_matrix.txt"
output_density = "Forskolin/Changing_PPI_count_BP_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## H2O2
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
pos_count_name = "H2O2/Changing_PPI_count_CC.txt"
output_pos_count = "H2O2/Changing_PPI_count_CC_matrix.txt"
output_density = "H2O2/Changing_PPI_count_CC_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
pos_count_name = "H2O2/Changing_PPI_count_BP.txt"
output_pos_count = "H2O2/Changing_PPI_count_BP_matrix.txt"
output_density = "H2O2/Changing_PPI_count_BP_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## HU
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
pos_count_name = "HU/Changing_PPI_count_CC.txt"
output_pos_count = "HU/Changing_PPI_count_CC_matrix.txt"
output_density = "HU/Changing_PPI_count_CC_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
pos_count_name = "HU/Changing_PPI_count_BP.txt"
output_pos_count = "HU/Changing_PPI_count_BP_matrix.txt"
output_density = "HU/Changing_PPI_count_BP_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)

## NaCl
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
pos_count_name = "NaCl/Changing_PPI_count_CC.txt"
output_pos_count = "NaCl/Changing_PPI_count_CC_matrix.txt"
output_density = "NaCl/Changing_PPI_count_CC_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
pos_count_name = "NaCl/Changing_PPI_count_BP.txt"
output_pos_count = "NaCl/Changing_PPI_count_BP_matrix.txt"
output_density = "NaCl/Changing_PPI_count_BP_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)


## Raffinose
# Cellular compartment
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_CC.txt"
pos_count_name = "Raffinose/Changing_PPI_count_CC.txt"
output_pos_count = "Raffinose/Changing_PPI_count_CC_matrix.txt"
output_density = "Raffinose/Changing_PPI_count_CC_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)
# Biological process
all_count_name = "~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/Pos_PPI_count_BP.txt"
pos_count_name = "Raffinose/Changing_PPI_count_BP.txt"
output_pos_count = "Raffinose/Changing_PPI_count_BP_matrix.txt"
output_density = "Raffinose/Changing_PPI_count_BP_density_matrix.txt"
transform_list_matrix(all_count_name, pos_count_name, output_pos_count, output_density)