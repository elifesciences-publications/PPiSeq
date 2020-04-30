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

setwd("~/Dropbox/PPiSeq_02/")
SD1_lineage = csvReader_T("Paper_data/Lineage_barcode_fitness_files/SD_PPI_barcodes_fitness_counts.csv")
SD2_lineage = csvReader_T("Paper_data/Lineage_barcode_fitness_files/SD2_PPI_barcodes_fitness_counts.csv")
SD1_mean = csvReader_T("Paper_data/PPI_mean_fitness_calling_files/SD_mean_fitness_positive.csv")
SD2_mean = csvReader_T("Paper_data/PPI_mean_fitness_calling_files/SD2_mean_fitness_positive.csv")

##### Normalize fitenss values for each barcode by range of positive and negative fitness values
Fitness_normalization = function(DMSO_lineage){
  DMSO_DHFR_Pos= DMSO_lineage[which(DMSO_lineage[,1] == "positive_DHFR"),]
  DMSO_DHFR_Neg= DMSO_lineage[which(DMSO_lineage[,1] == "negative_non_DHFR"),]
  
  DMSO_DHFR_Pos_mean = mean(as.numeric(DMSO_DHFR_Pos[,4])) 
  DMSO_DHFR_Neg_mean = mean(as.numeric(DMSO_DHFR_Neg[,4])) 
  
  DMSO_lineage_chosen = DMSO_lineage[, c(1,2,4)]
  DMSO_lineage_chosen[,3] = (as.numeric(DMSO_lineage_chosen[,3]) - DMSO_DHFR_Neg_mean)/(DMSO_DHFR_Pos_mean - DMSO_DHFR_Neg_mean)
  return(DMSO_lineage_chosen)
}

##### Put fitness values of the same PPI onto the same line
cluster_same_PPI = function(DMSO_lineage_norm_file, DMSO_multiple){
  PPI_multiple_RRS = DMSO_multiple[grep("Neg_PPI", DMSO_multiple[,1]),] 
  PPI_multiple_PRS = DMSO_multiple[grep("Pos_PPI", DMSO_multiple[,1]),] 
  PPI_multiple_pos = DMSO_multiple[grep("positive_DHFR", DMSO_multiple[,1]),]
  PPI_multiple_neg = DMSO_multiple[grep("negative_non_DHFR", DMSO_multiple[,1]),]
  PPI_multiple_control = rbind(PPI_multiple_RRS, PPI_multiple_PRS, PPI_multiple_pos, PPI_multiple_neg)
  PPI_multiple_select = DMSO_multiple[which(!DMSO_multiple[,1] %in% PPI_multiple_control[,1]),] 
  DMSO_lineage_norm= DMSO_lineage_norm_file[which(DMSO_lineage_norm_file[,1] %in% PPI_multiple_select[,1]), ] 
  
  PPI_unique= unique(PPI_multiple_select[,1])
  PPI_indiv_matrix= matrix(0, length(PPI_unique), 7)
  PPI_indiv_matrix[,1]= PPI_unique
  PPI_indiv_matrix[,2] = as.numeric(PPI_multiple_select[match(as.character(PPI_indiv_matrix[,1]), 
                                                              as.character(PPI_multiple_select[,1])),2])
  index = 0
  for (i in 1:length(PPI_unique)){
    if(as.numeric(PPI_indiv_matrix[i,2]) == 4){
      PPI_indiv_matrix[i, 4]= DMSO_lineage_norm[index + 4 -3, 3]
      PPI_indiv_matrix[i, 5]= DMSO_lineage_norm[index + 4 -2, 3]
      PPI_indiv_matrix[i, 6]= DMSO_lineage_norm[index + 4 -1, 3]
      PPI_indiv_matrix[i, 7]= DMSO_lineage_norm[index + 4, 3] 
      PPI_indiv_matrix[i,3] = mean(as.numeric(PPI_indiv_matrix[i,4:7]))
      index = index + 4
    }else if(as.numeric(PPI_indiv_matrix[i,2] == 3)){
      PPI_indiv_matrix[i, 4]= DMSO_lineage_norm[index + 3 -2, 3]
      PPI_indiv_matrix[i, 5]= DMSO_lineage_norm[index + 3 -1, 3]
      PPI_indiv_matrix[i, 6]= DMSO_lineage_norm[index + 3, 3]
      PPI_indiv_matrix[i,3] = mean(as.numeric(PPI_indiv_matrix[i,4:6]))
      index = index + 3
    }else if(as.numeric(PPI_indiv_matrix[i,2] == 2)){
      PPI_indiv_matrix[i, 4]= DMSO_lineage_norm[index + 2 -1, 3]
      PPI_indiv_matrix[i, 5]= DMSO_lineage_norm[index + 2, 3]
      PPI_indiv_matrix[i,3] = mean(as.numeric(PPI_indiv_matrix[i,4:5]))
      index = index + 2
    }
  }
  colnames(PPI_indiv_matrix)= c("PPI", "Barcodes_count", "Mean_fitness", "fit01", "fit02", "fit03", "fit04")
  return(PPI_indiv_matrix)
}

### normalize each lineage and take the mean fitness for lineages that tag the same PPI
SD1_lineage_norm = Fitness_normalization(SD1_lineage)
SD1_PPI_norm = cluster_same_PPI(SD1_lineage_norm, SD1_mean) # This function has removed the control strains
SD1_pos_PPI = SD1_mean[which(as.numeric(SD1_mean[,8]) == 1), 1]
SD1_PPI_norm_pos = SD1_PPI_norm[which(SD1_PPI_norm[,1] %in% SD1_pos_PPI),]

SD2_lineage_norm = Fitness_normalization(SD2_lineage)
SD2_PPI_norm = cluster_same_PPI(SD2_lineage_norm, SD2_mean) # This function has removed the control strains
SD2_pos_PPI = SD2_mean[which(as.numeric(SD2_mean[,8]) == 1), 1]
SD2_PPI_norm_pos = SD2_PPI_norm[which(SD2_PPI_norm[,1] %in% SD2_pos_PPI),]

all_pos = unique(c(SD1_PPI_norm_pos[,1], SD2_PPI_norm_pos[,1])) # 6185
intersect_pos = match_both_direction(SD1_PPI_norm_pos, SD2_PPI_norm_pos[,1]) # 3637
#intersect_pos_unique = mark_duplicates_fast(intersect_pos[,1]) # 3282

SD1_all_pos = SD1_PPI_norm[which(SD1_PPI_norm[,1] %in% all_pos),] # 6041
SD2_all_pos = SD2_PPI_norm[which(SD2_PPI_norm[,1] %in% all_pos),] # 6079


matrix_SD = matrix(NA, length(all_pos), 3)
matrix_SD[,1] = all_pos
matrix_SD[,2] = as.numeric(SD1_all_pos[match(all_pos, SD1_all_pos[,1]),3])
matrix_SD[,3] = as.numeric(SD2_all_pos[match(all_pos, SD2_all_pos[,1]),3])
matrix_SD_no_na = na.omit(matrix_SD) # 5935
cor(as.numeric(matrix_SD_no_na[,2]), as.numeric(matrix_SD_no_na[,3]), method="spearman") # 0.7321476
matrix_final = data.frame(matrix_SD_no_na[,1], as.numeric(matrix_SD_no_na[,2]), 
                          as.numeric(matrix_SD_no_na[,3]))
colnames(matrix_final) = c("PPI", "SD", "SD2")
csvWriter(matrix_final, "~/Dropbox/PPiSeq_02/Working_data_2/SD_replciates.csv")

matrix_final = dataFrameReader_T("~/Dropbox/PPiSeq_02/Working_data_2/SD_replciates.csv")

library(scales)
library(ggplot2)
ggplot() +
  geom_hex(aes(x= SD, y= SD2, fill = log10(..count..)), matrix_final, bins = 60)+
  scale_fill_gradient(low= "white", high = apple_colors[7])+

  geom_smooth(aes(x = seq(-0.2, 1.2, by = 0.2), y = seq(-0.2, 1.2, by = 0.2)), linetype =2,method='lm', se= FALSE, col= apple_colors[11], cex = 0.3)+
  #annotate("text", x = 0.1, y = 1.1, label = expression(paste("Spearman's ", italic(r), " = 0.73")), 
           #parse = TRUE, col = apple_colors[11], size = 2) +

  scale_y_continuous(name = "Fitness of replicate culture 2",
                     limits=c(-0.2, 1.2),
                     breaks=seq(-0.2,1.2, by =0.2),
                     labels = seq(-0.2,1.2, by= 0.2)) +
  scale_x_continuous(name = "Fitness of replicate culture 1", 
                     limits=c(-0.2, 1.2),
                     breaks=seq(-0.2,1.2, by =0.2),
                     labels = seq(-0.2,1.2, by= 0.2))+
  labs(fill = expression('Log'[10]* '(count)')) +  
  #theme(legend.position =c(0.85,0.3), legend.key=element_blank(),
        #legend.text=element_text(size=5),legend.title=element_text(size=6),
        #legend.key.size = unit(0.3, "cm")) +
  theme(legend.key=element_blank(),
        legend.text=element_text(size=5),legend.title=element_text(size=6),
        legend.key.size = unit(0.3, "cm")) +
  #guides(fill=guide_legend(title="Log10(Count)")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 6, color = "black"),
        axis.text.y.left = element_text(size = 6, color = "black"),
        axis.title=element_text(size= 6))

ggsave("~/Dropbox/PPiSeq_02/working_figure/Figure1/Figure1F_correlation_SD_replicates_hexagonlot.pdf", height = 2.2, width = 3)
