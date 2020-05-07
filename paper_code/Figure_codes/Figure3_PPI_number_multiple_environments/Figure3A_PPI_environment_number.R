### barplot to show how many of them have been reproted
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
all_PPI_matrix_final = csvReader_T("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter.csv")
reported_PPI = csvReader_T("Outsourced_datasets/BIOGRID/multiple_validated_PPI.csv")
PCA_lower = as.matrix(read.table("Outsourced_datasets/Tarassov_PPI_PPV_80.txt", header= T, sep = "\t"))
PCA_lower_PPI = paste(PCA_lower[,1], PCA_lower[,4], sep = "_")
PCA_lower_PPI_matrix = cbind(PCA_lower_PPI, rep(1, length(PCA_lower_PPI)))
PCA_lower_PPI_reported = match_both_direction(PCA_lower_PPI_matrix, reported_PPI[,1])
PCA_lower_PPI_unreported = PCA_lower_PPI_matrix[which(!PCA_lower_PPI_matrix[,1] %in% PCA_lower_PPI_reported[,1]),]

matrix_PPI_env_rep = matrix(0, 3, 9)
for(i in 1:9){
  all = all_PPI_matrix_final[which(as.numeric(all_PPI_matrix_final[,2]) == i),]
  all_reported_PCA_low = match_both_direction(all, PCA_lower_PPI_unreported[,1])
  all_reported_BioGrid = match_both_direction(all,reported_PPI[,1])
  all_reported = c(all_reported_PCA_low[,1], all_reported_BioGrid[,1])
  all_unreported = all[which(!all[,1] %in% all_reported),]
  matrix_PPI_env_rep[1,i] = nrow(all_unreported)
  matrix_PPI_env_rep[2,i] = nrow(all_reported_PCA_low)
  matrix_PPI_env_rep[3,i] = nrow(all_reported_BioGrid)
}
matrix_PPI_env_rep[3,] # 278 143 136 148 151 174 208 304 250
all_PPI_count = matrix_PPI_env_rep[1,] + matrix_PPI_env_rep[2,] + matrix_PPI_env_rep[3,]
all_PPI_count # 7724 1266  730  579  579  553  497  579  474
ratio_BioGrid = matrix_PPI_env_rep[3,]/all_PPI_count
ratio_BioGrid # 0.03599171 0.11295419 0.18630137 0.25561313 0.26079447 0.31464738 0.41851107 0.52504318 0.52742616
ratio_BioGrid_reported = c("3.6%", "11.3%", "18.6%", "25.6%", "26.1%", "31.5%", "41.9%", "52.5%", "52.7%")
ratio_PCA_low = matrix_PPI_env_rep[2,]/all_PPI_count
ratio_PCA_low # 0.03560331 0.11848341 0.20821918 0.25215889 0.31951641 0.35443038 0.38631791 0.35924007 0.38607595
ratio_PCA_low_overlapped = c("3.6%", "11.8%", "20.8%", "25.2%", "32.0%", "35.4%", "38.6%", "35.9%", "38.6%")

library(RColorBrewer)
col_chosen = apple_colors[c(5,3,7)]
pdf("Figures/Figure3/Figure3A_Number_environments_PPI_reproted.pdf", height = 4, width = 5)
barCenter = barplot(matrix_PPI_env_rep, horiz=F, beside=F, ylim=c(0,8000), #ylab="Number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
                    col= col_chosen,  ylab = "Number of PPIs", border=NA, cex.axis = 0.6)
legend("topright", legend=c("BioGRID + mDHFR-PCA(PPV > 98.2%)", "mDHFR-PCA (PPV > 80%)", "Previously unreported"), 
       fill=col_chosen[c(3,2,1)], cex = 0.6, bty="n", border=FALSE)
text(x= barCenter, y = all_PPI_count + 200, labels = ratio_PCA_low_overlapped, 
     cex=0.6, xpd = TRUE, col= col_chosen[2]) 
text(x= barCenter, y = all_PPI_count + 500, labels = ratio_BioGrid_reported, 
     cex=0.6, xpd = TRUE, col= col_chosen[3]) 
text(x= barCenter, y = -300, labels = as.character(1:9), xpd = TRUE, cex = 0.6)
text(median(barCenter), y = -800, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()






