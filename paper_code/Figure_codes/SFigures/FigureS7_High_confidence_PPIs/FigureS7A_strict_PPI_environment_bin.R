setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used fuctions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

### barplot to show how many of them have been reproted
all_PPI_matrix_final = csvReader_T("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter_strict_threshold.csv")
reported_PPI = csvReader_T("Outsourced_datasets/BIOGRID//multiple_validated_PPI.csv")
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
matrix_PPI_env_rep[3,] # 255 110  83  60  76  89 127 220 161
all_PPI_count = matrix_PPI_env_rep[1,] + matrix_PPI_env_rep[2,] + matrix_PPI_env_rep[3,]
all_PPI_count # 3287  571  283  220  207  208  273  388  344
ratio_BioGrid = matrix_PPI_env_rep[3,]/all_PPI_count
ratio_BioGrid # 0.07757834 0.19264448 0.29328622 0.27272727 0.36714976 0.42788462 0.46520147 0.56701031 0.46802326
ratio_BioGrid_reported = c("7.8%", "19.3%", "29.3%", "27.3%", "36.7%", "42.8%", "46.5%", "56.7%", "46.8%")
ratio_PCA_low = matrix_PPI_env_rep[2,]/all_PPI_count
ratio_PCA_low # 0.08001217 0.24518389 0.33922261 0.37727273 0.39613527 0.38942308 0.42124542 0.34278351 0.45639535
ratio_PCA_low_overlapped = c("8%", "24.5%", "33.9%", "37.7%", "39.6%", "38.9%", "42.1%", "34.3%", "45.6%")

library(RColorBrewer)
col_chosen = apple_colors[c(1,3,4)]
pdf("Figures/SFigures/SFigure7/FigureS7A_Number_environments_PPI_reproted.pdf", height = 5, width = 5)
par(mar = c(3,5,2,1))
barCenter = barplot(matrix_PPI_env_rep, horiz=F, beside=F, ylim=c(0,4000), ylab="Number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
                    col= col_chosen, axisnames=F, border=NA)
legend("topright", legend=c("BioGRID", "mDHFR-PCA (80% < PPV < 98.2%)", "Previously unreported"), 
       fill=col_chosen[c(3,2,1)], cex = 0.8, bty="n", border=FALSE)
text(x= barCenter, y = all_PPI_count + 100, labels = ratio_PCA_low_overlapped, 
     cex=0.7, xpd = TRUE, col= col_chosen[2]) 
text(x= barCenter, y = all_PPI_count + 200, labels = ratio_BioGrid_reported, 
     cex=0.7, xpd = TRUE, col= col_chosen[3]) 
text(x= barCenter, y = -100, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -300, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()






