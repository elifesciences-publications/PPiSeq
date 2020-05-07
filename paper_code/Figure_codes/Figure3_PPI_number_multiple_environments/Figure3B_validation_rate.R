setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
all_Tecan = csvReader_T("Growth_curve_validation_data/SD_environment_validation_summary.csv")
reported_PPI = csvReader_T("Outsourced_datasets/BIOGRID/multiple_validated_PPI.csv")
all_Tecan_rep = match_both_direction(all_Tecan, reported_PPI[,1]) # 145
all_Tecan_unrep = all_Tecan[which(!all_Tecan[,1] %in% all_Tecan_rep[,1]),] # 357
#split validated PPIs into different groups

rep_PPI_matrix = matrix(0, 4,9)
unrep_PPI_matrix = matrix(0, 4, 9)

for(i in 1:9){
  PPI_select = PPI_group[which(as.numeric(PPI_group[,2]) == i),1]
  reported_PPI_select = match_both_direction(all_Tecan_rep, PPI_select)
  if(length(reported_PPI_select) > 11){
    validate_PPI = length(which(as.numeric(reported_PPI_select[,11]) <= 0.05))
    non_val_PPI = length(which(as.numeric(reported_PPI_select[,11]) > 0.05))
    rep_PPI_matrix[1,i] = validate_PPI
    rep_PPI_matrix[2,i] = non_val_PPI
    rep_PPI_matrix[3,i] = validate_PPI + non_val_PPI
    rep_PPI_matrix[4,i] = rep_PPI_matrix[1,i]/rep_PPI_matrix[3,i]
  }else{
    validate_PPI = length(which(as.numeric(reported_PPI_select[11]) <= 0.05))
    non_val_PPI = length(which(as.numeric(reported_PPI_select[11]) > 0.05))
    rep_PPI_matrix[1,i] = validate_PPI
    rep_PPI_matrix[2,i] = non_val_PPI
    rep_PPI_matrix[3,i] = validate_PPI + non_val_PPI
    rep_PPI_matrix[4,i] = rep_PPI_matrix[1,i]/rep_PPI_matrix[3,i]
  }
  
  
  unrep_PPI_select = match_both_direction(all_Tecan_unrep, PPI_select)
  #unrep_PPI_select = all_Tecan_unrep[which(all_Tecan_unrep[,1] %in% PPI_select),]
  validate_PPI_unrep = length(which(as.numeric(unrep_PPI_select[,11]) <= 0.05))
  non_val_PPI_unrep = length(which(as.numeric(unrep_PPI_select[,11]) > 0.05))
  unrep_PPI_matrix[1,i] = validate_PPI_unrep
  unrep_PPI_matrix[2,i] = non_val_PPI_unrep
  unrep_PPI_matrix[3,i] = validate_PPI_unrep + non_val_PPI_unrep
  unrep_PPI_matrix[4,i] = unrep_PPI_matrix[1,i]/unrep_PPI_matrix[3,i]
}
colnames(rep_PPI_matrix) = c("Envir_1", "Envir_2", "Envir_3", "Envir_4", "Envir_5",
                             "Envir_6", "Envir_7", "Envir_8", "Envir_9")
rownames(rep_PPI_matrix) = c("Validate", "Nonvalidate", "Total", "Val_ratio")
csvWriter_rownames(rep_PPI_matrix, "Growth_curve_validation_data/Reported_validation_matrix_SD_merge.csv")

colnames(unrep_PPI_matrix) = c("Envir_1", "Envir_2", "Envir_3", "Envir_4", "Envir_5",
                               "Envir_6", "Envir_7", "Envir_8", "Envir_9")
rownames(unrep_PPI_matrix) = c("Validate", "Nonvalidate", "Total", "Val_ratio")
csvWriter_rownames(unrep_PPI_matrix, "Growth_curve_validation_data/Unreported_validation_matrix_SD_merge.csv")


#############################################
# Merge the reported and unreported PPI and receck the number
rep_PPI_matrix = dataFrameReader_T("Growth_curve_validation_data/Reported_validation_matrix_SD_merge.csv")
unrep_PPI_matrix = dataFrameReader_T("Growth_curve_validation_data/Unreported_validation_matrix_SD_merge.csv")

merge_validate = rep_PPI_matrix[1, 2:ncol(rep_PPI_matrix)] + unrep_PPI_matrix[1, 2:ncol(rep_PPI_matrix)]
merge_nonvalidate = rep_PPI_matrix[2, 2:ncol(rep_PPI_matrix)] + unrep_PPI_matrix[2, 2:ncol(rep_PPI_matrix)]
merge_sum = merge_validate + merge_nonvalidate

merge_ratio = merge_validate/merge_sum
merge_validate #  55      22      26      20      42      52      52      78      42
merge_sum      # 104      33      37      29      52      63      59      82      43
counts_label = c("55/104", "22/33", "26/37", "20/29", "42/52", "52/63",
                  "52/59", "78/82", "42/43")
library(RColorBrewer)
col_chosen = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
pdf("Figures/Figure3/Figure3B_Validation_bar_plot_merge_calling_all.pdf", width= 5, height=5)
par(mar = c(5,4,1,1)) 
barCenter = barplot(as.numeric(merge_ratio) *100, horiz=F, beside=F, ylim=c(0,100), ylab="Validation rate (%)",
                    space= c(0.15, 0.15,  0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15),axisnames=F, border=NA, 
                    col = col_chosen, cex.axis=0.8)
text(x= barCenter, y = as.numeric(merge_ratio)*100 + 2, labels = counts_label, cex=0.8, xpd = TRUE)
text(x = barCenter, y = -8, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

