#### Take the mean CV for each GO term across all GO terms
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
CC_variation = csvReader_T("Figure2_related_data_generated_by_python_script/Variation_CC_all_environments.csv")
BP_variation = csvReader_T("Figure2_related_data_generated_by_python_script/Variation_BP_all_environments_chosen.csv")
GO_BP_order = as.matrix(read.table("Figures/Figure2/GO_BP_order_chosen.txt",header = T, sep = "\t")) 
GO_CC_order = as.matrix(read.table("Figures/Figure2/GO_CC_order.txt",header = T, sep = "\t"))

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
colnames(matrix_CC) = c("CC","index","CV")
colnames(matrix_BP) = c("BP","index","CV")

pdf("Figures/Figure2/Figure2_CC_dynamics_primary.pdf", width= 4, height= 4)
barCenter = barplot(as.numeric(matrix_CC[,3]), horiz=T, beside=F, 
                    xlab="Mean CV",axisnames=F, border=NA, cex.lab = 0.5,
                    col = apple_colors[5], cex.axis = 0.5, ylab = "Cellular compartment")
text(y= barCenter, x = -0.1, labels = matrix_CC[,2], cex = 0.5, xpd = TRUE)
dev.off()

pdf("Figures/Figure2/Figure2_BP_dynamics_primary.pdf", width= 6, height=4)
barCenter = barplot(as.numeric(matrix_BP[,3]), horiz=F, beside=F,
                    ylab="Mean CV",axisnames=F, border=NA, cex.lab = 0.5,
                    col = apple_colors[5], cex.axis = 0.5)
text(x= barCenter, y = -0.15, labels = matrix_BP[,2], srt = 60, cex = 0.4, xpd = TRUE)
dev.off()

