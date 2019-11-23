##### Only plot average line for each group
col_purple = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")

mean_count_all = data.frame(col_mean_1, col_mean_2, col_mean_3, col_mean_4, col_mean_5,
                            col_mean_6, col_mean_7, col_mean_8, col_mean_9)
#pdf("Working_figure/Figure2/Figure2B_PPI_number_accumulation_environments_combine.pdf", width= 5, height =5)
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2B_cumulative_PPI_number/Figure2B_PPI_number_accumulation_environments_combine_average.pdf", width =5, height =5)
par(mar = c(4,4,2,4))
plot(1:9, col_mean_1, xlim = c(1,9), ylim = c(0,14000), type = "l",
     col = col_purple[1],lwd = 2,
     xlab = "Number of environments assayed",
     ylab = "Number of PPIs", main = NA, bty = "n")
axis(2, at= seq(0, 14000, by = 2000), labels = seq(0, 14000, by = 2000))
axis(1, at= 1:9, labels = 1:9)
for(i in 2:ncol(mean_count_all)){
        lines(i:9, mean_count_all[,i][i:9], col = col_purple[i], lwd = 2)
}
legend(1,14000, as.character(1:9),lty = 1, col = col_purple, ncol = 3, bty= "n")
text(4, 14500, "Minimum number of environments\nin which the PPI is observed", xpd = TRUE)
dev.off()

###### Check each PPI group (One, two, three, four, five, six, seven, eight, nine environment only PPIs)

setwd("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Accumulated_PPI/exact_number/")
matrix_count_1 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_1.csv")
matrix_count_2 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_2.csv")
matrix_count_3 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_3.csv")
matrix_count_4 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_4.csv")
matrix_count_5 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_5.csv")
matrix_count_6 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_6.csv")
matrix_count_7 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_7.csv")
matrix_count_8 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_8.csv")
matrix_count_9 = dataFrameReader_T("PPI_number_environment_accumulation_all_order_9.csv")

col_mean_1 = colMeans(matrix_count_1)
col_mean_2 = colMeans(matrix_count_2)
col_mean_3 = colMeans(matrix_count_3)
col_mean_4 = colMeans(matrix_count_4)
col_mean_5 = colMeans(matrix_count_5)
col_mean_6 = colMeans(matrix_count_6)
col_mean_7 = colMeans(matrix_count_7)
col_mean_8 = colMeans(matrix_count_8)
col_mean_9 = colMeans(matrix_count_9)

mean_count_all = data.frame(col_mean_1, col_mean_2, col_mean_3, col_mean_4, col_mean_5,
                            col_mean_6, col_mean_7, col_mean_8, col_mean_9)
#pdf("Working_figure/Figure2/Figure2B_PPI_number_accumulation_environments_combine.pdf", width= 5, height =5)
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2B_cumulative_PPI_number/Figure2B_PPI_number_accumulation_environments_differen_group_average.pdf", width =5, height =5)
par(mar = c(4,4,2,4))
plot(1:9, col_mean_1, xlim = c(1,9), ylim = c(0,7000), type = "l",
     col = col_purple[1],lwd = 2,
     xlab = "Number of environments assayed",
     ylab = "Number of PPIs", main = NA, bty = "n")
axis(2, at= seq(0,7000, by = 1000), labels = seq(0, 7000, by = 1000))
axis(1, at= 1:9, labels = 1:9)
for(i in 2:ncol(mean_count_all)){
        lines(i:9, mean_count_all[,i][i:9], col = col_purple[i], lwd = 2)
}
legend(1,7000, as.character(1:9),lty = 1, col = col_purple, ncol = 3, bty= "n")
text(4, 7300, "Number of environments\nin which the PPI is observed", xpd = TRUE)
dev.off()

###### Use a new way to check these