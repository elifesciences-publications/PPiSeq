# Barplot for PPI number in each environment and show number of PPIs in each bin (found in 1-9 environments)
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
count_summary = csvReader_T("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter.csv") # 12981

environment_matrix = matrix(0, 9, 9)
colnames(environment_matrix) = colnames(count_summary)[3:11]
rownames(environment_matrix) = as.character(1:9)
for(i in 3:11){
        PPI_chosen= count_summary[which(count_summary[,i] == "1"),]
        counts = as.data.frame(table(PPI_chosen[,2]))
        environment_matrix[,i-2] = counts$Freq
}
environment_matrix_order = environment_matrix[,c(1,9,2,3,7,5,6,4,8)]
col_purple = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
pdf("Figures/SFigures/SFigure4/FigureS4A_Number_PPIs_across_environments.pdf", height = 5, width = 5)
par(mar = c(4,4,2,1))
barCenter = barplot(environment_matrix_order, horiz=F, beside=F,  ylab="Number of PPIs",
                    space= c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), yaxt="n", ylim = c(0,7000),
                    col= col_purple, axisnames=F, border=NA)

axis(2,at = c(seq(0, 7000, by = 1000)))
legend(0, 7200, legend=as.character(rev(1:9)), title = c("Number of environments in which \n the PPI is detected"),
       fill=col_purple[rev(1:9)],  bty="n", border=FALSE, xpd = TRUE, ncol = 3)

text(x= barCenter, y = -150, cex =0.8,
     labels = rep(c("SD", "FK506", expression('H'[2]* 'O'[2]), "Hydroxyurea", "NaCl", 
                    "Forskolin",  "Raffinose", "Doxorubicin", "16 \u00B0C"),2), 
     srt= 45, adj = 1, xpd = TRUE)

dev.off()









