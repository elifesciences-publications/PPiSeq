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

#### Take the mean CV for each GO term across all GO terms
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/")
CC_variation = csvReader_T("Variation_CC_all_environments.csv")
BP_variation = csvReader_T("Variation_BP_all_environments.csv")
GO_BP_order = as.matrix(read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_BP_order_chosen.txt",
                                   header = T, sep = "\t")) 
GO_CC_order = as.matrix(read.table("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/GO_CC_order.txt",
                                   header = T, sep = "\t"))
GO_BP_matrix = cbind(GO_BP_order, as.character(1:59))
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
matrix_CC_order = matrix_CC[order(as.numeric(matrix_CC[,3]), decreasing = T),]
matrix_BP_order = matrix_BP[order(as.numeric(matrix_BP[,3]), decreasing = T),]
csvWriter(matrix_CC_order, "Variation_single_CC_order.csv")
csvWriter(matrix_BP_order, "Variation_single_BP_order.csv")


### plot the top 10 most dynamic GO terms
setwd("~/Dropbox/PPiseq_02/Working_data_2/PPI_pair_GO/environment/")
matrix_CC_order = csvReader_T("Variation_single_CC_order_color.csv")
matrix_BP_order = csvReader_T("Variation_single_BP_order_color.csv")
# CC: 1: others, 2: chromosome, 3: nucleolus,  
# BP: 1:others, 2: Transcription, 3: Ribosome regulation, 4: RNA: processing, 5:Translation,
col_chosen = c("#8E8E93", "#386cb0", "#ffff99", "#fdc086", "#beaed4" )
pdf("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/Figure2C_CC_dynamics_ordered.pdf", width= 4, height=2)
par(mar = c(2,2,1,1))

barCenter = barplot(as.numeric(matrix_CC_order[,3]), horiz=F, beside=F, ylim=c(0,1.5), 
                    ylab="Coefficient variration",axisnames=F, border=NA, cex.lab = 0.5,
                    col = col_chosen[as.numeric(matrix_CC_order[,4])], cex.axis = 0.5)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x= barCenter, y = -0.08, labels = matrix_CC_order[,2], cex = 0.5, xpd = TRUE)
legend(19,1.2, legend= c("Chromosome", "Nucleolus"), fill = col_chosen[c(2,3)], bty = "n", 
       border = FALSE, xpd = TRUE, cex = 0.6)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

pdf("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/Figure2D_BP_dynamics_ordered.pdf", width= 5, height=2)
par(mar = c(1,2,0.5,0.5))

barCenter = barplot(as.numeric(matrix_BP_order[,3]), horiz=F, beside=F, ylim=c(0,3), 
                    ylab="Coefficient variration",axisnames=F, border=NA, cex.lab = 0.5,
                    col = col_chosen[as.numeric(matrix_BP_order[,4])], cex.axis = 0.5)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x= barCenter, y = -0.15, labels = matrix_BP_order[,2], srt = 60, cex = 0.4, xpd = TRUE)
legend(25,3.0, legend= c("Transcription", "Ribosome regulation", "RNA processing",
                         "Translation"), fill = col_chosen[2:5], bty = "n", 
       border = FALSE, xpd = TRUE, cex = 0.6)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

### Split the data into two parts
pdf("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/Figure2D_BP_dynamics_ordered_1.pdf", width= 4.5, height=2)
par(mar = c(1,2,0.5,0.5))

barCenter = barplot(as.numeric(matrix_BP_order[1:30,3]), horiz=F, beside=F, ylim=c(0,3), 
                    ylab="Coefficient variration",axisnames=F, border=NA, cex.lab = 0.5,
                    col = col_chosen[as.numeric(matrix_BP_order[1:30,4])], cex.axis = 0.5)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x= barCenter, y = -0.15, labels = matrix_BP_order[1:30,2],  cex = 0.5, xpd = TRUE)
legend(25,3.0, legend= c("Transcription", "Ribosome regulation", "RNA processing",
                         "Translation"), fill = col_chosen[2:5], bty = "n", 
       border = FALSE, xpd = TRUE, cex = 0.6)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

pdf("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/Figure2D_BP_dynamics_ordered_2.pdf", width= 3, height=2)
par(mar = c(1,2,0.5,0.5))

barCenter = barplot(as.numeric(matrix_BP_order[31:59,3]), horiz=F, beside=F, ylim=c(0,3), 
                    ylab="Coefficient variration",axisnames=F, border=NA, cex.lab = 0.5,
                    col = col_chosen[as.numeric(matrix_BP_order[31:59,4])], cex.axis = 0.5)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x= barCenter, y = -0.15, labels = matrix_BP_order[31:59,2],  cex = 0.5, xpd = TRUE)
#legend(30,3.0, legend= c("Transcription", "Ribosome regulation", "RNA processing",
                         #"Translation"), fill = col_chosen[2:5], bty = "n", border = FALSE, xpd = TRUE, cex = 0.6)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()