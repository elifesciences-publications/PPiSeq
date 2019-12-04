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
CC_variation = dataFrameReader_T("Variation_CC_all_environments.csv")
BP_variation = dataFrameReader_T("Variation_BP_all_environments.csv")
CC_unique = unique(CC_variation$rowv)
BP_unique = unique(BP_variation$rowv)
CC_CV = rep(0, length(CC_unique))
BP_CV = rep(0, length(BP_unique))
for(i in 1:length(CC_unique)){
        index = which(CC_variation$rowv == CC_unique[i])
        CC_CV[i] = mean(CC_variation[index, 3], na.rm = T)
}
for(i in 1:length(BP_unique)){
        index = which(BP_variation$rowv == BP_unique[i])
        BP_CV[i] = mean(BP_variation[index, 3], na.rm = T)
}

matrix_CC = cbind(as.character(CC_unique), CC_CV)
matrix_BP = cbind(as.character(BP_unique), BP_CV)
colnames(matrix_CC) = c("CC", "CV")
colnames(matrix_BP) = c("BP", "CV")
matrix_CC_order = matrix_CC[order(matrix_CC[,2], decreasing = T),]
matrix_BP_order = matrix_BP[order(matrix_BP[,2], decreasing = T),]
csvWriter(matrix_CC_order, "Variation_single_CC_order.csv")
csvWriter(matrix_BP_order, "Variation_single_BP_order.csv")

### plot the top 10 most dynamic GO terms
pdf("~/Dropbox/PPiseq_02/Working_figure/Figure2_PPI_enrichment_GO/Figure2C_CC_dynamics_ordered.pdf", width= 5.5, height=5)
par(mar = c(7,4,1,1))
barCenter = barplot(as.numeric(matrix_CC_order[,2]), horiz=F, beside=F, ylim=c(0,1.5), 
                    ylab="Coefficient variration",
                    #space= c(0.15, 0.15,  0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
                    axisnames=F, border=NA, col = apple_colors[5], cex.axis=0.8)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x = barCenter, y = -0.05, labels = as.character(c("A", "B", "C", "D", "E", "F", "G",
                                                       "H", "I", "G", "K", "L", "M", "N",
                                                       "O", "P", "Q", "R", "S", "T", "U",
                                                       "V")), cex = 0.8, xpd = TRUE)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

pdf("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/FigureS4_PPI_enrichment_GO/Top10_dynamic_BPs.pdf", width= 6, height=6)
par(mar = c(10,5,1,1))
barCenter = barplot(as.numeric(matrix_BP_order[1:10,2]), horiz=F, beside=F, ylim=c(0,1.5), 
                    ylab="Coefficient variration",
                    #space= c(0.15, 0.15,  0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15),
                    axisnames=F, border=NA, col = apple_colors[5], cex.axis=0.7)
#text(x= barCenter, y = as.numeric(merge_ratio)*100, labels = counts_label, cex=0.8, xpd = TRUE)
text(x = barCenter, y = -0.05, labels = as.character(matrix_BP_order[1:10,1]), 
     adj = 1, srt =45, xpd = TRUE, cex = 0.6)
#text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()