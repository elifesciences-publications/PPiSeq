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

# Figure2E coannotation rate for different groups of PPIs
setwd("~/Dropbox/PPiSeq_02/")
## Write out all PPIs covered by PPiSeq in all environments
PPI_all = csvReader_T("Working_data/Positive_PPI_environment/All_PPI_environments_normalized_fit.csv")
PPI_all_matrix = cbind(PPI_all[,1], rep(1, nrow(PPI_all)))
colnames(PPI_all_matrix) = c("PPI", "no_meaning")
csvWriter(PPI_all_matrix, "Working_data/Positive_PPI_environment/All_PPI_for_coannotation.csv")

### Run python code to check the annotation for positive PPIs and all PPIs

#### Create a matrix for the coannotation rate
PPI_coannotation_MF = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_environment_coannotation_MF.csv")
PPI_coannotation_BP = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_environment_coannotation_BP.csv")
PPI_coannotation_CC = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_environment_coannotation_CC.csv")
coannotation_matrix = matrix(0, 4, 9)
coannotation_matrix[1,] = as.character(1:9)
for(i in 1:9){
        PPI_select_CC = PPI_coannotation_CC[which(as.numeric(PPI_coannotation_CC[,2]) == i),]
        PPI_select_BP = PPI_coannotation_BP[which(as.numeric(PPI_coannotation_BP[,2]) == i),]
        PPI_select_MF = PPI_coannotation_MF[which(as.numeric(PPI_coannotation_MF[,2]) == i),]
        PPI_co_CC = length(which(PPI_select_CC[,3] == "1"))
        PPI_co_BP = length(which(PPI_select_BP[,3] == "1"))
        PPI_co_MF = length(which(PPI_select_MF[,3] == "1"))
        PPI_all_CC = nrow(PPI_select_CC)
        PPI_all_BP = nrow(PPI_select_BP)
        PPI_all_MF = nrow(PPI_select_MF)
        coannotation_matrix[2,i] = PPI_co_CC/PPI_all_CC
        coannotation_matrix[3,i] = PPI_co_BP/PPI_all_BP
        coannotation_matrix[4,i] = PPI_co_MF/PPI_all_MF
        
}

# Check the coannotation rate for all PPIs (neg + pos)
PPI_all_MF = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_all_neg_pos_coannotation_MF.csv")
PPI_all_BP = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_all_neg_pos_coannotation_BP.csv")
PPI_all_CC = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_all_neg_pos_coannotation_CC.csv")

PPI_all_co_CC = length(which(PPI_all_CC[,3] == "1")) # 1075009
PPI_all_co_BP = length(which(PPI_all_BP[,3] == "1")) # 139624
PPI_all_co_MF = length(which(PPI_all_MF[,3] == "1")) # 119164

all_coannotation_CC = PPI_all_co_CC/nrow(PPI_all_CC) # 0.6748457
all_coannotation_BP = PPI_all_co_BP/nrow(PPI_all_BP) # 0.08765011
all_coannotation_MF = PPI_all_co_MF/nrow(PPI_all_MF) # 0.07480618

coannotation_matrix = coannotation_matrix[2:4,]
co_rate = as.numeric(c(coannotation_matrix[,1], coannotation_matrix[,2], coannotation_matrix[,3],
            coannotation_matrix[,4], coannotation_matrix[,5], coannotation_matrix[,6],
            coannotation_matrix[,7], coannotation_matrix[,8], coannotation_matrix[,9]))

pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2E_coannotation_rate_bar_plot.pdf", width= 6, height=5)
barCenter = barplot(co_rate*100, horiz=F, beside=F, ylim=c(0,100), ylab="Fraction co-annotated (%)",
                    space= c(0.4, 0.08, 0.08, 0.4, 0.08, 0.08, 0.4, 0.08, 0.08, 
                             0.4, 0.08, 0.08, 0.4, 0.08, 0.08, 0.4, 0.08, 0.08,
                             0.4, 0.08, 0.08, 0.4, 0.08, 0.08, 0.4, 0.08, 0.08),
                    col= apple_colors[c(5,3,7)] , axisnames=F, border=NA, cex.axis=0.8)
lines(1:32, rep(all_coannotation_CC*100, 32), lty = 2, col = apple_colors[11], lwd = 1.5)
lines(1:32, rep(all_coannotation_BP*100, 32), lty = 2, col = apple_colors[11], lwd = 1.5)
lines(1:32, rep(all_coannotation_MF*100, 32), lty = 2, col = apple_colors[11], lwd = 1.5)
text(c(33,33,33), c(all_coannotation_CC*100, (all_coannotation_BP*100+2), (all_coannotation_MF*100 -2)),
     labels = c("CC", "BP", "MF"), col = apple_colors[c(5,3,7)], xpd = TRUE, cex = 0.8)
legend(-0.5,110, legend=c("Cellular compartment (CC)", "Biological process (BP)", "Molecular function (MF)"),
       fill=apple_colors[c(5,3,7)], cex=0.8, bty="n",border=FALSE, xpd = TRUE)

env_num_loc = rep(0, 9)
for(i in 1:9){
        env_num_loc[i] = mean(barCenter[(3*i-2):(3*i)])
}
text(x = env_num_loc, y = -8, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()
