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
# first extract all the negative PP-pairs from all environments
DMSO = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv")[,1]
H2O2 = csvReader_T("Paper_data/H2O2_mean_fitness_positive.csv")[,1]
HU = csvReader_T("Paper_data/Hydroxyurea_mean_fitness_positive.csv")[,1]
Dox = csvReader_T("Paper_data/Doxorubicin_mean_fitness_positive.csv")[,1]
Forskolin = csvReader_T("Paper_data/Forskolin_mean_fitness_positive.csv")[,1]
Raffinose = csvReader_T("Paper_data/Raffinose_mean_fitness_positive.csv")[,1]
NaCl = csvReader_T("Paper_data/NaCl_mean_fitness_positive.csv")[,1]
cold = csvReader_T("Paper_data/Cold_16C_mean_fitness_positive.csv")[,1]
FK506 = csvReader_T("Paper_data/FK506_mean_fitness_positive.csv")[,1]

all_PPI = unique(c(DMSO, H2O2, HU, Dox, Forskolin, Raffinose, NaCl, cold, FK506)) # 1593177
all_PPI_dup = mark_duplicates_fast(all_PPI) # 1469024
pos_PPI = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
all_PPI_pos = match_both_direction(all_PPI_dup, pos_PPI[,1]) # 13829
all_PPI_neg = all_PPI_dup[which(!all_PPI_dup[,1] %in% all_PPI_pos[,1]),] # 1455195
all_PPI_neg[,2] = 1 # make the count be 1
csvWriter(all_PPI_neg, "Working_data/Positive_PPI_environment/Neg_PPIs_all_environment.csv")

### check the annotation for these negative PPIs

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

# Check the coannotation rate for negative PPIs
PPI_negative_MF = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_neg_coannotation_MF.csv")
PPI_negative_BP = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_neg_coannotation_BP.csv")
PPI_negative_CC = csvReader_T("Working_data/Positive_PPI_environment/GO_coannotation/PPI_neg_coannotation_CC.csv")

PPI_neg_co_CC = length(which(PPI_negative_CC[,3] == "1")) # 971461
PPI_neg_co_BP = length(which(PPI_negative_BP[,3] == "1")) # 133768
PPI_neg_co_MF = length(which(PPI_negative_MF[,3] == "1")) # 179303

neg_coannotation_CC = PPI_neg_co_CC/nrow(PPI_negative_CC) # 0.6675813
neg_coannotation_BP = PPI_neg_co_BP/nrow(PPI_negative_BP) # 0.09192445
neg_coannotation_MF = PPI_neg_co_MF/nrow(PPI_negative_MF) # 0.1232158

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
lines(1:32, rep(neg_coannotation_CC*100, 32), lty = 2, col = apple_colors[11], lwd = 1.5)
lines(1:32, rep(neg_coannotation_BP*100, 32), lty = 2, col = apple_colors[11], lwd = 1.5)
lines(1:32, rep(neg_coannotation_MF*100, 32), lty = 2, col = apple_colors[11], lwd = 1.5)
text(c(33,33,33), c(neg_coannotation_CC*100, (neg_coannotation_BP*100-1), (neg_coannotation_MF*100 + 1)),
     labels = c("CC", "BP", "MF"), col = apple_colors[11], xpd = TRUE, cex = 0.8)
legend(-0.5,110, legend=c("Cellular compartment (CC)", "Biological process (BP)", "Molecular function (MF)"),
       fill=apple_colors[c(5,3,7)], cex=0.8, bty="n",border=FALSE, xpd = TRUE)

env_num_loc = rep(0, 9)
for(i in 1:9){
        env_num_loc[i] = mean(barCenter[(3*i-2):(3*i)])
}
text(x = env_num_loc, y = -8, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is observed", xpd = TRUE)
dev.off()
