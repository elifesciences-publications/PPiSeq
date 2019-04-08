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
## Figure2C A barplot to show the validation rate
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/")
### Combine all the chosen PPIs 
d_1_1 = csvReader_T("Diploid_01_01_both_T5.csv")
d_1_2 = csvReader_T("Diploid_01_02_both_T14.csv")
d_1_3 = csvReader_T("Diploid_01_03_both_T15.csv")
d_1_4 = csvReader_T("Diploid_01_04_both_T2.csv")
d_1_5 = csvReader_T("Diploid_01_05_both_T5.csv")
d_1_6 = csvReader_T("Diploid_01_06_both_T7.csv")
d_1_7 = csvReader_T("Diploid_01_07_both_T14.csv")

d_2_1 = csvReader_T("Diploid_02_03_04_01_both_T15.csv")
d_2_2 = csvReader_T("Diploid_02_03_04_02_T7_T2.csv")
d_2_3 = csvReader_T("Diploid_02_03_04_03_both_T5.csv")
d_2_4 = csvReader_T("Diploid_02_03_04_04_both_T7.csv")
d_2_5 = csvReader_T("Diploid_02_03_04_05_both_T14.csv")
d_2_6 = csvReader_T("Diploid_02_03_04_06_both_T15.csv")
d_2_7 = csvReader_T("Diploid_02_03_04_07_both_T2.csv")
d_2_8 = csvReader_T("Diploid_02_03_04_08_both_T5.csv")
d_2_9 = csvReader_T("Diploid_02_03_04_09_both_T7.csv")
d_2_10 = csvReader_T("Diploid_02_03_04_10_both_T14.csv")

d_5_1 = csvReader_T("Diploid_05_01_both_T15.csv")
d_5_2 = csvReader_T("Diploid_05_02_both_T2.csv")
d_5_3 = csvReader_T("Diploid_05_03_both_T5.csv")
d_5_4 = csvReader_T("Diploid_05_04_both_T7.csv")
d_5_5 = csvReader_T("Diploid_05_05_both_T14.csv")
d_5_6 = csvReader_T("Diploid_05_06_both_T15.csv")
d_5_7 = csvReader_T("Diploid_05_07_both_T2.csv")

all_Tecan = rbind(d_1_1, d_1_2, d_1_3, d_1_4, d_1_5, d_1_6, d_1_7,
                  d_2_1, d_2_2, d_2_3, d_2_4, d_2_5, d_2_6, d_2_7,
                  d_2_8, d_2_9, d_2_10, d_5_1, d_5_2, d_5_3, d_5_4, 
                  d_5_5, d_5_6, d_5_7)
reported_PPI = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/multiple_validated_PPI.csv")
all_Tecan_rep = match_both_direction(all_Tecan, reported_PPI[,1]) # 185
all_Tecan_unrep = all_Tecan[which(!all_Tecan[,1] %in% all_Tecan_rep[,1]),] # 510
#split validated PPIs into different groups
PPI_group = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
length(which(all_Tecan[,1] %in% PPI_group[,1])) # 695 all in

rep_PPI_matrix = matrix(0, 4,9)
unrep_PPI_matrix = matrix(0, 4, 9)

for(i in 1:9){
  PPI_select = PPI_group[which(as.numeric(PPI_group[,2]) == i),1]
  reported_PPI_select = all_Tecan_rep[which(all_Tecan_rep[,1] %in% PPI_select),]
  validate_PPI = length(which(as.numeric(reported_PPI_select[,11]) <= 0.05))
  non_val_PPI = length(which(as.numeric(reported_PPI_select[,11]) > 0.05))
  rep_PPI_matrix[1,i] = validate_PPI
  rep_PPI_matrix[2,i] = non_val_PPI
  rep_PPI_matrix[3,i] = validate_PPI + non_val_PPI
  rep_PPI_matrix[4,i] = rep_PPI_matrix[1,i]/rep_PPI_matrix[3,i]
  
  unrep_PPI_select = all_Tecan_unrep[which(all_Tecan_unrep[,1] %in% PPI_select),]
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
csvWriter_rownames(rep_PPI_matrix, "Reported_validation_matrix.csv")

colnames(unrep_PPI_matrix) = c("Envir_1", "Envir_2", "Envir_3", "Envir_4", "Envir_5",
                               "Envir_6", "Envir_7", "Envir_8", "Envir_9")
rownames(unrep_PPI_matrix) = c("Validate", "Nonvalidate", "Total", "Val_ratio")
csvWriter_rownames(unrep_PPI_matrix, "Unreported_validation_matrix.csv")

### Make barplot to show the percentage
### put the reported and unreported on to the same figure
ratio_rep = rep_PPI_matrix[4,]
ratio_unrep = unrep_PPI_matrix[4,]
ratio_all = c(ratio_rep[1], ratio_unrep[1], ratio_rep[2], ratio_unrep[2], 
              ratio_rep[3], ratio_unrep[3], ratio_rep[4], ratio_unrep[4],
              ratio_rep[5], ratio_unrep[5], ratio_rep[6], ratio_unrep[6],
              ratio_rep[7], ratio_unrep[7], ratio_rep[8], ratio_unrep[8],
              ratio_rep[9], ratio_unrep[9])
rep_PPI_matrix[1,] #2       5      19      17      10      14      33      34      22 
rep_PPI_matrix[3,] #7      10      21      20      10      21      38      36      22 
unrep_PPI_matrix[1,]#69      44      43      37      37      33      27      42      35 
unrep_PPI_matrix[3,]# 132      89      57      41      41      38      34      42      36 

counts_label = c("2/7", "69/132", "5/10", "44/89", "19/21", "43/57",
                 "17/20", "37/41", "10/10", "37/41", "14/21", "33/38",
                 "33/38", "27/34", "34/36", "42/42", "22/22", "35/36")
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2C_Validation_bar_plot.pdf", width= 6, height=5)
barCenter = barplot(ratio_all*100, horiz=F, beside=F, ylim=c(0,100), ylab="Validation rate (%)",
                    space= c(0.4, 0.08, 0.4, 0.08, 0.4, 0.08, 0.4, 0.08, 0.4, 0.08,
                             0.4, 0.08, 0.4, 0.08, 0.4, 0.08, 0.4, 0.08),
                    col= apple_colors[c(5,3)] , axisnames=F, border=NA, cex.axis=0.8)
legend(-0.5,120, legend=c("Previously reported", "Previously unreported"),fill=apple_colors[c(5,3)], cex=0.8, bty="n",
       border=FALSE, xpd = TRUE)
text(x= barCenter, y = ratio_all*100 + 2, labels = counts_label, cex=0.5, xpd = TRUE)
env_num_loc = rep(0, 9)
for(i in 1:9){
  env_num_loc[i] = mean(barCenter[(2*i-1):(2*i)])
}
text(x = env_num_loc, y = -8, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

