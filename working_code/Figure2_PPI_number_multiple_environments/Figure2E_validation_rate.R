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
                  d_5_5, d_5_6, d_5_7) # 695
reported_PPI = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/multiple_validated_PPI.csv")
all_Tecan_rep = match_both_direction(all_Tecan, reported_PPI[,1]) # 185
all_Tecan_unrep = all_Tecan[which(!all_Tecan[,1] %in% all_Tecan_rep[,1]),] # 510
#split validated PPIs into different groups
PPI_group = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
length(which(all_Tecan[,1] %in% PPI_group[,1])) # 646 all in

rep_PPI_matrix = matrix(0, 4,10)
unrep_PPI_matrix = matrix(0, 4, 10)

for(i in 1:10){
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
                             "Envir_6", "Envir_7", "Envir_8", "Envir_9", "Envir_10")
rownames(rep_PPI_matrix) = c("Validate", "Nonvalidate", "Total", "Val_ratio")
csvWriter_rownames(rep_PPI_matrix, "Reported_validation_matrix.csv")

colnames(unrep_PPI_matrix) = c("Envir_1", "Envir_2", "Envir_3", "Envir_4", "Envir_5",
                               "Envir_6", "Envir_7", "Envir_8", "Envir_9", "Envir_10")
rownames(unrep_PPI_matrix) = c("Validate", "Nonvalidate", "Total", "Val_ratio")
csvWriter_rownames(unrep_PPI_matrix, "Unreported_validation_matrix.csv")

### Make barplot to show the percentage
### put the reported and unreported on to the same figure
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/")
rep_PPI_matrix = dataFrameReader_T("Reported_validation_matrix.csv")
unrep_PPI_matrix = dataFrameReader_T("Unreported_validation_matrix.csv")
ratio_rep = rep_PPI_matrix[4,-1]
ratio_unrep = unrep_PPI_matrix[4,-1]
ratio_all = as.numeric(c(ratio_rep[1], ratio_unrep[1], ratio_rep[2], ratio_unrep[2], 
              ratio_rep[3], ratio_unrep[3], ratio_rep[4], ratio_unrep[4],
              ratio_rep[5], ratio_unrep[5], ratio_rep[6], ratio_unrep[6],
              ratio_rep[7], ratio_unrep[7], ratio_rep[8], ratio_unrep[8],
              ratio_rep[9], ratio_unrep[9], ratio_rep[10], ratio_unrep[10]))
rep_PPI_matrix[1,] # 2       3       9      16      15       7      18      29      42       14
rep_PPI_matrix[3,] # 5       7      12      18      17       7      27      32      44       14
unrep_PPI_matrix[1,]#35      28      31      45      30      35      33      33      40       27
unrep_PPI_matrix[3,]#65      54      54      65      37      43      38      39      40       28
counts_label = c("2/5", "35/65", "3/17", "28/54", "9/12", "31/54",
                 "16/18", "45/65", "15/17", "30/37", "7/7", "35/43",
                 "18/27", "33/38", "29/32", "33/39", "42/44", "40/40",
                 "14/14", "27/28")
library(RColorBrewer)
#col_chosen = brewer.pal(3,"Dark2")[1:2]
col_chosen = c("#d73027","#4575b4")
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2E_Validation_bar_plot.pdf", width= 6, height=5)
barCenter = barplot(ratio_all*100, horiz=F, beside=F, ylim=c(0,100), ylab="Validation rate (%)",
                    space= c(0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15,
                             0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15),
                    col= col_chosen , axisnames=F, border=NA, cex.axis=0.8)
legend(-0.5,120, legend=c("Previously reported", "Previously unreported"),fill=col_chosen, cex=0.8, bty="n",
       border=FALSE, xpd = TRUE)
text(x= barCenter, y = ratio_all*100 + 2, labels = counts_label, cex=0.5, xpd = TRUE)
env_num_loc = rep(0, 9)
for(i in 1:10){
  env_num_loc[i] = mean(barCenter[(2*i-1):(2*i)])
}
text(x = env_num_loc, y = -8, labels = as.character(1:10), xpd = TRUE)
text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

