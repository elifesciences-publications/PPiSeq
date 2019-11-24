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

### Make barplot to show the percentage
### put the reported and unreported on to the same figure
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/")
rep_PPI_matrix = dataFrameReader_T("Reported_validation_matrix_SD_merge.csv")
unrep_PPI_matrix = dataFrameReader_T("Unreported_validation_matrix_SD_merge.csv")
ratio_rep = rep_PPI_matrix[4,-1]
ratio_unrep = unrep_PPI_matrix[4,-1]
ratio_all = as.numeric(c(ratio_rep[1], ratio_unrep[1], ratio_rep[2], ratio_unrep[2], 
                         ratio_rep[3], ratio_unrep[3], ratio_rep[4], ratio_unrep[4],
                         ratio_rep[5], ratio_unrep[5], ratio_rep[6], ratio_unrep[6],
                         ratio_rep[7], ratio_unrep[7], ratio_rep[8], ratio_unrep[8],
                         ratio_rep[9], ratio_unrep[9]))


rep_PPI_matrix[1,] #   0       0       5       7      14      15      23      43      16
rep_PPI_matrix[3,] #   1       1       6       9      19      18      27      44      16
unrep_PPI_matrix[1,]# 9      22      22      13      28      38      30      37      26
unrep_PPI_matrix[3,]#13      32      31      20      33      45      32      38      27
counts_label = c("0/1", "9/13", "0/1", "22/32", "5/6", "22/31",
                 "7/9", "13/20", "14/19", "28/33", "15/18", "38/45",
                 "23/27", "30/32", "43/44", "37/38", "16/16", "26/27")
library(RColorBrewer)
#col_chosen = brewer.pal(3,"Dark2")[1:2]
col_chosen = c("#d73027","#4575b4")
pdf("~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/FigureS4_validation_rate/FigureS4A_Validation_bar_plot_merge_reported_unreported.pdf", width= 5.5, height=5)
barCenter = barplot(ratio_all*100, horiz=F, beside=F, ylim=c(0,100), ylab="Validation rate (%)",
                    space= c(0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15,
                             0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15),
                    col= col_chosen , axisnames=F, border=NA, cex.axis=0.8)
legend(-0.5,120, legend=c("Previously reported", "Previously unreported"),fill=col_chosen, cex=0.8, bty="n",
       border=FALSE, xpd = TRUE)
text(x= barCenter, y = ratio_all*100 + 2, labels = counts_label, cex=0.5, xpd = TRUE)
env_num_loc = rep(0, 9)
for(i in 1:9){
  env_num_loc[i] = mean(barCenter[(2*i-1):(2*i)])
}
text(x = env_num_loc, y = -8, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()
