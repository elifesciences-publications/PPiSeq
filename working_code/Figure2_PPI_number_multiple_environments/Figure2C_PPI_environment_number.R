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

##################################################
# Figure2B make a pie plot to show the number of detected environments for each PPI 

# Or make a barplot to show how many of them have been reproted
setwd("~/Dropbox/PPiSeq_02/")
all_PPI_matrix_final = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
reported_PPI = csvReader_T("Working_data/multiple_validated_PPI.csv")
PCA_lower = as.matrix(read.table("Working_data/Tarassov_PPI_PPV_80.txt", header= T, sep = "\t"))
min(as.numeric(PCA_lower[,ncol(PCA_lower)])) # 80.60
PCA_lower_PPI = paste(PCA_lower[,1], PCA_lower[,4], sep = "_")
PCA_lower_PPI_matrix = cbind(PCA_lower_PPI, rep(1, length(PCA_lower_PPI)))
PCA_lower_PPI_reported = match_both_direction(PCA_lower_PPI_matrix, reported_PPI[,1]) # 3392
PCA_lower_PPI_unreported = PCA_lower_PPI_matrix[which(!PCA_lower_PPI_matrix[,1] %in% PCA_lower_PPI_reported[,1]),] # 6838

matrix_PPI_env_rep = matrix(0, 3, 9)
for(i in 1:9){
  all = all_PPI_matrix_final[which(as.numeric(all_PPI_matrix_final[,2]) == i),]
  all_reported_PCA_low = match_both_direction(all, PCA_lower_PPI_unreported[,1])
  all_reported_BioGrid = match_both_direction(all,reported_PPI[,1])
  all_reported = c(all_reported_PCA_low[,1], all_reported_BioGrid[,1])
  all_unreported = all[which(!all[,1] %in% all_reported),]
  matrix_PPI_env_rep[1,i] = nrow(all_unreported)
  matrix_PPI_env_rep[2,i] = nrow(all_reported_PCA_low)
  matrix_PPI_env_rep[3,i] = nrow(all_reported_BioGrid)
}
matrix_PPI_env_rep[3,] # 312 157 121 127 185 272 314 222 181
all_PPI_count = matrix_PPI_env_rep[1,] + matrix_PPI_env_rep[2,] + matrix_PPI_env_rep[3,]
all_PPI_count # 7812 1395  760  698  688  884  632  445  336
ratio_BioGrid = matrix_PPI_env_rep[3,]/all_PPI_count
ratio_BioGrid # 0.03993856 0.11254480 0.15921053 0.18194842 0.26889535 0.30769231 0.49683544 0.49887640 0.53869048
ratio_BioGrid_reported = c("4.0%", "11.3%", "15.9%", "18.2%", "26.9%", "30.8%", "49.7%", "49.9%", "53.9%")
ratio_PCA_low = matrix_PPI_env_rep[2,]/all_PPI_count
ratio_PCA_low # 0.03609831 0.11469534 0.16315789 0.24785100 0.28633721 0.34728507 0.33860759 0.40224719 0.35416667
ratio_PCA_low_overlapped = c("3.6%", "11.5%", "16.3%", "24.8%", "28.6%", "34.7%", "33.9%", "40.2%", "35.4%")

#col_purple = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
library(RColorBrewer)
#col_chosen = brewer.pal(4, "Set1")[2:4]
col_chosen = c("#4575b4","#fdae61","#d73027")
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2C_Number_environments_PPI_reproted.pdf", height = 5, width = 5)
#pdf("~/Desktop/Figure2B_Number_environments_PPI_reproted.pdf", height = 5, width = 5)
par(mar = c(3,4,2,1))
barCenter = barplot(matrix_PPI_env_rep, horiz=F, beside=F, ylim=c(0,10000), ylab="Number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
                    col= col_chosen, axisnames=F, border=NA)
legend("topright", legend=c("BioGRID", "Marginal PCA", "Previously unreported"), 
       fill=col_chosen[c(3,2,1)], bty="n", border=FALSE)
text(x= barCenter, y = all_PPI_count + 150, labels = ratio_PCA_low_overlapped, 
     cex=0.7, xpd = TRUE, col= col_chosen[2]) 
text(x= barCenter, y = all_PPI_count + 400, labels = ratio_BioGrid_reported, 
     cex=0.7, xpd = TRUE, col= col_chosen[3]) 
text(x= barCenter, y = -300, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -800, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()






## Or make a pie plot for environment_number
environment_count = data.frame(table(environment_number))
## Environment number per positive PPI

library(ggplot2)
colfunc <- colorRampPalette(apple_colors[c(5,3,7)])
ggplot(environment_count, aes(x = factor(1), y = Freq, fill= environment_number)) +
  geom_bar(width =1, size = 0.1, color = "white", stat = "identity") + 
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 1.8)+
  coord_polar("y") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values=colfunc(9),
                    name = "Number of environments",
                    breaks = 1:9,
                    labels = 1:9)

ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2B_PPI_environment_distribution.pdf", width =5 , height = 5)

