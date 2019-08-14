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
all_PPI_matrix_final = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary_combine_SD.csv")
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
matrix_PPI_env_rep[3,] # 351 159 114 121 167 216 244 344 211
all_PPI_count = matrix_PPI_env_rep[1,] + matrix_PPI_env_rep[2,] + matrix_PPI_env_rep[3,]
all_PPI_count # 9465 1685  816  656  661  715  585  650  423
ratio_BioGrid = matrix_PPI_env_rep[3,]/all_PPI_count
ratio_BioGrid # 0.03708399 0.09436202 0.13970588 0.18445122 0.25264750 0.30209790 0.41709402 0.52923077 0.49881797
ratio_BioGrid_reported = c("3.7%", "9.4%", "14%", "18.4%", "25.3%", "30.2%", "41.7%", "52.9%", "49.9%")
ratio_PCA_low = matrix_PPI_env_rep[2,]/all_PPI_count
ratio_PCA_low # 0.03222398 0.09139466 0.16911765 0.23170732 0.26021180 0.33426573 0.35897436 0.36923077 0.39007092
ratio_PCA_low_overlapped = c("3.2%", "9.1%", "16.9%", "23.2%", "26%", "33.4%", "35.9%", "36.9%", "39%")

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

