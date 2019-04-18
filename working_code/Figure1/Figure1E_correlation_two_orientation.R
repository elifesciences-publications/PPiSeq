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
######### Figure 1E: scatter plot to show the correlation between mean fitness values of the same PPI with two directions
setwd("~/Dropbox/PPiSeq_02/")
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv") # 1459163
DMSO_pos = DMSO_mean[which(DMSO_mean[,7] == "1"),] # 5211
# first find these PPIs that have two orientations in our data (check function of mark_duplicates_fast)
DMSO_pos_two = mark_duplicates_fast(DMSO_pos[,1]) # 4819  The second column are the opposite orientation (0 if none)
pos_PPI_two_direction = DMSO_pos_two[which(DMSO_pos_two[,2] != 0),] # 392
fitness_two= matrix(0, nrow(pos_PPI_two_direction), 2)
for(i in 1:nrow(fitness_two)){
  fitness_two[i,1] = DMSO_pos[which(DMSO_pos[,1] == pos_PPI_two_direction[i,1]),3]
  fitness_two[i,2] = DMSO_pos[which(DMSO_pos[,1] == pos_PPI_two_direction[i,2]),3]
}

#Transform the data into a data.frame
fit_01 = as.numeric(fitness_two[,1])
fit_02 = as.numeric(fitness_two[,2])
PPI = pos_PPI_two_direction[,1]
fitness_two_data = data.frame(PPI, fit_01, fit_02)
cor(fitness_two_data[,2], fitness_two_data[,3], method = "pearson") # 0.575
cor(fitness_two_data[,2], fitness_two_data[,3], method = "spearman") # 0.514

library(ggplot2)
ggplot() +
  geom_point(aes(x= fit_01, y= fit_02), col = apple_colors[7], alpha = 0.5) +
  geom_smooth(aes(x = seq(0.2, 1, by= 0.2), y = seq(0.2, 1, by = 0.2)), 
              method='lm',se = FALSE, color = apple_colors[11], cex = 0.4, linetype = 2) +
  scale_y_continuous(name = "Fitness of ORF2-DHFR[1,2] X ORF1-DHFR[3]",
                     limits=c(0.2, 1),
                     breaks=seq(0.2,1, by =0.2),
                     labels = seq(0.2,1, by= 0.2)) +
  scale_x_continuous(name = "Fitness of ORF1-DHFR[1,2] X ORF2-DHFR[3]", 
                     limits=c(0.2, 1),
                     breaks=seq(0.2,1, by =0.2),
                     labels = seq(0.2,1, by= 0.2)) +
  
  annotate("text", x = 0.4, y = 0.95, label = expression(paste("Spearman's ", italic(r), " = 0.51")), 
           parse = TRUE, size = 5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black", hjust =1),
        axis.text.y.left = element_text(size = 10, color = "black"))

ggsave("~/Dropbox/PPiSeq_02/working_figure/Figure1/Figure1E_correlation_fitness_two_directions.pdf", width = 4, height = 4)
