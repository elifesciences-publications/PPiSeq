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

### Figure S2A_check the fitness distribution for PPIs that are detected in different environments
# I have generated normalized fitness values in each environment. 
# In another code, I get the positive PPIs out and put fitness values in different environments in the same row

setwd("~/Dropbox/PPiSeq_02/")
PPI_fit_norm = csvReader_T("Working_data/Positive_PPI_environment/Pos_PPI_normalized_fit.csv") # 14164
PPI_dup = mark_duplicates_fast(PPI_fit_norm[,1]) #13430
matrix = matrix(NA, nrow(PPI_dup), 2* ncol(PPI_fit_norm))
for(i in 1:nrow(PPI_dup)){
  if (PPI_dup[i,2] != "0"){
    matrix[i,1:10] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
    matrix[i,11:20] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,2]),]
  }else{
    matrix[i,1:10] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
  }
}
colnames(matrix) = c("PPI", "DMSO", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", 
                     "NaCl", "16C", "FK506", "PPI_opposite", "DMSO", "H2O2", "HU",
                     "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
mean_fitness = rep(0, nrow(matrix))
for(i in 1:length(mean_fitness)){
  mean_fitness[i] = mean(na.omit(as.numeric(matrix[i, c(2:10, 12:20)])))
}

matrix_final = cbind(matrix[,1], mean_fitness, matrix[,2:ncol(matrix)])
csvWriter(matrix_final, "Working_data/Positive_PPI_environment/Normalzied_fitness_PPI_all.csv")

## Make violin plot for PPIs that are detected in different number of environments
PPI_fit = csvReader_T("Working_data/Positive_PPI_environment/Normalzied_fitness_PPI_all.csv")
PPI_env_count = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
nrow(PPI_env_count)
nrow(PPI_fit)
mean_fit = as.numeric(PPI_fit[match(PPI_env_count[,1], PPI_fit[,1]),2])
PPI_env_fit = data.frame(PPI_env_count[,1], as.character(PPI_env_count[,2]), mean_fit)
colnames(PPI_env_fit) = c("PPI", "Environment", "Fitness")
PPI_env_fit$Environment = factor(PPI_env_fit$Environment, levels = as.character(1:9))
## Get the number of PPIs for each environment
table(PPI_env_fit$Environment)
# 1:8326; 2:1356; 3:655; 4:535; 5:555; 6:664; 7:570; 8:619; 9:513
### All PPIs
count = c("8326", "1356", "655", "535", "555", "664", "570", "619", "513")

library(ggplot2)

ggplot(PPI_env_fit, aes(x = Environment, y = Fitness, group = Environment))+
  #geom_boxplot(col = apple_colors[5])+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), col= apple_colors[5])+
  #geom_point(col = apple_colors[3])

  stat_summary(aes(x = Environment, y = Fitness, group = Environment), PPI_env_fit,
               fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
  scale_x_discrete(limits = c("1","2", "3","4", "5", "6", "7", "8", "9"))+
  scale_y_continuous(name = "Mean fitness of a PPI", 
                     limits=c(0, 1.2),
                     breaks = seq(0,1.2, by =0.2),
                     labels = seq(0,1.2, by= 0.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"))+
  xlab("Number of environments in which a PPI is observed") + 
  annotate("text", x = 1:9,  y = rep(0.1, 9), label = count)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2A_mean_fit_PPI_different_group/SFigure2_Mean_fitness_PPI_different_envrionments_violin.pdf", width= 5, height = 5)

### Tease out the reported PPIs
reported_PPI = csvReader_T("Working_data/multiple_validated_PPI.csv")
PPI_env_count_reported = match_both_direction(PPI_env_count, reported_PPI[,1]) # 1853
mean_fit = as.numeric(PPI_fit[match(PPI_env_count_reported[,1], PPI_fit[,1]),2]) # 1853
PPI_env_fit = data.frame(PPI_env_count_reported[,1], as.character(PPI_env_count_reported[,2]), mean_fit)
colnames(PPI_env_fit) = c("PPI", "Environment", "Fitness")
PPI_env_fit$Environment = factor(PPI_env_fit$Environment, levels = as.character(1:9))

table(PPI_env_fit$Environment)
# 1:308; 2:144; 3:112; 4:115; 5:148; 6:212; 7:254; 8:338; 9:222
### All PPIs
count = c("308", "144", "112", "115", "148", "212", "254", "338", "222")

### Reported PPIs
library(ggplot2)
ggplot(PPI_env_fit, aes(x = Environment, y = Fitness, group = Environment))+
  #geom_boxplot(col = apple_colors[3])+
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), col= apple_colors[5])+
  #geom_point(col = apple_colors[3])
  
  stat_summary(aes(x = Environment, y = Fitness, group = Environment), PPI_env_fit,
               fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
  scale_x_discrete(limits = c("1","2", "3","4", "5", "6", "7", "8", "9"))+
  scale_y_continuous(name = "Mean fitness of a PPI", 
                     limits=c(0, 1.2),
                     breaks = seq(0,1.2, by =0.2),
                     labels = seq(0,1.2, by= 0.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"))+
  xlab("Number of environments in which a PPI is observed") +
  annotate("text", x = 1:9,  y = rep(0.1, 9), label = count)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2A_mean_fit_PPI_different_group/SFigure2_Mean_fitness_PPI_different_envrionments_violin_reported.pdf", width= 5, height = 5)

