### Estimate the standard error of the fitness for the same protein protein pair
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

DMSO_fit = dataFrameReader_T("PPI_mean_fitness_calling_files/SD_mean_fitness_positive.csv")
library(scales)
Create_sliding_matrix = function(DMSO_fit){
  #### Create a matrix that contains the mean, median, and CI of sd in each bin of fitness
  bin_fit = seq(-0.2, 1, by = 0.05)
  matrix_mean_CI = data.frame(bin_fit, rep(0,length(bin_fit)), rep(0,length(bin_fit)),rep(0,length(bin_fit)), rep(0, length(bin_fit)))
  sd = DMSO_fit$sd[which(DMSO_fit$Mean_fitness <= -0.2 & DMSO_fit$Mean_fitness > -0.25)]
  matrix_mean_CI[1,2] = mean(sd)
  matrix_mean_CI[1,3] = median(sd)
  sem_sd= sd(sd)/(length(sd)^0.5)
  matrix_mean_CI[1,4] = mean(sd) - qnorm(0.975)*sem_sd
  matrix_mean_CI[1,5] = mean(sd) + qnorm(0.975)*sem_sd
  for(i in 2:length(bin_fit)){
    sd = DMSO_fit$sd[which(DMSO_fit$Mean_fitness <= bin_fit[i] & DMSO_fit$Mean_fitness > bin_fit[i-1])]
    matrix_mean_CI[i,2] = mean(sd)
    matrix_mean_CI[i,3] = median(sd)
    sem_sd= (sd(sd)/(length(sd))^0.5)
    matrix_mean_CI[i,4] = mean(sd) - qnorm(0.975)*sem_sd
    matrix_mean_CI[i,5] = mean(sd) + qnorm(0.975)*sem_sd
  }
  colnames(matrix_mean_CI) = c("Fitness", "Mean", "Median", "Lower", "Upper")
  return(matrix_mean_CI)
}

matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)

library(ggplot2)
ggplot() +
        #geom_hex(aes(x = Mean_fitness, y = sd, fill = log10(..count..)), DMSO_fit,
                 #bins = 60, size = 0.05) +
        #scale_fill_gradient(low= "white", high = apple_colors[8])+
        geom_line(aes(x = Fitness, y = Median), matrix_mean_CI,col = apple_colors[5], size = 1.2)+
        geom_line(aes(x = Fitness, y = Median), matrix_mean_CI_pos,col = apple_colors[7], size =1.2)+
        #geom_ribbon(aes(x = Fitness, ymin= Lower, ymax = Upper), matrix_mean_CI, alpha = 0.2) +
    
        scale_y_continuous(name = "Standard deviation of fitness",
                     limits=c(0, 0.12),
                     breaks= seq(0,0.12, by = 0.02),
                     labels =seq(0,0.12, by = 0.02)) +
  
        scale_x_continuous(name = "Estimated mean fitness of each protein pair", 
                           limits=c(-0.2, 1.0),
                           breaks=seq(-0.2, 1.0, by =0.2),
                           labels = seq(-0.2, 1.0, by= 0.2)) +
        labs(fill = expression('Log'[10]* '(count)')) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.key=element_blank(), legend.position =c(0.8,0.7)) +
        theme(axis.text.x = element_text(size = 8, color = "black"),
              axis.text.y.left = element_text(size = 8, color = "black"),
              axis.title=element_text(size=8))
ggsave("Figures/Figure1/Figure1C_SD_standard_error_fitness.pdf", width = 2.8, height = 2)
