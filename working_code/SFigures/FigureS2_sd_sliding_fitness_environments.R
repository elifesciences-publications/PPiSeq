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

### Estimate the standard error of the fitness for the same protein protein pair
setwd("~/Dropbox/PPiSeq_02/")

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

library(ggplot2)
plot_hexagon_sd = function(matrix_mean_CI, matrix_mean_CI_pos, output){
  ggplot() +
    #geom_hex(aes(x = Mean_fitness, y = sd, fill = log10(..count..)), DMSO_fit,
             #bins = 60, size = 0.05) +
    #scale_fill_gradient(low= "white", high = apple_colors[8])+
    geom_line(aes(x = Fitness, y = Median), matrix_mean_CI,col = apple_colors[5], size = 1.2)+
    geom_line(aes(x = Fitness, y = Median), matrix_mean_CI_pos,col = apple_colors[7], size =1.2)+
    #geom_ribbon(aes(x = Fitness, ymin= Lower, ymax = Upper), matrix_mean_CI, alpha = 0.2) +
    
    #geom_line(aes(x = Fitness, y = Upper), matrix_CI, linetype = "dashed", col = apple_colors[11])+
    #scale_color_manual(name = "",values = col_values) +
    scale_y_continuous(name = "Standard deviation of fitness",
                       limits=c(0, 0.3),
                       breaks= seq(0,0.3, by = 0.05),
                       labels =seq(0,0.3, by = 0.05)) +
    
    scale_x_continuous(name = "Estimated mean fitness of each protein pair", 
                       limits=c(-0.2, 1.0),
                       breaks=seq(-0.2, 1.0, by =0.2),
                       labels = seq(-0.2, 1.0, by= 0.2)) +
    #labs(fill = expression('Log'[10]* '(count)')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key=element_blank(), legend.position =c(0.8,0.7))+
  
    theme(axis.text.x = element_text(size = 7, color = "black"),
          axis.text.y.left = element_text(size = 7, color = "black"),
          axis.title=element_text(size=7))
    ggsave(output, width = 2.5, height = 2.5)
}


### SD2
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/SD2_mean_fitness_positive.csv")
matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/SD2_sd_sliding_fitness.pdf"
plot_hexagon_sd(matrix_mean_CI, matrix_mean_CI_pos, output)


### FK506
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/FK506_mean_fitness_positive.csv")
matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/FK506_sd_sliding_fitness.pdf"
plot_hexagon_sd(matrix_mean_CI, matrix_mean_CI_pos, output)

### H2O2
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/H2O2_mean_fitness_positive.csv")
matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/H2O2_sd_sliding_fitness.pdf"
plot_hexagon_sd(matrix_mean_CI, matrix_mean_CI_pos, output)

### Hydroxyurea
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/Hydroxyurea_mean_fitness_positive.csv")
matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/Hydroxyurea_sd_sliding_fitness.pdf"
ggplot() +
  #geom_hex(aes(x = Mean_fitness, y = sd, fill = log10(..count..)), DMSO_fit,
           #bins = 60, size = 0.05) +
  #scale_fill_gradient(low= "white", high = apple_colors[8])+
  geom_line(aes(x = Fitness, y = Median), matrix_mean_CI,col = apple_colors[5], size = 1.2)+
  geom_line(aes(x = Fitness, y = Median), matrix_mean_CI_pos,col = apple_colors[7], size =1.2)+
  #geom_ribbon(aes(x = Fitness, ymin= Lower, ymax = Upper), matrix_mean_CI, alpha = 0.2) +
  
  #geom_line(aes(x = Fitness, y = Upper), matrix_CI, linetype = "dashed", col = apple_colors[11])+
  #scale_color_manual(name = "",values = col_values) +
  scale_y_continuous(name = "Standard deviation of fitness",
                     limits=c(0, 0.3),
                     breaks= seq(0,0.3, by = 0.05),
                     labels =seq(0,0.3, by = 0.05)) +
  
  scale_x_continuous(name = "Estimated mean fitness of each protein pair", 
                     limits=c(-0.2, 0.8),
                     breaks=seq(-0.2, 0.8, by =0.2),
                     labels = seq(-0.2, 0.8, by= 0.2)) +
  labs(fill = expression('Log'[10]* '(count)')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank()) +
  theme(axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y.left = element_text(size = 7, color = "black"),
        axis.title=element_text(size=7))
ggsave(output, width = 2.5, height = 2.5)

### NaCl
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/NaCl_mean_fitness_positive.csv")
matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/NaCl_sd_sliding_fitness.pdf"
ggplot() +
  #geom_hex(aes(x = Mean_fitness, y = sd, fill = log10(..count..)), DMSO_fit,
           #bins = 60, size = 0.05) +
  #scale_fill_gradient(low= "white", high = apple_colors[8])+
  geom_line(aes(x = Fitness, y = Median), matrix_mean_CI,col = apple_colors[5], size = 1.2)+
  geom_line(aes(x = Fitness, y = Median), matrix_mean_CI_pos,col = apple_colors[7], size =1.2)+
  #geom_ribbon(aes(x = Fitness, ymin= Lower, ymax = Upper), matrix_mean_CI, alpha = 0.2) +
  
  #geom_line(aes(x = Fitness, y = Upper), matrix_CI, linetype = "dashed", col = apple_colors[11])+
  #scale_color_manual(name = "",values = col_values) +
  scale_y_continuous(name = "Standard deviation of fitness",
                     limits=c(0, 0.3),
                     breaks= seq(0,0.3, by = 0.05),
                     labels =seq(0,0.3, by = 0.05)) +
  
  scale_x_continuous(name = "Estimated mean fitness of each protein pair", 
                     limits=c(-0.2, 0.5),
                     breaks=seq(-0.2, 0.5, by =0.2),
                     labels = seq(-0.2, 0.5, by= 0.2)) +
  labs(fill = expression('Log'[10]* '(count)')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank()) +
  theme(axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y.left = element_text(size = 7, color = "black"),
        axis.title=element_text(size=7))
ggsave(output, width = 2.5, height = 2.5)

### Forskolin
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/Forskolin_mean_fitness_positive.csv")
Create_sliding_matrix_Forskolin = function(DMSO_fit){
  #### Create a matrix that contains the mean, median, and CI of sd in each bin of fitness
  bin_fit = seq(-0.15, 1, by = 0.05)
  matrix_mean_CI = data.frame(bin_fit, rep(0,length(bin_fit)), rep(0,length(bin_fit)),rep(0,length(bin_fit)), rep(0, length(bin_fit)))
  sd = DMSO_fit$sd[which(DMSO_fit$Mean_fitness <= -0.15 & DMSO_fit$Mean_fitness > -0.2)]
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

matrix_mean_CI = Create_sliding_matrix_Forskolin(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix_Forskolin(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/Forskolin_sd_sliding_fitness.pdf"
plot_hexagon_sd(matrix_mean_CI, matrix_mean_CI_pos, output)

### Raffinose
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/Raffinose_mean_fitness_positive.csv")
matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/Raffinose_sd_sliding_fitness.pdf"
#plot_hexagon_sd(matrix_mean_CI, matrix_mean_CI_pos, output)
ggplot() +
        #geom_hex(aes(x = Mean_fitness, y = sd, fill = log10(..count..)), DMSO_fit,
                # bins = 60, size = 0.05) +
        #scale_fill_gradient(low= "white", high = apple_colors[8])+
        geom_line(aes(x = Fitness, y = Median), matrix_mean_CI,col = apple_colors[5], size = 1.2)+
        geom_line(aes(x = Fitness, y = Median), matrix_mean_CI_pos,col = apple_colors[7], size =1.2)+
        #geom_ribbon(aes(x = Fitness, ymin= Lower, ymax = Upper), matrix_mean_CI, alpha = 0.2) +
        
        #geom_line(aes(x = Fitness, y = Upper), matrix_CI, linetype = "dashed", col = apple_colors[11])+
        #scale_color_manual(name = "",values = col_values) +
        scale_y_continuous(name = "Standard deviation of fitness",
                           limits=c(0, 0.5),
                           breaks= seq(0,0.5, by = 0.1),
                           labels =seq(0,0.5, by = 0.1)) +
        
        scale_x_continuous(name = "Estimated mean fitness of each protein pair", 
                           limits=c(-0.2, 1),
                           breaks=seq(-0.2, 1, by =0.2),
                           labels = seq(-0.2, 1, by= 0.2)) +
        labs(fill = expression('Log'[10]* '(count)')) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.key=element_blank()) +
        theme(axis.text.x = element_text(size = 7, color = "black"),
              axis.text.y.left = element_text(size = 7, color = "black"),
              axis.title=element_text(size=7))
ggsave(output, width = 2.5, height = 2.5)


### Doxorubicin
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/Doxorubicin_mean_fitness_positive.csv")
matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/Doxorubicin_sd_sliding_fitness.pdf"
ggplot() +
  #geom_hex(aes(x = Mean_fitness, y = sd, fill = log10(..count..)), DMSO_fit,
           #bins = 60, size = 0.05) +
  #scale_fill_gradient(low= "white", high = apple_colors[8])+
  geom_line(aes(x = Fitness, y = Median), matrix_mean_CI,col = apple_colors[5], size = 1.2)+
  geom_line(aes(x = Fitness, y = Median), matrix_mean_CI_pos,col = apple_colors[7], size =1.2)+
  #geom_ribbon(aes(x = Fitness, ymin= Lower, ymax = Upper), matrix_mean_CI, alpha = 0.2) +
  
  #geom_line(aes(x = Fitness, y = Upper), matrix_CI, linetype = "dashed", col = apple_colors[11])+
  #scale_color_manual(name = "",values = col_values) +
  scale_y_continuous(name = "Standard deviation of fitness",
                     limits=c(0, 0.3),
                     breaks= seq(0,0.3, by = 0.05),
                     labels =seq(0,0.3, by = 0.05)) +
  
  scale_x_continuous(name = "Estimated mean fitness of each protein pair", 
                     limits=c(-0.2, 0.6),
                     breaks=seq(-0.2, 0.6, by =0.2),
                     labels = seq(-0.2, 0.6, by= 0.2)) +
  labs(fill = expression('Log'[10]* '(count)')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank()) +
  theme(axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y.left = element_text(size = 7, color = "black"),
        axis.title=element_text(size=7))
ggsave(output, width = 2.5, height = 2.5)


### 16C
DMSO_fit = dataFrameReader_T("Paper_data/PPI_mean_fitness_calling_files/Cold_16C_mean_fitness_positive.csv")
matrix_mean_CI = Create_sliding_matrix(DMSO_fit)
DMSO_fit_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
matrix_mean_CI_pos = Create_sliding_matrix(DMSO_fit_pos)
output = "Working_figure/SFigures/paper/FigureS2_sd_sliding_fit_environments/Cold_16C_sd_sliding_fitness.pdf"

ggplot() +
  #geom_hex(aes(x = Mean_fitness, y = sd, fill = log10(..count..)), DMSO_fit,
           #bins = 60, size = 0.05) +
  #scale_fill_gradient(low= "white", high = apple_colors[8])+
  geom_line(aes(x = Fitness, y = Median), matrix_mean_CI,col = apple_colors[5], size = 1.2)+
  geom_line(aes(x = Fitness, y = Median), matrix_mean_CI_pos,col = apple_colors[7], size =1.2)+
  #geom_ribbon(aes(x = Fitness, ymin= Lower, ymax = Upper), matrix_mean_CI, alpha = 0.2) +
  
  #geom_line(aes(x = Fitness, y = Upper), matrix_CI, linetype = "dashed", col = apple_colors[11])+
  #scale_color_manual(name = "",values = col_values) +
  scale_y_continuous(name = "Standard deviation of fitness",
                     limits=c(0, 0.3),
                     breaks= seq(0,0.3, by = 0.05),
                     labels =seq(0,0.3, by = 0.05)) +
  
  scale_x_continuous(name = "Estimated mean fitness of each protein pair", 
                     limits=c(-0.2, 0.3),
                     breaks=seq(-0.2, 0.3, by =0.1),
                     labels = seq(-0.2, 0.3, by= 0.1)) +
  labs(fill = expression('Log'[10]* '(count)')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank()) +
  theme(axis.text.x = element_text(size = 7, color = "black"),
        axis.text.y.left = element_text(size = 7, color = "black"),
        axis.title=element_text(size=7))
ggsave(output, width = 2.5, height = 2.5)




