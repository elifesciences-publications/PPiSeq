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
DMSO_fit = dataFrameReader_T("Paper_data/SD_mean_fitness_positive.csv")

#DMSO_fit_2 = dataFrameReader_T("Paper_data/SD_mean_fitness_positive.csv")

#DMSO_fit$SD[which(DMSO_fit$SD == 0)]= 1e-8
#mean_sd = mean(DMSO_fit$SD) # 0.06912734
#sd_sd = sd(DMSO_fit$SD) # 0.05069146
#error_sd = qnorm(0.975)*sd_sd/sqrt(nrow(DMSO_fit)) # 8.224907e-05
#error_sd = qnorm(0.975)*sd_sd/sqrt(4) # 0.04967672
#pdf("Working_figure/SFigures/Figure1_related/FigureSX_SD_standard_error_estimation_related_Figure1C/QQ_plot_sd_SD2.pdf")
#qqnorm(DMSO_fit$SD, frame = FALSE, pch = 1)
#qqline(DMSO_fit$SD, col = apple_colors[5], lwd = 2)
#dev.off()
#pdf("Working_figure/SFigures/Figure1_related/FigureSX_SD_standard_error_estimation_related_Figure1C/Histogram_sd_SD2.pdf")
#hist(DMSO_fit$SD, breaks = seq(0, 0.95, by = 0.01), xlab = "Standard error of fitness",
#     ylab = "Frequency", xlim = c(0, 0.4), main = "SD")
#dev.off()


#col_values = c("2" = "#fc8d62", "3" = "#8da0cb", "4" = "#e78ac3", 
#               "10" = apple_colors[7], "100" = apple_colors[4])

length(which(DMSO_fit$Positive == 0)) # 1443413
length(which(DMSO_fit$Positive != 0)) # 5286

#DMSO_fit$Positive[which(DMSO_fit$Positive == 0)] = "Negative"
#DMSO_fit$Positive[which(DMSO_fit$Positive == 1)] = "Positive"

library(scales)
#col_neg = alpha(apple_colors[5], 0.1) 
#col_pos = alpha(apple_colors[7], 0.5)
#col_values = c("Positive" = col_pos, "Negative" = col_neg)

bin_fit = seq(-0.2, 1, by = 0.1)
matrix_mean_CI = data.frame(bin_fit, rep(0,length(bin_fit)), rep(0,length(bin_fit)),rep(0,length(bin_fit)), rep(0, length(bin_fit)))
sd = DMSO_fit$SD[which(DMSO_fit$Mean_fitness <= -0.2 & DMSO_fit$Mean_fitness > -0.3)]
matrix_mean_CI[1,2] = mean(sd)
matrix_mean_CI[1,3] = median(sd)
sem_sd= sd(sd)/(length(sd)^0.5)
matrix_mean_CI[1,4] = mean(sd) - qnorm(0.975)*sem_sd
matrix_mean_CI[1,5] = mean(sd) + qnorm(0.975)*sem_sd
for(i in 2:length(bin_fit)){
  sd = DMSO_fit$SD[which(DMSO_fit$Mean_fitness <= bin_fit[i] & DMSO_fit$Mean_fitness > bin_fit[i-1])]
  matrix_mean_CI[i,2] = mean(sd)
  matrix_mean_CI[i,3] = median(sd)
  sem_sd= (sd(sd)/(length(sd))^0.5)
  matrix_mean_CI[i,4] = mean(sd) - qnorm(0.975)*sem_sd
  matrix_mean_CI[i,5] = mean(sd) + qnorm(0.975)*sem_sd
}
colnames(matrix_mean_CI) = c("Fitness", "Mean", "Median", "Lower", "Upper")
#colnames(matrix_CI) = c("Fitness", "Lower", "Upper")
#scale_y_continuous(name = "Standard deviation of fitness",
                   #limits=c(1e-8, 1),
                   #trans = log10_trans(),
                   #breaks= trans_breaks("log10", function(x) 10^x),
                   #labels = trans_format("log10", math_format(10^.x)))+

##### Make a similar matrix for the positive PPIs
DMSO_pos = DMSO_fit[which(DMSO_fit$Positive != 0),]
bin_fit = seq(-0.2, 1, by = 0.1)
matrix_mean_CI_pos = data.frame(bin_fit, rep(0,length(bin_fit)), rep(0,length(bin_fit)),rep(0,length(bin_fit)), rep(0, length(bin_fit)))
sd = DMSO_pos$SD[which(DMSO_pos$Mean_fitness <= -0.2 & DMSO_pos$Mean_fitness > -0.3)]
matrix_mean_CI_pos[1,2] = mean(sd)
matrix_mean_CI_pos[1,3] = median(sd)
sem_sd= sd(sd)/(length(sd)^0.5)
matrix_mean_CI_pos[1,4] = mean(sd) - qnorm(0.975)*sem_sd
matrix_mean_CI_pos[1,5] = mean(sd) + qnorm(0.975)*sem_sd
for(i in 2:length(bin_fit)){
        sd = DMSO_pos$SD[which(DMSO_pos$Mean_fitness <= bin_fit[i] & DMSO_pos$Mean_fitness > bin_fit[i-1])]
        matrix_mean_CI_pos[i,2] = mean(sd)
        matrix_mean_CI_pos[i,3] = median(sd)
        sem_sd= (sd(sd)/(length(sd))^0.5)
        matrix_mean_CI_pos[i,4] = mean(sd) - qnorm(0.975)*sem_sd
        matrix_mean_CI_pos[i,5] = mean(sd) + qnorm(0.975)*sem_sd
}
colnames(matrix_mean_CI_pos) = c("Fitness", "Mean", "Median", "Lower", "Upper")


library(ggplot2)
ggplot() +
        geom_hex(aes(x = Mean_fitness, y = SD, fill = log10(..count..)), DMSO_fit,
                 bins = 60, size = 0.05) +
        scale_fill_gradient(low= "white", high = apple_colors[5])+
        #scale_fill_gradientn(colours = c("white", apple_colors[10], apple_colors[3]))+
        #geom_smooth(aes(x = Mean_fitness, y = SD), DMSO_fit, method= 'lm',  
                    #size = 1, se =FALSE, col= apple_colors[6])+
        #geom_line(aes(x = Fitness, y = Mean), matrix_mean_CI, size = 0.2,col = apple_colors[11])+
        geom_line(aes(x = Fitness, y = Median), matrix_mean_CI,col = apple_colors[11],
                  linetype = "dashed")+
        geom_line(aes(x = Fitness, y = Median), matrix_mean_CI_pos,col = apple_colors[7],
                  linetype = "dashed")+
        #geom_ribbon(aes(x = Fitness, ymin= Lower, ymax = Upper), matrix_mean_CI, alpha = 0.2) +
    
        #geom_line(aes(x = Fitness, y = Upper), matrix_CI, linetype = "dashed", col = apple_colors[11])+
        #scale_color_manual(name = "",values = col_values) +
        scale_y_continuous(name = "Standard deviation of fitness",
                     limits=c(0, 0.5),
                     breaks= seq(0,0.5, by = 0.1),
                     labels =seq(0,0.5, by = 0.1)) +
  
        scale_x_continuous(name = "Mean fitness of each protein pair", 
                           limits=c(-0.4, 1.0),
                           breaks=seq(-0.4, 1.0, by =0.2),
                           labels = seq(-0.4, 1.0, by= 0.2)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.key=element_blank(), legend.position =c(0.9,0.8)) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("Working_figure/Figure1/Figure1D_standard_error_fitness_estimation_Reported_wide.pdf", width = 5, height = 5)
#ggsave("~/Desktop/Figure1F_standard_error_fitness_estimation_Reported_wide.pdf", width = 7, height = 5)
#ggsave("Working_figure/SFigures/Figure1_related/FigureSX_SD_standard_error_estimation_related_Figure1C/Hexagon_plot_standard_deviation_SD_mean.pdf", width = 7, height = 5)
