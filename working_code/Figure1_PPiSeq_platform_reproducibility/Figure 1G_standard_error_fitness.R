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
DMSO_fit = dataFrameReader_T("Paper_data/DMSO_mean_fitness_positive.csv")
mean_sd = mean(DMSO_fit$SD) # 0.06912734
sd_sd = sd(DMSO_fit$SD) # 0.05069146
#error_sd = qnorm(0.975)*sd_sd/sqrt(nrow(DMSO_fit)) # 8.224907e-05
error_sd = qnorm(0.975)*sd_sd/sqrt(4) # 0.04967672
pdf("Working_figure/Figure1/Figure1G_standard_error_estimation/QQ_plot_sd_DMSO.pdf")
qqnorm(DMSO_fit$SD, frame = FALSE, pch = 1)
qqline(DMSO_fit$SD, col = apple_colors[5], lwd = 2)
dev.off()
pdf("Working_figure/Figure1/Figure1G_standard_error_estimation/Histogram_sd_DMSO.pdf")
hist(DMSO_fit$SD, breaks = seq(0, 0.95, by = 0.01), xlab = "Standard error of fitness",
     ylab = "Frequency", xlim = c(0, 0.4), main = "SD")
dev.off()

CI_sd_lower = mean_sd - error_sd # 0.01945062
CI_sd_upper = mean_sd + error_sd # 0.11880419

#col_values = c("2" = "#fc8d62", "3" = "#8da0cb", "4" = "#e78ac3", 
               #"10" = apple_colors[7], "100" = apple_colors[4])

DMSO_fit$Positive[which(DMSO_fit$Positive == 0)] = "Negative"
DMSO_fit$Positive[which(DMSO_fit$Positive == 1)] = "Positive"
col_values = c("Positive" = apple_colors[7], "Negative" = apple_colors[5])

matrix_CI = data.frame(seq(-0.8, 1.2, by = 0.2),rep(CI_sd_lower,11), rep(CI_sd_upper, 11))
colnames(matrix_CI) = c("Fitness", "Lower", "Upper")
library(ggplot2)
ggplot() +
        geom_point(aes(x = Mean_fitness, y = SD, col = as.character(Positive)), DMSO_fit, alpha = 0.1) +
        geom_line(aes(x = Fitness, y = Lower), matrix_CI, linetype = "dashed", col = apple_colors[11])+
        geom_line(aes(x = Fitness, y = Upper), matrix_CI, linetype = "dashed", col = apple_colors[11])+
        scale_color_manual(name = "",values = col_values) +
        scale_y_continuous(name = "Standard deviation of fitness",
                           limits=c(0, 1),
                           breaks=seq(0,1, by =0.1),
                           labels = seq(0,1, by= 0.1)) +
        scale_x_continuous(name = "Mean fitness of each protein pair", 
                           limits=c(-0.8, 1.2),
                           breaks=seq(-0.8, 1.2, by =0.2),
                           labels = seq(-0.8, 1.2, by= 0.2)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.key=element_blank()) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))
#ggsave("~/Desktop/Figure1F_standard_error_fitness_estimation_Reported_wide.pdf", width = 7, height = 5)
ggsave("Working_figure/Figure1/Figure1G_standard_error_estimation/Figure1G_standard_deviation_SD.pdf", width = 7, height = 5)
