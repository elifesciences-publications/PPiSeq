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


#############
setwd("~/Dropbox/PPiSeq_02/")
homo_plot_matrix = dataFrameReader_T("Working_data/homo_dimer/PPI_homo_dimer_correlation_abundance.csv")

self_interacting = dataFrameReader_T("Working_data/homo_dimer/Self_protein_matrix.csv")
sum(self_interacting[which(self_interacting$PPI_number >= 5),2])
length(self_interacting[which(self_interacting$PPI_number >= 5),2]) # 197

self_interacting_large = as.character(self_interacting[which(self_interacting$PPI_number >= 10 
                                                             & abs(self_interacting$Mean_cor) >= 0.5
                                                          ),1])

#self_interacting_large = as.character(self_interacting[which(abs(self_interacting$Mean_cor) >= 0.2 
#& self_interacting$PPI_number >=20),1])
homo_plot_matrix_large = homo_plot_matrix[which(as.character(homo_plot_matrix[,2]) %in% self_interacting_large),] #2488

homo_plot_matrix_chosen = homo_plot_matrix[which(as.character(homo_plot_matrix[,2]) == "YML017W"),]

length(which(homo_plot_matrix_large$Correlation > 0 & homo_plot_matrix_large$Difference_amount > 0)) #779
length(which(homo_plot_matrix_large$Correlation > 0 & homo_plot_matrix_large$Difference_amount < 0)) #346
length(which(homo_plot_matrix_large$Correlation < 0 & homo_plot_matrix_large$Difference_amount < 0)) #37
length(which(homo_plot_matrix_large$Correlation < 0 & homo_plot_matrix_large$Difference_amount > 0)) #61
diff = homo_plot_matrix$Difference_amount # min: -10.45546; max: 5.863211
CC = homo_plot_matrix$Correlation # -0.9259513; max:0.9903858
library(scales)
col_chosen = alpha(apple_colors[5], 0.1)
#pdf("Working_figure/Figure5/Figure5B_abudance_change_versus_correlation_homo_dimer.pdf", width = 5, height =5)
#pdf("~/Desktop/Figure5B_abudance_change_versus_correlation_homo_dimer.pdf", width = 5, height =5)
#plot(homo_plot_matrix$Difference_amount, homo_plot_matrix$Correlation, xlim = c(-6, 6), ylim = c(-1,1),
     #type = "p", col = col_chosen, pch = 16)
#dev.off()
library(ggplot2)
ggplot() +
  #geom_hex(aes(x = Difference_amount, y = Correlation), homo_plot_matrix_large, size = 0.05, bins = 60)+
  #scale_fill_gradientn(colours = c(apple_colors[5], apple_colors[10], apple_colors[7]))+
  geom_point(aes(x = Difference_amount, y = Correlation), homo_plot_matrix_large, alpha = 0.5, col = apple_colors[5])+
  geom_vline(xintercept = 0, col = apple_colors[7])+
  geom_hline(yintercept = 0, col = apple_colors[7])+
  scale_y_continuous(name = "Correlation coefficient with homo-dimer",
                     limits=c(-1, 1),
                     breaks= seq(-1,1, by = 0.5),
                     labels =seq(-1,1, by = 0.5)) +
  
  scale_x_continuous(name = "Difference of protein abundance with homo-dmer", 
                     limits=c(-6, 6),
                     breaks=seq(-6, 6, by =2),
                     labels = seq(-6, 6, by= 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank()) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("Working_figure/Figure5/FigureS5C_abundance_change_versus_correlation_homo_dimer_>10_PPI_>0.5_dot.pdf", width = 6, height =5)

a = homo_plot_matrix_chosen$Correlation
b = homo_plot_matrix_chosen$Difference_amount
ggplot() +
  #geom_hex(aes(x = Difference_amount, y = Correlation), homo_plot_matrix_large, size = 0.05, bins = 60)+
  #scale_fill_gradientn(colours = c(apple_colors[5], apple_colors[10], apple_colors[7]))+
  geom_point(aes(x = Difference_amount, y = Correlation), homo_plot_matrix_chosen, alpha = 0.5, col = apple_colors[5])+
  geom_vline(xintercept = 0, col = apple_colors[10])+
  geom_hline(yintercept = 0, col = apple_colors[10])+
  scale_y_continuous(name = "Correlation coefficient with homo-dimer",
                     limits=c(-1, 1),
                     breaks= seq(-1,1, by = 0.5),
                     labels =seq(-1,1, by = 0.5)) +
  
  scale_x_continuous(name = "Difference of abundance from homo-dimer (Log2)", 
                     limits=c(-6, 6),
                     breaks=seq(-6, 6, by =2),
                     labels = seq(-6, 6, by= 2)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key=element_blank()) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("Working_figure/Figure5/Figure5C_abundance_change_versus_correlation_homo_dimer_Psp2_dot.pdf", width = 4.5, height =4.5)


