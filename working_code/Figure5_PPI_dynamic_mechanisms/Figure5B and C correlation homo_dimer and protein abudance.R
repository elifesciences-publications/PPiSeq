
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

### Function to extract homo-dimers
extract_self_PPI = function(DMSO_PPI){
  self_PPI = extract_repeat_PPI(as.character(DMSO_PPI[,1])) # 448
  DMSO_self_PPI = DMSO_PPI[which(DMSO_PPI[,1] %in% self_PPI),]
  return(DMSO_self_PPI)
}

#A function to extract specific PPIs
check_specific_protein = function(PPI, Gene_Carbon){
  PPI_chosen = "0"
  protein_pair = split_string_vector(PPI[,1])
  if(length(Gene_Carbon) > 1){
    for(i in 1:nrow(protein_pair)){
      if(protein_pair[i,1] %in% Gene_Carbon | protein_pair[i,2] %in% Gene_Carbon){
        PPI_chosen = c(PPI_chosen, PPI[i,1])
      }
    }  
  }else {
    for(i in 1:nrow(protein_pair)){
      if(protein_pair[i,1] == Gene_Carbon | protein_pair[i,2] == Gene_Carbon){
        PPI_chosen = c(PPI_chosen, PPI[i,1])
      }
    }  
  }
  
  PPI_chosen = PPI_chosen[2:length(PPI_chosen)]
  return(PPI_chosen)
}
##### Create a matrix to show the PPI, homo_dimer, correlation, protein_abundance, homo-dimer_abundance
setwd("~/Dropbox/PPiSeq_02/")
PPI_fit= csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
homo_dimer = extract_self_PPI(PPI_fit) # 213
homo_dimer_protein = split_string_vector(homo_dimer[,1])[,1] # 213
protein_abundance = csvReader_T("~/Dropbox/PPiseq_02/Working_data/protein_abundance/table_S4.csv")

homo_dimer_matrix = c("0", "0", 0, 0, 0) 
for(i in 1:length(homo_dimer_protein)){
  protein_chosen = homo_dimer_protein[i]
  PPI = check_specific_protein(PPI_fit, protein_chosen)
  if (length(PPI) > 1){
    PPI_fit_chosen = PPI_fit[which(PPI_fit[,1] %in% PPI),]
    PPI_self = extract_self_PPI(PPI_fit_chosen)
    PPI_fit_others = PPI_fit_chosen[which(!PPI_fit_chosen[,1] %in% PPI_self[1]),]
    cor = rep(0, nrow(PPI_fit_others))
    for(j in 1:nrow(PPI_fit_others)){
      cor[j] = cor(as.numeric(PPI_self[4:12]),as.numeric(PPI_fit_others[j,4:12]))
    }
    protein_others_matrix = split_string_vector(PPI_fit_others[,1])
    protein_others = rep("0", nrow(protein_others_matrix))
    for(i in 1:nrow(protein_others_matrix)){
      protein_pair = protein_others_matrix[i,]
      protein_others[i] = protein_pair[which(protein_pair != protein_chosen)]
    }
    protein_others_abundance = log2(as.numeric(protein_abundance[match(protein_others, protein_abundance[,1]),5]))
    if(protein_chosen %in% protein_abundance[,1]){
      homo_dimer_abundance = log2(as.numeric(protein_abundance[which(protein_abundance[,1] == protein_chosen),5]))
    }else{
      homo_dimer_abudance = NA
    }
    
    matrix_combine = cbind(protein_others, rep(protein_chosen, length(protein_others)), cor, 
                           protein_others_abundance, rep(homo_dimer_abundance, length(protein_others)))
    homo_dimer_matrix = rbind(homo_dimer_matrix, matrix_combine)
  }
  
}

homo_dimer_matrix = homo_dimer_matrix[2:nrow(homo_dimer_matrix),] # 7370
homo_dimer_matrix_omit = na.omit(homo_dimer_matrix) # 7312
abundance_change = as.numeric(homo_dimer_matrix_omit[,4]) - as.numeric(homo_dimer_matrix_omit[,5])
homo_dimer_matrix_final = cbind(homo_dimer_matrix_omit, abundance_change)
colnames(homo_dimer_matrix_final) = c("Protein", "Homo-dimer", "Correlation", 
                                      "Protein_amount", "Homo_amount", "Difference_amount")
csvWriter(homo_dimer_matrix_final, "Working_data/homo_dimer/PPI_homo_dimer_correlation_abundance.csv")

homo_plot_matrix = dataFrameReader_T("Working_data/homo_dimer/PPI_homo_dimer_correlation_abundance.csv")

self_interacting = dataFrameReader_T("Working_data/homo_dimer/Self_protein_matrix.csv")
self_interacting_large = as.character(self_interacting[which(self_interacting$PPI_number >5 
                                                             & abs(self_interacting$Mean_cor) >= 0.4
                                                          ),1])
#self_interacting_large = as.character(self_interacting[which(abs(self_interacting$Mean_cor) >= 0.2 
#& self_interacting$PPI_number >=20),1])
homo_plot_matrix_large = homo_plot_matrix[which(as.character(homo_plot_matrix[,2]) %in% self_interacting_large),] #2488
length(which(homo_plot_matrix_large$Correlation > 0 & homo_plot_matrix_large$Difference_amount > 0)) #1379
length(which(homo_plot_matrix_large$Correlation > 0 & homo_plot_matrix_large$Difference_amount < 0)) #812
length(which(homo_plot_matrix_large$Correlation < 0 & homo_plot_matrix_large$Difference_amount < 0)) #119
length(which(homo_plot_matrix_large$Correlation < 0 & homo_plot_matrix_large$Difference_amount > 0)) #136
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
ggsave("~/Desktop/Figure5B_abudance_change_versus_correlation_homo_dimer_>5_PPI_>0.4_dot.pdf", width = 7, height =5)

########## Correlation between homo-dimer correaltion with ratio of larger proteins
homo_plot_matrix = csvReader_T("Working_data/homo_dimer/PPI_homo_dimer_correlation_abundance.csv")
self_interacting = csvReader_T("Working_data/homo_dimer/Self_protein_matrix.csv") # 213
pos_number = rep(0, nrow(self_interacting))
for(i in 1:nrow(self_interacting)){
  index = which(homo_plot_matrix[,2] == self_interacting[i,1])
  pos_number[i] = length(which(as.numeric(homo_plot_matrix[index,6]) > 0))
}

self_interacting_large = cbind(self_interacting, pos_number)
ratio = as.numeric(self_interacting_large[,8])/as.numeric(self_interacting_large[,2])
self_interacting_final = cbind(self_interacting_large, ratio)
csvWriter(self_interacting_final, "Working_data/homo_dimer/Self_protein_matrix_comprehensive.csv")
self_interacting_final= self_interacting_final[which(as.numeric(self_interacting_final[,2]) >= 5),]
cor(as.numeric(self_interacting_final[,3]), as.numeric(self_interacting_final[,9])) # 0.35
library(scales)
pdf("~/Desktop/Correlation_coefficient_versus_large_protein_ratio_>5_PPI.pdf", width = 5, height =5)
plot(as.numeric(self_interacting_final[,3]), as.numeric(self_interacting_final[,9]), type = "p",
     col= alpha(apple_colors[5], 0.5), xlab = "Correlation coefficient with homo-dimer (>= 5 PPI)", pch = 16,
     ylab = "Ratio of PPIs with the homo-dimer as minor protein")
dev.off()

