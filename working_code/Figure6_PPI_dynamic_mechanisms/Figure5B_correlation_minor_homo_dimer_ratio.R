
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

########## Correlation between homo-dimer correaltion with ratio of larger proteins
homo_plot_matrix = csvReader_T("Working_data/homo_dimer/PPI_homo_dimer_correlation_abundance.csv")
self_interacting = csvReader_T("Working_data/homo_dimer/Self_protein_matrix.csv") # 213
pos_number = rep(0, nrow(self_interacting))
pos_high_cor_number = rep(0, nrow(self_interacting))
pos_low_cor_number = rep(0, nrow(self_interacting))
neg_high_cor_number = rep(0, nrow(self_interacting))
neg_low_cor_number = rep(0, nrow(self_interacting))
for(i in 1:nrow(self_interacting)){
  index = which(homo_plot_matrix[,2] == self_interacting[i,1])
  pos_number[i] = length(which(as.numeric(homo_plot_matrix[index,6]) > 0))
  pos_high_cor_number[i] = length(which(as.numeric(homo_plot_matrix[index,6]) > 0 &
                                          as.numeric(homo_plot_matrix[index,3]) > 0))
  pos_low_cor_number[i] = length(which(as.numeric(homo_plot_matrix[index,6]) > 0 &
                                          as.numeric(homo_plot_matrix[index,3]) < 0))
  neg_high_cor_number[i] = length(which(as.numeric(homo_plot_matrix[index,6]) < 0 &
                                          as.numeric(homo_plot_matrix[index,3]) > 0))
  neg_low_cor_number[i] = length(which(as.numeric(homo_plot_matrix[index,6]) < 0 &
                                          as.numeric(homo_plot_matrix[index,3]) < 0))
}

self_interacting_large = cbind(self_interacting, pos_number, pos_high_cor_number, 
                               pos_low_cor_number, neg_high_cor_number, neg_low_cor_number)
ratio_pos = as.numeric(self_interacting_large[,8])/as.numeric(self_interacting_large[,2])
ratio_pos_high_cor = as.numeric(self_interacting_large[,9])/as.numeric(self_interacting_large[,2])
self_interacting_final = cbind(self_interacting_large, ratio_pos, ratio_pos_high_cor)
csvWriter(self_interacting_final, "Working_data/homo_dimer/Self_protein_matrix_comprehensive.csv")
self_interacting_final = csvReader_T("Working_data/homo_dimer/Self_protein_matrix_comprehensive.csv")

self_interacting_final= self_interacting_final[which(as.numeric(self_interacting_final[,2]) >= 10),] # 152
cor(as.numeric(self_interacting_final[,3]), as.numeric(self_interacting_final[,13]), method = "spearman") # 0.3495059
cor(as.numeric(self_interacting_final[,3]), as.numeric(self_interacting_final[,14]), method = "spearman") # 0.5911601
sum(as.numeric(self_interacting_final[,2]))
library(scales)
pdf("Working_figure/Figure5/Figure 5B Correlation_coefficient_versus_large_protein_ratio_>=10_PPI.pdf", width = 5, height =5)
plot(as.numeric(self_interacting_final[,3]), as.numeric(self_interacting_final[,13]), type = "p", bty = "n",
     col= alpha(apple_colors[5], 0.5), xlab = "Mean correlation coefficient with homo-dimer", pch = 16, xlim = c(-0.6, 0.8),
     ylim = c(0, 1), ylab = "Homo-dimer as the less abundant protein (ratio)")
text(-0.2, 0.8, expression(paste("Spearman's ", italic(r), " = 0.35")))
dev.off()

pdf("Working_figure/Figure5/Figure 5B Correlation_coefficient_versus_large_protein_high-cor_ratio_>=10_PPI.pdf", width = 5, height =5)
plot(as.numeric(self_interacting_final[,3]), as.numeric(self_interacting_final[,14]), type = "p", bty = "n",
     col= alpha(apple_colors[5], 0.5), xlab = "Mean correlation coefficient with homo-dimer", pch = 16, xlim = c(-0.6, 0.8),
     ylim = c(0, 1), ylab = "Homo-dimer as the less abundant protein (ratio)")
text(-0.2, 0.95, expression(paste("Spearman's ", italic(r), " = 0.59")))
dev.off()
