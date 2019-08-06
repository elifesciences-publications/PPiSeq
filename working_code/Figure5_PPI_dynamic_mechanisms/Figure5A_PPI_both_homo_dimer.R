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

setwd("~/Dropbox/PPiSeq_02/")
PPI_fit= csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
homo_dimer = extract_self_PPI(PPI_fit) # 213
homo_dimer_protein = split_string_vector(homo_dimer[,1])[,1] # 213
homo_dimer[,1] = homo_dimer_protein
PPI_fit_non = PPI_fit[which(!PPI_fit[,1] %in% homo_dimer[,1]),] #12120

protein_abundance = csvReader_T("~/Dropbox/PPiseq_02/Working_data/protein_abundance/table_S4.csv")
protein_median_abundance = log2(as.numeric(protein_abundance[,5]))
PPI_fit_both = rep("0", ncol(PPI_fit))
for(i in 1:nrow(PPI_fit_non)){
  PPP = split_string(PPI_fit_non[i,1])
  if(PPP[1] %in% homo_dimer_protein & PPP[2] %in% homo_dimer_protein){
    PPI_fit_both = rbind(PPI_fit_both, PPI_fit_non[i,])
  }
}
PPI_fit_both= PPI_fit_both[2:nrow(PPI_fit_both),]
nrow(PPI_fit_both)

#Create a matrix to show the correlation and protein abundance
PPI_fit_both_cor = matrix(0, nrow(PPI_fit_both), 7)
PPI_fit_both_cor[,1] = PPI_fit_both[,1]
for(i in 1:nrow(PPI_fit_both_cor)){
  PPP = split_string(PPI_fit_both_cor[i,1])
  fit_PPP = as.numeric(PPI_fit_both[i,4:12])
  fit_dimer_1 = as.numeric(homo_dimer[which(homo_dimer[,1] == PPP[1]), 4:12])
  fit_dimer_2 = as.numeric(homo_dimer[which(homo_dimer[,1] == PPP[2]), 4:12])
  PPI_fit_both_cor[i,2] = cor(fit_PPP, fit_dimer_1, method = "spearman")
  PPI_fit_both_cor[i,3] = cor(fit_PPP, fit_dimer_2, method = "spearman")
  PPI_fit_both_cor[i,4] = protein_median_abundance[which(protein_abundance[,1] == PPP[1])]
  PPI_fit_both_cor[i,5] = protein_median_abundance[which(protein_abundance[,1] == PPP[2])]
  if(PPI_fit_both_cor[i,4] < PPI_fit_both_cor[i,5]){
    PPI_fit_both_cor[i,6] = "Minor"
    PPI_fit_both_cor[i,7] = "Major"
  }else if(PPI_fit_both_cor[i,4] == PPI_fit_both_cor[i,5]){
    PPI_fit_both_cor[i,6] = "Equal"
    PPI_fit_both_cor[i,7] = "Equal"
  }else{
    PPI_fit_both_cor[i,6] = "Major"
    PPI_fit_both_cor[i,7] = "Minor"
  }
}
csvWriter(PPI_fit_both_cor, "Working_data/homo_dimer/PPI_both_homo_dimer_cor.csv")

setwd("~/Dropbox/PPiSeq_02/")
PPI_fit_both_cor = csvReader_T("Working_data/homo_dimer/PPI_both_homo_dimer_cor.csv")
PPI_fit_both_cor_omit = na.omit(PPI_fit_both_cor) # 1351

Minor_cor_1 = PPI_fit_both_cor_omit[which(PPI_fit_both_cor_omit[,6] == "Minor"),2] # 752
Minor_cor_2 = PPI_fit_both_cor_omit[which(PPI_fit_both_cor_omit[,7] == "Minor"),3] # 598
minor_all = c(Minor_cor_1, Minor_cor_2)
minor_matrix = data.frame(as.numeric(minor_all), rep("Less abundant", length(minor_all)))
colnames(minor_matrix) = c("Correlation", "Label")

Major_cor_1 = PPI_fit_both_cor_omit[which(PPI_fit_both_cor_omit[,6] == "Major"),2] # 598
Major_cor_2 = PPI_fit_both_cor_omit[which(PPI_fit_both_cor_omit[,7] == "Major"),3] # 752
major_all = c(Major_cor_1, Major_cor_2)
major_matrix = data.frame(as.numeric(major_all), rep("More abundant", length(major_all)))
colnames(major_matrix) = c("Correlation", "Label")

t.test(as.numeric(minor_all), as.numeric(major_all)) # p-value = 1.8e-06

Equal_cor_1= PPI_fit_both_cor_omit[which(PPI_fit_both_cor_omit[,6] == "Equal"),2] # 1
Equal_cor_2= PPI_fit_both_cor_omit[which(PPI_fit_both_cor_omit[,6] == "Equal"),3] # 1

plot_matrix = rbind(minor_matrix, major_matrix)

library("ggplot2")
ggplot()+
  geom_violin(aes(x = Label, y = Correlation), plot_matrix, 
              draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE)  +
  #geom_boxplot(aes(x = PPI, y = fitness, group = PPI, col = color), bar_plot_data_control, 
  #show.legend = FALSE)+
  
  #scale_x_discrete(limits = c("ORF x Null","DHFR(-)", "LCL3 x SSM4","SSS1 x YOP1", 
                              #"CKA1 x CKB1", "HNM1 x KEX1", "GNP1 x SND3", "VOA1 x VPH1", 
                              #"DHFR(+)")) +
  
  stat_summary(aes(x = Label, y = Correlation), plot_matrix,
               fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1) +
  #scale_color_manual(name = "", breaks = c("Negative PPI", 'Positive PPI'),
                     #values  = apple_colors[c(5,7)]) +
  #scale_fill_manual(name = "", breaks = c( "Negative PPI", 'Positive PPI'),
                    #values  = apple_colors[c(5,7)])+
  scale_y_continuous(name = "Correlation coefficient between each PPI and homo-dimer", 
                     limits=c(-1, 1),
                     breaks = seq(-1,1, by =0.2),
                     labels = seq(-1,1, by= 0.2))+
  #guides(color = guide_legend(override.aes = list(size = 2, alpha = 0.5)))+
  theme(legend.key = element_blank(), legend.position = c(0.8,0.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.title.x = element_blank(),axis.text.y.left = element_text(size = 10, color = "black")) + 
  theme(text = element_text(size=12))+ theme(plot.margin = unit(c(0.2,0.1,0.4,0.5), "cm"))
ggsave("Working_figure/Figure5/Figure5A_PPI_both_homo_dimer_cor.pdf", width= 4.5, height = 4.5)
