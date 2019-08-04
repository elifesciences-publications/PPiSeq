########################### Extract all the homo-dimers and then check their correlation with PPIs containing that Protein
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

setwd("~/Dropbox/PPiSeq_02/")
PPI_fit_csv = dataFrameReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
# A function to extract self-interacting PPIs
extract_self_PPI = function(DMSO_PPI){
  self_PPI = extract_repeat_PPI(as.character(DMSO_PPI[,1])) # 448
  DMSO_self_PPI = DMSO_PPI[which(DMSO_PPI[,1] %in% self_PPI),]
  return(DMSO_self_PPI)
}
homo_dimer = extract_self_PPI(PPI_fit_csv) # 213

homo_dimer_protein = split_string_vector(homo_dimer[,1])[,1]

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
### Construct a matrix to show the correlation between homo_dimers and PPIs containing the self-interacting PPIs
PPI_fit_csv = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
# A function to extract self-interacting PPIs
self_protein_matrix = matrix(0, length(homo_dimer_protein), 7)
self_protein_matrix[,1] = homo_dimer_protein
for(i in 1:length(homo_dimer_protein)){
  PPI = check_specific_protein(PPI_fit_csv, homo_dimer_protein[i])
  if (length(PPI) == 1){
    self_protein_matrix[i,2:7] = 0
  }else{
    PPI_fit = PPI_fit_csv[which(PPI_fit_csv[,1] %in% PPI),]
    PPI_self = extract_self_PPI(PPI_fit)
    PPI_fit_others = PPI_fit[which(!PPI_fit[,1] %in% PPI_self[1]),]
    self_protein_matrix[i,2] = nrow(PPI_fit_others)
    cor = rep(0, nrow(PPI_fit_others))
    for(j in 1:nrow(PPI_fit_others)){
      cor[j] = cor(as.numeric(PPI_self[4:12]),as.numeric(PPI_fit_others[j,4:12]))
    }
    self_protein_matrix[i,3] = mean(cor)
    self_protein_matrix[i,4] = median(cor)
    self_protein_matrix[i,5] = sd(cor)
    self_protein_matrix[i,6] = min(cor)
    self_protein_matrix[i,7] = max(cor)
  }
  
}
colnames(self_protein_matrix) = c("Protein", "PPI_number", "Mean_cor", "Median_cor",
                                  "SD_cor", "Min_cor", "Max_cor")
csvWriter(self_protein_matrix, "Working_data/homo_dimer/Self_protein_matrix.csv")

self_protein_matrix = dataFrameReader_T("Working_data/homo_dimer/Self_protein_matrix.csv")
pdf("Working_Figure/Figure5/Figure5A_Homo_dimer_Number_PPI_versus_correlation.pdf", width= 5, height = 5)
plot(self_protein_matrix[,2], self_protein_matrix[,3], type = "p", cex =1,
     xlab = "Number of PPIs containing a self-interacting protein", col = apple_colors[5],
     ylab = "Mean correlation coefficient")
#text(as.numeric(self_protein_matrix[,2]), as.numeric(self_protein_matrix[,3]),
     #labels = self_protein_matrix[,1], cex=0.4)
dev.off()

pdf("Working_Figure/Figure5/Homo_dimer_Number_PPI_versus_sd.pdf", width= 5, height = 5)
plot(self_protein_matrix[,2], self_protein_matrix[,5], type = "p", cex =1, 
     xlab = "Number of PPIs containing a self-interacting protein", col = apple_colors[5],
     ylab = "Standard deviation of correlation coefficient")
#text(as.numeric(self_protein_matrix[,2]), as.numeric(self_protein_matrix[,3]),
#labels = self_protein_matrix[,1], cex=0.4)
dev.off()

self_protein_matrix[which(self_protein_matrix$SD_cor > 0.7),]

pdf("Working_Figure/Figure5/Homo_dimer_mean_correlation_versus_sd.pdf", width= 5, height = 5)
plot(self_protein_matrix[,3], self_protein_matrix[,5], type = "p", cex =1,
     xlab = "Mean correlation coefficient", col = apple_colors[5],
     ylab = "Standard deviation of correlation coefficient")
#text(as.numeric(self_protein_matrix[,2]), as.numeric(self_protein_matrix[,3]),
#labels = self_protein_matrix[,1], cex=0.4)
dev.off()



mean_cor = as.numeric(self_protein_matrix[,3]) # min: -0.4212758, max: 0.7209423,
pdf("Working_Figure/Figure5/Figure5B_Histogram_correlation.pdf", width= 5, height = 5)
hist(as.numeric(self_protein_matrix[,3]), breaks = seq(-0.5, 0.75, by = 0.05), 
     xlab= "Mean correlation coefficient", main = NA)
dev.off()
#### Heatmap to show the fitness of homo_dimer across different environments
Cor = self_protein_matrix[,3]
Cor[which(Cor <= -0.4)] = -0.4
Cor[which(Cor <= -0.2 & Cor > -0.4)] = -0.2
Cor[which(Cor > -0.2 & Cor <= 0)] = -0.1
Cor[which(Cor > 0 & Cor < 0.2)] = 0.1
Cor[which(Cor >= 0.2 & Cor < 0.4)] = 0.2
Cor[which(Cor >= 0.4)] = 0.4

homo_dimer = data.frame(homo_dimer, Cor)
homo_dimer_heatmap = homo_dimer[,c(4,8,12,10,9,6,5,7,11)]
colnames(homo_dimer_heatmap) = c("SD", "Forskolin", "FK506", "NaCl", "Raffinose", "Hydroxyurea",  
                                 "H2O2", "Doxorubicin", "16 \u00B0C")
rownames(homo_dimer_heatmap) = as.character(homo_dimer[,1])
row_ann = data.frame(Correlation = as.character(homo_dimer$Cor))
my_colour = list(Correlation = c("-0.4" = "#ffff33","-0.2" = "#ff7f00", "-0.1" = "#984ea3",
                                 "0.1" = "#4daf4a", "0.2"="#377eb8", "0.4" = "#e41a1c"))
row.names(row_ann) = row.names(homo_dimer_heatmap)

hc = hclust(dist(homo_dimer_heatmap), method = "complete")
homo_dimer_heatmap_order = homo_dimer_heatmap[hc$order,]
col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=100)
library(pheatmap)
fit_heatmap = pheatmap(homo_dimer_heatmap_order, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames=FALSE,
                       annotation_colors = my_colour, annotation_row = row_ann, 
                       show_colnames=T, col = color_scale, fontsize_row =5)

save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(fit_heatmap, "Working_Figure/Figure5/Figure5C_homo_dimer_heatmap.pdf")

#### Make heatmap for YML017W
PPI_fit_csv = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
PPI_PSP2 = check_specific_protein(PPI_fit_csv, "YML017W")
PPI_fit_data = dataFrameReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
PPI_PSP2_fit = PPI_fit_data[which(PPI_fit_data[,1] %in% PPI_PSP2),]
which(PPI_PSP2_fit[,1] == "YML017W_YML017W")

PPI_PSP2_heatmap = PPI_PSP2_fit[,c(4,8,12,10,9,6,5,7,11)]
colnames(PPI_PSP2_heatmap) = c("SD", "Forskolin", "FK506", "NaCl", "Raffinose", "Hydroxyurea",  
                                      "H2O2", "Doxorubicin", "16 \u00B0C")
rownames(PPI_PSP2_heatmap) = PPI_PSP2_fit[,1]
hc_others = hclust(dist(PPI_PSP2_heatmap), method = "complete")
PPI_PSP2_heatmap_order = PPI_PSP2_heatmap[hc_others$order,]

col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=100)

library(pheatmap)
fit_heatmap = pheatmap(PPI_PSP2_heatmap_order, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames=TRUE,
                       show_colnames=T, col = color_scale, fontsize_row =5)

save_pheatmap_pdf <- function(x, filename, width=6, height=5) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#save_pheatmap_pdf(fit_heatmap, "Working_figure/Figure3/Figure3B_carbonhydrate_transport_fitness_environment_primary.pdf")
save_pheatmap_pdf(fit_heatmap, "Working_Figure/Figure5/Figure5D_PPI_PSP2_heatmap.pdf")

##############################
#check distribution of correlation coefficients between each homo-dimer and homo-dimer containing PPIs
setwd("~/Dropbox/PPiSeq_02/")
PPI_fit_csv = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
# A function to extract self-interacting PPIs
extract_self_PPI = function(DMSO_PPI){
  self_PPI = extract_repeat_PPI(as.character(DMSO_PPI[,1])) # 448
  DMSO_self_PPI = DMSO_PPI[which(DMSO_PPI[,1] %in% self_PPI),]
  return(DMSO_self_PPI)
}
homo_dimer = extract_self_PPI(PPI_fit_csv) # 213

homo_dimer_protein = split_string_vector(homo_dimer[,1])[,1]

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
#### Make distribution of correlation coefficients between each homodimer and homodimer containing PPIs
setwd("~/Dropbox/PPiSeq_02/working_data/homo_dimer/distribution_correlation_20/")
for(i in 1:length(homo_dimer_protein)){
  PPI = check_specific_protein(PPI_fit_csv, homo_dimer_protein[i])
  if (length(PPI) >= 20){
    PPI_fit = PPI_fit_csv[which(PPI_fit_csv[,1] %in% PPI),]
    PPI_self = extract_self_PPI(PPI_fit)
    PPI_fit_others = PPI_fit[which(!PPI_fit[,1] %in% PPI_self[1]),]
    cor = rep(0, nrow(PPI_fit_others))
    for(j in 1:nrow(PPI_fit_others)){
      cor[j] = cor(as.numeric(PPI_self[4:12]),as.numeric(PPI_fit_others[j,4:12]))
    }
    max = max(cor)
    min = min(cor)
    pdf(paste(homo_dimer_protein[i], "correlation coefficient distribution.pdf", sep = " "), width= 5, height=5)
    hist(cor, breaks = seq((min -0.1),(max + 0.1), by = 0.05), xlab = "Correlation coefficient with homo-dimer",
         ylab = "Frequency", xlim = c((min - 0.1), (max + 0.1)), main = homo_dimer_protein[i])
    dev.off()
  }
  
}

### check the homodimer's localization
setwd("~/Dropbox/PPiSeq_02/")
self_protein_matrix = csvReader_T("Working_data/homo_dimer/Self_protein_matrix.csv")
self_protein_matrix = self_protein_matrix[which(as.numeric(self_protein_matrix[,2]) >=20),]
self_protein_matrix = self_protein_matrix[order(as.numeric(self_protein_matrix[,3]), decreasing = T),]
GO_slim = as.matrix(read.table("Working_data/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
GO_slim_homo = GO_slim[which(GO_slim[,1] %in% self_protein_matrix[,1]),]
GO_slim_homo_CC = GO_slim_homo[which(GO_slim_homo[,4] == "C"),] # 726
cc_unique = unique(GO_slim_homo_CC[,5]) # 23

homo_CC_matrix = matrix(0, nrow(self_protein_matrix), (length(cc_unique) + 2))
homo_CC_matrix[,1] = self_protein_matrix[,1]
homo_CC_matrix[,2] = as.numeric(self_protein_matrix[,3])
colnames(homo_CC_matrix) = c("Protein", "Mean_cor",cc_unique)

for(i in 1:nrow(homo_CC_matrix)){
  index = which(GO_slim_homo_CC[,1] == homo_CC_matrix[i,1])
  CC_chosen = GO_slim_homo_CC[index,5]
  homo_CC_matrix[i,which(colnames(homo_CC_matrix) %in% CC_chosen)] = 1
}
homo_CC_matrix[1,]

csvWriter(homo_CC_matrix,"Working_data/homo_dimer/homo_dimer_>=20_localization.csv")

homo_CC = dataFrameReader_T("Working_data/homo_dimer/homo_dimer_>=20_localization.csv")
homo_CC_heatmap = homo_CC[,2:ncol(homo_CC)]
colnames(homo_CC_heatmap) = c("Cor",cc_unique)
rownames(homo_CC_heatmap) = as.character(homo_CC[,1])

col_chosen = c(apple_colors[5], "white",apple_colors[3],apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=100)

library(pheatmap)
fit_heatmap = pheatmap(homo_CC_heatmap, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames=TRUE,
                       show_colnames=T, fontsize_row = 3,col = color_scale, border_color = apple_colors[9])

save_pheatmap_pdf <- function(x, filename, width=6, height=7) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(fit_heatmap, "Working_data/homo_dimer/homo-dimer_>=20_localization.pdf")
#save_pheatmap_pdf(fit_heatmap, "~/Desktop/homo-dimer_>=20_localization.pdf")

############
##Grab PPIs that contain YML017W and show the correlation and localization
setwd("~/Dropbox/PPiSeq_02/")
PPI_fit_csv = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
GO_slim = as.matrix(read.table("Working_data/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
Check_CC_homodimer= function(protein_chosen, csv_name, pdf_name){
  PPI = check_specific_protein(PPI_fit_csv, protein_chosen)
  PPI_fit = PPI_fit_csv[which(PPI_fit_csv[,1] %in% PPI),]
  PPI_self = extract_self_PPI(PPI_fit)
  PPI_fit_others = PPI_fit[which(!PPI_fit[,1] %in% PPI_self[1]),]
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
  PPI_chosen_matrix = cbind(protein_others, cor)
  PPI_chosen_matrix = rbind(c(protein_chosen, 1), PPI_chosen_matrix)
  PPI_chosen_matrix= PPI_chosen_matrix[order(as.numeric(PPI_chosen_matrix[,2]), decreasing = T),]
  
  protein_CC_matrix = matrix(0, nrow(PPI_chosen_matrix), (length(cc_unique) + 2))
  protein_CC_matrix[,1] = PPI_chosen_matrix[,1]
  protein_CC_matrix[,2] = as.numeric(PPI_chosen_matrix[,2])
  colnames(protein_CC_matrix) = c("Protein", "Mean_cor",cc_unique)
  
  GO_slim_CC = GO_slim[which(GO_slim[,4] == "C"),] 
  
  for(i in 1:nrow(protein_CC_matrix)){
    index = which(GO_slim_CC[,1] == protein_CC_matrix[i,1])
    CC_chosen = GO_slim_CC[index,5]
    protein_CC_matrix[i,which(colnames(protein_CC_matrix) %in% CC_chosen)] = 1
  }
  
  csvWriter(protein_CC_matrix,csv_name)
  
  homo_CC = dataFrameReader_T(csv_name)
  homo_CC_heatmap = homo_CC[,2:ncol(homo_CC)]
  colnames(homo_CC_heatmap) = c("Cor",cc_unique)
  rownames(homo_CC_heatmap) = as.character(homo_CC[,1])
  
  col_chosen = c(apple_colors[5], "white",apple_colors[3],apple_colors[7])
  color_scale = colorRampPalette(col_chosen)(n=nrow(homo_CC))
  
  library(pheatmap)
  fit_heatmap = pheatmap(homo_CC_heatmap, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames=TRUE,
                         show_colnames=T, fontsize_row = 5,col = color_scale, border_color = apple_colors[9])
  
  save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
    pdf(filename, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  save_pheatmap_pdf(fit_heatmap, heatmap_name )
  #save_pheatmap_pdf(fit_heatmap, "~/Desktop/homo-dimer_>=20_localization.pdf")
  
}
csv_name = "Working_data/homo_dimer/localization/PPI_YDL222C_localization.csv"
heatmap_name = "Working_data/homo_dimer/localization/PPI_YDL222C_correlation_localization.pdf"
Check_CC_homodimer("YDL222C", csv_name, heatmap_name)

### check the abundance level in these proteins
setwd("~/Dropbox/PPiseq_02/")
PSP2_matrix = csvReader_T("Working_data/homo_dimer/localization/PPI_YDL222C_localization.csv")
protein_abundance = csvReader_T("~/Dropbox/PPiseq_02/Working_data/protein_abundance/table_S4.csv")
PSP2_matrix[1,]
protein_abundance[1,]
PSP2_abundance = log2(as.numeric(protein_abundance[match(PSP2_matrix[,1], protein_abundance[,1]),5]))

PSP2_cor_abundance = cbind(PSP2_matrix[,1:2], PSP2_abundance)
#pdf("Working_data/homo_dimer/correlation_abundance_homo_dimer_correlation_PSP2.pdf", width = 5, height = 5)
pdf("~/Desktop/correlation_abundance_homo_dimer_correlation_YDL222C.pdf", width = 5, height = 5)
plot(as.numeric(PSP2_cor_abundance[,2]), as.numeric(PSP2_cor_abundance[,3]),main = "YCR093W (0.668476868)",
     xlab = "Correlation coefficient with homo-dimer", ylab = "Median molecules per cell (Log2)")
lines(c(-1, -0.5, 0, 0.5, 1), rep(as.numeric(PSP2_cor_abundance[1,3]), 5), col = "red")
dev.off()
