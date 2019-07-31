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

pdf("Working_Figure/Figure5/Figure5A_Homo_dimer_Number_PPI_versus_correlation.pdf", width= 5, height = 5)
plot(self_protein_matrix[,2], self_protein_matrix[,3], type = "p", cex =1,
     xlab = "Number of PPIs containing a self-interacting protein", col = apple_colors[5],
     ylab = "Mean correlation coefficient")
#text(as.numeric(self_protein_matrix[,2]), as.numeric(self_protein_matrix[,3]),
     #labels = self_protein_matrix[,1], cex=0.4)
dev.off()
mean_cor = as.numeric(self_protein_matrix[,3]) # min: -0.4212758, max: 0.7209423,
pdf("Working_Figure/Figure5/Figure5B_Histogram_correlation.pdf", width= 5, height = 5)
hist(as.numeric(self_protein_matrix[,3]), breaks = seq(-0.5, 0.75, by = 0.05), 
     xlab= "Mean correlation coefficient", main = NA, )
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

