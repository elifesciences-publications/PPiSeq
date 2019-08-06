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

#### Use Psp2 containing PPIs as an example to explain cellular localization change underlying PPI dynamics
#### Focus on PPIs in which Psp2 as less abundant protein
setwd("~/Dropbox/PPiSeq_02/")
homo_plot_matrix = csvReader_T("Working_data/homo_dimer/PPI_homo_dimer_correlation_abundance.csv")

PPI_psp2_matrix = homo_plot_matrix[which(homo_plot_matrix[,2] == "YML017W" & as.numeric(homo_plot_matrix[,6]) > 0),]
PPI_psp2_matrix_order = PPI_psp2_matrix[order(as.numeric(PPI_psp2_matrix[,3]), decreasing = T),]
psp2_dimer = c("YML017W", "YML017W", 1, 12.29548, 12.29548, 0)
PPI_psp2_all = rbind(psp2_dimer, PPI_psp2_matrix_order)

GO_slim = as.matrix(read.table("Working_data/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
GO_slim_homo = GO_slim[which(GO_slim[,1] %in% PPI_psp2_all[,1]),]
GO_slim_homo_CC = GO_slim_homo[which(GO_slim_homo[,4] == "C"),] # 726
cc_unique = unique(GO_slim_homo_CC[,5]) # 10


### Make a heatmap to show the cellular localization of these proteins
protein_CC_matrix = matrix(0, nrow(PPI_psp2_all), (length(cc_unique) + 2))
protein_CC_matrix[,1] = PPI_psp2_all[,1]
protein_CC_matrix[,2] = as.numeric(PPI_psp2_all[,3])
colnames(protein_CC_matrix) = c("Protein", "Mean_cor",cc_unique)

for(i in 1:nrow(protein_CC_matrix)){
  index = which(GO_slim_homo_CC[,1] == protein_CC_matrix[i,1])
  CC_chosen = GO_slim_homo_CC[index,5]
  protein_CC_matrix[i,which(colnames(protein_CC_matrix) %in% CC_chosen)] = 1
}

csvWriter(protein_CC_matrix,"Working_data/homo_dimer/Psp2_PPIs_minor_localization.csv")
setwd("~/Dropbox/PPiSeq_02/")
psp2_CC_matrix = dataFrameReader_T("Working_data/homo_dimer/Psp2_PPIs_minor_localization.csv")

psp2_CC_heatmap = psp2_CC_matrix[,3:ncol(psp2_CC_matrix)]
colnames(psp2_CC_heatmap) = c(cc_unique)
rownames(psp2_CC_heatmap) = as.character(psp2_CC_matrix[,1])
#row_ann = data.frame(Correlation = as.character(psp2_CC_matrix$Mean_cor))
row_ann = data.frame(Correlation = psp2_CC_matrix$Mean_cor)

col_chosen = c(apple_colors[5], "white",apple_colors[3],apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=15)
my_colour = list(Correlation = color_scale)
row.names(row_ann) = row.names(psp2_CC_heatmap)
library(pheatmap)
col_chosen = c(apple_colors[c(9,3)])
#color_scale = colorRampPalette(col_chosen)(n=100)

fit_heatmap = pheatmap(psp2_CC_heatmap, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames=TRUE,
                       show_colnames=T, annotation_colors = my_colour, treeheight_col = 0, border_color = "white",
                       annotation_row = row_ann, fontsize_row = 5, fontsize = 6, fontsize_col = 6, col = col_chosen, 
                       legend = FALSE)

save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(fit_heatmap, "Working_figure/Figure5/Figure5E_Psp2_PPI_CC.pdf", width = 5, height =3)

###### Plot the PPI fitness across different environments
setwd("~/Dropbox/PPiSeq_02/")
homo_plot_matrix = csvReader_T("Working_data/homo_dimer/PPI_homo_dimer_correlation_abundance.csv")

PPI_psp2_matrix = homo_plot_matrix[which(homo_plot_matrix[,2] == "YML017W" & as.numeric(homo_plot_matrix[,6]) > 0),]
PPI_psp2_matrix_order = PPI_psp2_matrix[order(as.numeric(PPI_psp2_matrix[,3]), decreasing = T),]
psp2_dimer = c("YML017W", "YML017W", 1, 12.29548, 12.29548, 0)
PPI_psp2_all = rbind(psp2_dimer, PPI_psp2_matrix_order)

PPI_fit_csv = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")

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

PPI_psp2_fit = check_specific_protein(PPI_fit_csv, "YML017W")
PPI_psp2_fit_split = split_string_vector(PPI_psp2_fit)
other_protein = rep(0, length(PPI_psp2_fit))
for(i in 1:length(PPI_psp2_fit)){
  protein_pair = PPI_psp2_fit_split[i,]
  if(protein_pair[1] == protein_pair[2]){
    other_protein[i] = protein_pair[1]
  }else{
    other_protein[i] = protein_pair[which(protein_pair != "YML017W")]
  }
}

PPI_psp2_all_PPI = PPI_psp2_fit[match(PPI_psp2_all[,1], other_protein)]
PPI_psp2_all_PPI_fit = PPI_fit_csv[match(PPI_psp2_all_PPI, PPI_fit_csv[,1]),]

PPI_psp2_all_final = cbind(PPI_psp2_all, PPI_psp2_all_PPI_fit)
csvWriter(PPI_psp2_all_final, "Working_data/homo_dimer/Psp2_PPIs_fitness_across_environments.csv")

PPI_psp2_fit_heatmap = dataFrameReader_T("Working_data/homo_dimer/Psp2_PPIs_fitness_across_environments.csv")

col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=100)

psp2_fit_heatmap = PPI_psp2_fit_heatmap[,10:ncol(PPI_psp2_fit_heatmap)]
colnames(psp2_fit_heatmap) = c("SD", "H2O2", "Hydroxyurea", "Doxorubicin", "Forskolin",
                               "Raffinose", "NaCl",  "16 \u00B0C", "FK506")
rownames(psp2_fit_heatmap) = as.character(PPI_psp2_fit_heatmap$Protein)
#row_ann = data.frame(Correlation = as.character(psp2_CC_matrix$Mean_cor))
row_ann = data.frame(Correlation = PPI_psp2_fit_heatmap$Correlation)

col_chosen = c(apple_colors[5], "white",apple_colors[3],apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=15)
my_colour = list(Correlation = color_scale)
row.names(row_ann) = row.names(psp2_fit_heatmap)
library(pheatmap)
col_chosen_fit = c(apple_colors[5], "#e7d4e8",apple_colors[7])
col_scale_fit = colorRampPalette(col_chosen_fit)(n=40)

fit_heatmap = pheatmap(psp2_fit_heatmap, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames=TRUE,
                       show_colnames=T, annotation_colors = my_colour, treeheight_col = 0, border_color = "white",
                       annotation_row = row_ann, fontsize_row = 5, fontsize = 6, fontsize_col = 6, col = col_scale_fit
                      )
save_pheatmap_pdf(fit_heatmap, "Working_figure/Figure5/Figure5D_Psp2_PPI_fitness_environments.pdf", width = 5.5, height =3)
#save_pheatmap_pdf(fit_heatmap, "~/Desktop/Figure5G_Psp2_PPI_fitness_environments.pdf", width = 5, height =2)
