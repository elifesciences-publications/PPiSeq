########################### Make a heatmap to show the normalized fitness values for Carbohydrate transport PPIs
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

# Input the normalized fitness values for all PPIs in each environment
setwd("~/Dropbox/PPiSeq_02/")
PPI_fit_csv = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
# A function to extract specific PPIs
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

PPI_HXT5 = check_specific_protein(PPI_fit_csv, "YHR096C") # 45

PPI_fit = dataFrameReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
PPI_HXT5_fit = PPI_fit_csv[which(PPI_fit_csv[,1] %in% PPI_HXT5),]
protein_HXT5 = unique(c(split_string_vector(PPI_HXT5)[,1],split_string_vector(PPI_HXT5)[,2]))
protein_others= protein_HXT5[which(protein_HXT5 != "YHR096C")]#44
PPI_HXT5_others = check_specific_protein(PPI_fit_csv, protein_others) # 3044
PPI_others_fit = PPI_fit_csv[which(PPI_fit_csv[,1] %in% PPI_HXT5_others),]
which(PPI_HXT5_fit[,1] == "YHR096C_YHR096C") # 45
PPI_HXT5_final = rbind(PPI_HXT5_fit[45,], PPI_HXT5_fit[1:44,])
group_HXT5 = rep("First", nrow(PPI_HXT5_final))
cor_HXT5 = rep(0, nrow(PPI_HXT5_final))
for(i in 1:nrow(PPI_HXT5_final)){
  cor_HXT5[i] = cor(as.numeric(PPI_HXT5_final[1,4:12]), as.numeric(PPI_HXT5_final[i,4:12]))
}
cor_others = rep(0, nrow(PPI_others_fit))
for(i in 1:nrow(PPI_others_fit)){
  cor_others[i] = cor(as.numeric(PPI_HXT5_final[1,4:12]), as.numeric(PPI_others_fit[i,4:12]))
}

min(cor_others) # -0.9182119
max(cor_others) # 0.9903858
mean(cor_others) #0.1071557
min(cor_HXT5) # -0.6948885
max(cor_HXT5) # 1
mean(cor_HXT5) # 0.7271435
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure5/Histogram_HXT5_PPIs_different_groups.pdf", width = 5, height = 5)
hist(cor_others, breaks = seq(-1, 1, by = 0.05), xlab = "Correlation coefficient with HXT5-HXT5",
     ylab = "Frequency", main = NA, bty = "n", col = rgb(0,0,1,alpha = 0.5))
hist(cor_HXT5, breaks = seq(-1, 1, by = 0.05), col = rgb(1,0,0, alpha = 0.5), add = T)
legend("topleft", c( "Second", "First"), fill = c(rgb(0,0,1,alpha = 0.5),rgb(1,0,0, alpha = 0.5)),
       bty="n")
dev.off()

t.test(cor_HXT5, cor_others)


PPI_HXT5_heatmap = PPI_HXT5_final[,c(4,8,12,10,9,6,5,7,11)]
colnames(PPI_HXT5_heatmap) = c("SD", "Forskolin", "FK506", "NaCl", "Raffinose", "Hydroxyurea",  
                                 "H2O2", "Doxorubicin", "16 \u00B0C")
rownames(PPI_HXT5_heatmap) = as.character(PPI_HXT5_final[,1])
#csvWriter(PPI_carbon_heatmap, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/environment/carbonhydrate_transport_network/PPI_carbohydrate_transport_heatmap.csv")

hc = hclust(dist(PPI_HXT5_heatmap), method = "complete")
PPI_HXT5_heatmap_order = PPI_HXT5_heatmap[hc$order,]
col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=100)

library(pheatmap)
fit_heatmap = pheatmap(PPI_HXT5_heatmap_order, cluster_rows = FALSE, cluster_cols = TRUE, show_rownames=TRUE,
                       show_colnames=T, col = color_scale, fontsize_row =5)

save_pheatmap_pdf <- function(x, filename, width=6, height=5) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
#save_pheatmap_pdf(fit_heatmap, "Working_figure/Figure3/Figure3B_carbonhydrate_transport_fitness_environment_primary.pdf")
save_pheatmap_pdf(fit_heatmap, "Working_Figure/Figure5/PPI_HXT5_heatmap.pdf")

##### Futher reduce the PPI number by take the first 12 colums
PPI_HXT5_fit_chosen = PPI_HXT5_heatmap_order[1:12,]
PPI_HXT5_chosen = rownames(PPI_HXT5_fit_chosen)
protein_HXT5_chosen = unique(c(split_string_vector(PPI_HXT5_chosen)[,1], split_string_vector(PPI_HXT5_chosen)[,2]))
protein_others= protein_HXT5[which(protein_HXT5_chosen != "YHR096C")]#11
PPI_HXT5_others = check_specific_protein(PPI_fit_csv, protein_others) # 857
PPI_HXT5_others_fit = PPI_fit[which(PPI_fit[,1] %in% PPI_HXT5_others),]

PPI_HXT5_others_heatmap = PPI_HXT5_others_fit[,c(4,8,12,10,9,6,5,7,11)]
colnames(PPI_HXT5_others_heatmap) = c("SD", "Forskolin", "FK506", "NaCl", "Raffinose", "Hydroxyurea",  
                               "H2O2", "Doxorubicin", "16 \u00B0C")
rownames(PPI_HXT5_others_heatmap) = PPI_HXT5_others_fit[,1]
hc_others = hclust(dist(PPI_HXT5_others_heatmap), method = "complete")
PPI_HXT5_others_heatmap_order = PPI_HXT5_others_heatmap[hc_others$order,]

Group_1 = rep("First", nrow(PPI_HXT5_heatmap_order))
Group_2 = rep("Second", nrow(PPI_HXT5_others_heatmap_order))


PPI_HXT5_all_heatmap = rbind(PPI_HXT5_heatmap_order, PPI_HXT5_others_heatmap_order)

row_ann = data.frame(Order = as.character(c(Group_1, Group_2)))
#my_colour = list(Protein = c("HXT1" = "#66c2a5", "HXT3" = "#fc8d62", "HXT5" = "#8da0cb", 
#"HXT7" = "#e78ac3", "Others" ="#CECED2"))

my_colour = list(Order = c("First" = "#1b9e77", "Second" = "#d95f02"))
row.names(row_ann) = row.names(PPI_HXT5_all_heatmap)

fit_heatmap = pheatmap(PPI_HXT5_all_heatmap, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames=FALSE,
                       annotation_colors = my_colour, annotation_row = row_ann, 
                       show_colnames=T, col = color_scale)


save_pheatmap_pdf(fit_heatmap, "~/Desktop/PPI_HXT5_all_heatmap_cluster_row.pdf")