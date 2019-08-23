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
PPI_fit = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_neg_zero.csv")
GO_slim = as.matrix(read.table("Working_data/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
Gene_carbon = unique(GO_slim[which(GO_slim[,6] == "GO:0008643"), 1])
Gene_transcription = unique(GO_slim[which(GO_slim[,6] == "GO:0006352"), 1])
Gene_translation = unique(GO_slim[which(GO_slim[,6] == "GO:0006413")])
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

PPI_carbon = check_specific_protein(PPI_fit, Gene_carbon) # 284
PPI_transcription = check_specific_protein(PPI_fit, Gene_transcription) # 731
PPI_translation = check_specific_protein(PPI_fit, Gene_translation) #658

# Make a heatmap for each groups
HXT = c("YHR094C","YDR345C","YHR096C", "YDR342C") # (HXT1, HXT3, HXT5, HXT7)
PPI_fit = dataFrameReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_neg_zero.csv")
PPI_carbon_fit = PPI_fit[which(as.character(PPI_fit[,1]) %in% PPI_carbon),] #284
PPI_carbon_order = sort(PPI_carbon_fit[,1]) #284
PPI_carbon_fit_order = PPI_carbon_fit[match(PPI_carbon_order, PPI_carbon_fit[,1]),]
PPI_HXT1 = PPI_carbon_fit[grep("YHR094C", PPI_carbon_fit[,1]),] #71
Group = rep("HXT1", nrow(PPI_HXT1))
PPI_HXT1 = data.frame(PPI_HXT1,Group )
PPI_HXT3 = PPI_carbon_fit[grep("YDR345C", PPI_carbon_fit[,1]),] #68
Group = rep("HXT3", nrow(PPI_HXT3))
PPI_HXT3 = data.frame(PPI_HXT3, Group)
PPI_HXT5 = PPI_carbon_fit[grep("YHR096C", PPI_carbon_fit[,1]),] #45
Group = rep("HXT5", nrow(PPI_HXT5))
PPI_HXT5 = data.frame(PPI_HXT5, Group)
PPI_HXT7 = PPI_carbon_fit[grep("YDR342C", PPI_carbon_fit[,1]),] #74
Group = rep("HXT7", nrow(PPI_HXT7))
PPI_HXT7 = data.frame(PPI_HXT7, Group)
PPI_HXT = rbind(PPI_HXT1, PPI_HXT3, PPI_HXT5, PPI_HXT7) # 258
PPI_HXT= PPI_HXT[!duplicated(PPI_HXT[,1]),] #255
PPI_other = PPI_carbon_fit_order[which(!as.character(PPI_carbon_fit_order[,1]) %in% as.character(PPI_HXT[,1])),] #29
Group = rep("Others", nrow(PPI_other))
PPI_other = data.frame(PPI_other, Group)
PPI_carbon_final = rbind(PPI_HXT, PPI_other) 
## add another annotations (New versus reported)
PPI_reported = dataFrameReader_T("Working_data/multiple_validated_PPI.csv")
PPI_carbon_final_reported = match_both_direction(PPI_carbon_final, PPI_reported[,1])
Reported = rep("No", nrow(PPI_carbon_final))
Reported[which(PPI_carbon_final[,1] %in% PPI_carbon_final_reported[,1])] = "Yes" # 43
PPI_carbon_final = data.frame(PPI_carbon_final, Reported)
## add another annotations (Environment number)
PPI_count = dataFrameReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
count = PPI_count[match(PPI_carbon_final[,1], PPI_count[,1]),2]
PPI_carbon_final = data.frame(PPI_carbon_final, count)

# Cluster PPIs with each HXT protein, respectively and then get the order for the final matrix
PPI_HXT1_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT1"),] # 71
PPI_HXT3_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT3"),] # 68
PPI_HXT5_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT5"),] # 43
PPI_HXT7_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "HXT7"),] # 73
PPI_others_matrix = PPI_carbon_final[which(PPI_carbon_final$Group == "Others"),] # 29

cluster_order = function(PPI_HXT1_matrix){
        PPI_fitness = PPI_HXT1_matrix[,4:12]
        hc = hclust(dist(PPI_fitness), method = "complete")
        PPI_HXT1_matrix_order = PPI_HXT1_matrix[hc$order,]
        return(PPI_HXT1_matrix_order)

}
PPI_HXT1_matrix_order = cluster_order(PPI_HXT1_matrix)
PPI_HXT3_matrix_order = cluster_order(PPI_HXT3_matrix)
PPI_HXT5_matrix_order = cluster_order(PPI_HXT5_matrix)
PPI_HXT7_matrix_order = cluster_order(PPI_HXT7_matrix)
PPI_others_matrix_order = cluster_order(PPI_others_matrix)

PPI_carbon_final_order = rbind(PPI_HXT1_matrix_order, PPI_HXT3_matrix_order, 
                               PPI_HXT5_matrix_order, PPI_HXT7_matrix_order,
                               PPI_others_matrix_order)
# Take the order the same with Figure2A
PPI_carbon_heatmap = PPI_carbon_final_order[,c(4,8,10,12,7,9,5,11,6)]
colnames(PPI_carbon_heatmap) = c("SD", "Forskolin","NaCl", "FK506", "Doxorubicin", 
                                 "Raffinose",  "H2O2",  "16 \u00B0C", "Hydroxyurea")
rownames(PPI_carbon_heatmap) = as.character(PPI_carbon_final[,1])
#csvWriter(PPI_carbon_heatmap, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_pair_GO/environment/carbonhydrate_transport_network/PPI_carbohydrate_transport_heatmap.csv")

'''
row_ann = data.frame(Protein = as.character(PPI_carbon_final$Group), 
                     Reported = as.character(PPI_carbon_final$Reported),
                     Environment = as.character(PPI_carbon_final$count))
my_colour = list(
  Protein = c("HXT1" = "#7b3294", "HXT3" = "#c2a5cf", "HXT5" = "#d01c8b", "HXT7" = "#a6dba0", "Others" ="#008837"),
  Reported = c("Yes" = apple_colors[10], "No"=  apple_colors[8]),
  Environment = c("1" = "#4575b4", "2" = "#74add1", "3" = "#abd9e9", "4" = "#e0f3f8", "5" = "#ffffbf",
                  "6" = "#fee090", "7" = "#fdae61", "8" = "#f46d43", "9" = "#d73027")
)
'''


row_ann = data.frame(Protein = as.character(PPI_carbon_final_order$Group))
#my_colour = list(Protein = c("HXT1" = "#66c2a5", "HXT3" = "#fc8d62", "HXT5" = "#8da0cb", 
                             #"HXT7" = "#e78ac3", "Others" ="#CECED2"))

my_colour = list(Protein = c("HXT1" = "#1b9e77", "HXT3" = "#d95f02", "HXT5" = "#7570b3", 
                             "HXT7" = "#e7298a", "Others" ="#CECED2"))
row.names(row_ann) = row.names(PPI_carbon_heatmap)
col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=100)

library(pheatmap)
fit_heatmap = pheatmap(PPI_carbon_heatmap, cluster_rows = FALSE, cluster_cols = T, show_rownames=FALSE,
                       annotation_colors = my_colour,
                       labels_col  = c("SD", "Forskolin","NaCl", "FK506", "Doxorubicin", 
                                     "Raffinose",  expression('H'[2]* 'O'[2]),  "16 \u00B0C", "Hydroxyurea"),
                       annotation_row = row_ann, show_colnames=T, col = color_scale)

save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(fit_heatmap, "Working_figure/Figure4/Figure4A_carbonhydrate_transport_fitness_environment_primary_cluster.pdf")


#################################################################################


##### Transcription correlated PPIs
setwd("~/Dropbox/PPiSeq_02/")
PPI_fit = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
GO_slim = as.matrix(read.table("Working_data/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
Gene_carbon = unique(GO_slim[which(GO_slim[,6] == "GO:0008643"), 1])
Gene_transcription = unique(GO_slim[which(GO_slim[,6] == "GO:0006352"), 1])
Gene_translation = unique(GO_slim[which(GO_slim[,6] == "GO:0006413")])
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

#PPI_carbon = check_specific_protein(PPI_fit, Gene_carbon) # 284
PPI_transcription = check_specific_protein(PPI_fit, Gene_transcription) # 731
PPI_translation = check_specific_protein(PPI_fit, Gene_translation) #658

PPI_fit = dataFrameReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
PPI_transcription_fit = PPI_fit[which(as.character(PPI_fit[,1]) %in% PPI_transcription),] #731
PPI_transcription_order = sort(PPI_transcription_fit[,1]) #731
PPI_transcription_fit_order = PPI_transcription_fit[match(PPI_transcription_order, PPI_transcription_fit[,1]),]
PPI_MED11 = PPI_transcription_fit[grep("YMR112C", PPI_transcription_fit[,1]),] #521
Group = rep("MED11", nrow(PPI_MED11))
PPI_MED11 = data.frame(PPI_MED11,Group )
PPI_other = PPI_transcription_fit_order[which(!PPI_transcription_fit_order[,1] %in% PPI_MED11[,1]),] # 210
Group = rep("Others", nrow(PPI_other))
PPI_other = data.frame(PPI_other, Group)

PPI_transcription_final = rbind(PPI_MED11, PPI_other)
## add another annotations (New versus reported)
PPI_reported = dataFrameReader_T("Working_data/multiple_validated_PPI.csv")
PPI_transcription_final_reported = match_both_direction(PPI_transcription_final, PPI_reported[,1])
Reported = rep("No", nrow(PPI_transcription_final))
Reported[which(PPI_transcription_final[,1] %in% PPI_transcription_final_reported[,1])] = "Yes" # 43
PPI_transcription_final = data.frame(PPI_transcription_final, Reported)
## add another annotations (Environment number)
PPI_count = dataFrameReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
count = PPI_count[match(PPI_transcription_final[,1], PPI_count[,1]),2]
PPI_transcription_final = data.frame(PPI_transcription_final, count)

# Take the order the same with Figure2A
PPI_transcription_heatmap = PPI_transcription_final[,c(4,8,12,10,9,6,5,7,11)]
colnames(PPI_transcription_heatmap) = c("SD", "Forskolin", "FK506", "NaCl", "Raffinose", "Hydroxyurea",  
                                 "H2O2", "Doxorubicin", "16 \u00B0C")
rownames(PPI_transcription_heatmap) = as.character(PPI_transcription_final[,1])
row_ann = data.frame(Protein = as.character(PPI_transcription_final$Group), 
                     Reported = as.character(PPI_transcription_final$Reported),
                     Environment = as.character(PPI_transcription_final$count))
my_colour = list(
  Protein = c("MED11" = "#7b3294",  "Others" ="#008837"),
  Reported = c("Yes" = apple_colors[10], "No"=  apple_colors[8]),
  Environment = c("1" = "#4575b4", "2" = "#74add1", "3" = "#abd9e9", "4" = "#e0f3f8", "5" = "#ffffbf",
                  "6" = "#fee090", "7" = "#fdae61", "8" = "#f46d43", "9" = "#d73027")
)
row.names(row_ann) = row.names(PPI_transcription_heatmap)
col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=100)
library(pheatmap)
fit_heatmap = pheatmap(PPI_transcription_heatmap, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames=FALSE,
                       annotation_colors = my_colour,
                       annotation_row = row_ann, show_colnames=T, col = color_scale)
save_pheatmap_pdf <- function(x, filename, width=5, height=5) {
        pdf(filename, width = width, height = height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
}
save_pheatmap_pdf(fit_heatmap, "Working_figure/Figure4/Figure4_transcription_fitness_environment_primary.pdf")

##### Translation correlated PPIs
PPI_fit = dataFrameReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
PPI_translation_fit = PPI_fit[which(as.character(PPI_fit[,1]) %in% PPI_translation),] #658
PPI_translation_order = sort(PPI_translation_fit[,1]) #658
PPI_translation_fit_order = PPI_translation_fit[match(PPI_translation_order, PPI_translation_fit[,1]),]
PPI_DED1 = PPI_translation_fit[grep("YOR204W", PPI_translation_fit[,1]),] #521
Group = rep("DED1", nrow(PPI_DED1))
PPI_DED1 = data.frame(PPI_DED1,Group )
PPI_other = PPI_translation_fit_order[which(!PPI_translation_fit_order[,1] %in% PPI_DED1[,1]),] # 210
Group = rep("Others", nrow(PPI_other))
PPI_other = data.frame(PPI_other, Group)

PPI_translation_final = rbind(PPI_DED1, PPI_other) 
## add another annotations (New versus reported)
PPI_reported = dataFrameReader_T("Working_data/multiple_validated_PPI.csv")
PPI_translation_final_reported = match_both_direction(PPI_translation_final, PPI_reported[,1])
Reported = rep("No", nrow(PPI_translation_final))
Reported[which(PPI_translation_final[,1] %in% PPI_translation_final_reported[,1])] = "Yes" # 48
PPI_translation_final = data.frame(PPI_translation_final, Reported)
## add another annotations (Environment number)
PPI_count = dataFrameReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
count = PPI_count[match(PPI_translation_final[,1], PPI_count[,1]),2]
PPI_translation_final = data.frame(PPI_translation_final, count)

# Take the order the same with Figure2A
PPI_translation_heatmap = PPI_translation_final[,c(4,8,12,10,9,6,5,7,11)]
colnames(PPI_translation_heatmap) = c("SD", "Forskolin", "FK506", "NaCl", "Raffinose", "Hydroxyurea",  
                                        "H2O2", "Doxorubicin", "16 \u00B0C")
rownames(PPI_translation_heatmap) = as.character(PPI_translation_final[,1])
row_ann = data.frame(Protein = as.character(PPI_translation_final$Group), 
                     Reported = as.character(PPI_translation_final$Reported),
                     Environment = as.character(PPI_translation_final$count))
my_colour = list(
  Protein = c("DED1" = "#7b3294",  "Others" ="#008837"),
  Reported = c("Yes" = apple_colors[10], "No"=  apple_colors[8]),
  Environment = c("1" = "#4575b4", "2" = "#74add1", "3" = "#abd9e9", "4" = "#e0f3f8", "5" = "#ffffbf",
                  "6" = "#fee090", "7" = "#fdae61", "8" = "#f46d43", "9" = "#d73027")
)
row.names(row_ann) = row.names(PPI_translation_heatmap)
col_chosen = c(apple_colors[5], "#e7d4e8",apple_colors[7])
color_scale = colorRampPalette(col_chosen)(n=100)
library(pheatmap)
fit_heatmap = pheatmap(PPI_translation_heatmap, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames=FALSE,
                       annotation_colors = my_colour,
                       annotation_row = row_ann, show_colnames=T, col = color_scale)

save_pheatmap_pdf(fit_heatmap, "Working_figure/Figure4/Figure4_translation_fitness_environment_primary.pdf")



