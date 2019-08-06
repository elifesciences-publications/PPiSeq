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
PPI_fit_csv = dataFrameReader_T("Working_data/Positive_PPI_environment/All_PPI_environments_normalized_fit.csv")
PPI_pos = dataFrameReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
# remove the potential promiscuous proteins
promiscuous_protein = dataFrameReader_T("Working_data/Promiscuous_PPIs/Promiscuous_protein_summary.csv")
promiscuous_protein = promiscuous_protein[which(promiscuous_protein[,2] >= 2), 1]
promiscuous_homo_dimer = paste(promiscuous_protein, promiscuous_protein, sep = "_")

PPI_pos = PPI_pos[which(PPI_pos$DMSO == 1),] # 4645
# A function to extract self-interacting PPIs
extract_self_PPI = function(DMSO_PPI){
  self_PPI = extract_repeat_PPI(as.character(DMSO_PPI[,1])) # 448
  DMSO_self_PPI = DMSO_PPI[which(DMSO_PPI[,1] %in% self_PPI),]
  return(DMSO_self_PPI)
}
homo_dimer_all = extract_self_PPI(PPI_fit_csv) # 484
homo_dimer_all = homo_dimer_all[which(!homo_dimer_all[,1] %in% promiscuous_homo_dimer),] # 482
homo_dimer_pos_index = extract_self_PPI(PPI_pos) # 213; 159 for DMSO

homo_dimer_pos = homo_dimer_all[which(homo_dimer_all[,1] %in% homo_dimer_pos_index[,1]),] # 213; 159 for DMSO
homo_dimer_non = homo_dimer_all[which(!homo_dimer_all[,1] %in% homo_dimer_pos_index[,1]),] # 269; 323 for DMSO
Pos = rep("1", nrow(homo_dimer_pos))
homo_dimer_pos = cbind(homo_dimer_pos, Pos)
Pos = rep("0", nrow(homo_dimer_non))
homo_dimer_non = cbind(homo_dimer_non, Pos)
homo_dimer_final = rbind(homo_dimer_pos, homo_dimer_non)

homo_dimer_final[,1] = split_string_vector(as.character(homo_dimer_final[,1]))[,1]
homo_dimer_final_median = rep(0, nrow(homo_dimer_final))
for(i in 1:nrow(homo_dimer_final)){
  #homo_dimer_final_median[i] = median(as.numeric(homo_dimer_final[i, 2:(ncol(homo_dimer_final)-1)]), na.rm = T)
  homo_dimer_final_median[i] = as.numeric(homo_dimer_final[i, 2])
}

protein_abundance = csvReader_T("~/Dropbox/PPiseq_02/Working_data/protein_abundance/table_S4.csv")
homo_dimer_abundance = log2(as.numeric(protein_abundance[match(homo_dimer_final[,1], protein_abundance[,1]),5]))

homo_dimer_abundance_comp = cbind(homo_dimer_final[,1], homo_dimer_final_median, homo_dimer_abundance, 
                                  as.character(homo_dimer_final$Pos)) # 484
homo_dimer_abundance_comp = na.omit(homo_dimer_abundance_comp) # 478; 444 for DMSO
cor(as.numeric(homo_dimer_abundance_comp[,2]), as.numeric(homo_dimer_abundance_comp[,3]), method = "spearman") # 0.3330563; 0.3562464 for DMSO

homo_dimer_abundance_comp_pos = homo_dimer_abundance_comp[which(homo_dimer_abundance_comp[,4] == "1"),] # 213; 159 for DMSO
cor(as.numeric(homo_dimer_abundance_comp_pos[,2]), as.numeric(homo_dimer_abundance_comp_pos[,3]), method = "spearman") # 0.3389362 ; 0.456619 for DMSO pos

homo_dimer_abundance_comp_neg = homo_dimer_abundance_comp[which(homo_dimer_abundance_comp[,4] == "0"),] # 265; 283 for DMSO
cor(as.numeric(homo_dimer_abundance_comp_neg[,2]), as.numeric(homo_dimer_abundance_comp_neg[,3]), method = "spearman") # 0.1124867; 0.1704747 for DMSO
colnames(homo_dimer_abundance_comp) = c("Protein", "Fitness", "Number","Positive" )
csvWriter(homo_dimer_abundance_comp, "Working_data/homo_dimer/Homo_dimer_fitness_versus_protein_abundance.csv")

library(scales)
col_chosen = alpha(apple_colors[c(5,7)], 0.3)
pdf("Working_data/homo_dimer/Homo_dimer_fitness_versus_protein_abundance_DMSO.pdf", width=5, height = 5)
plot(as.numeric(homo_dimer_abundance_comp_pos[1,2]), homo_dimer_abundance_comp_pos[1,3], type = "p", xlim = c(-0.4, 1.2),
     xlab = "Fitness of homo-dimer in SD", ylab = "Number of molecules per cell (log2)", axes = FALSE,
     col = col_chosen[2], pch = 16)
axis(1, seq(-0.3, 1.2, by = 0.3), as.character(seq(-0.3, 1.2, by = 0.3)))
axis(2, seq(6, 18, by = 2), as.character(seq(6, 18, by = 2)))
for(i in 2: nrow(homo_dimer_abundance_comp_pos)){
  points(as.numeric(homo_dimer_abundance_comp_pos[i,2]), homo_dimer_abundance_comp_pos[i,3], col = col_chosen[2], pch = 16)
}
for(j in 1: nrow(homo_dimer_abundance_comp_neg)){
  points(as.numeric(homo_dimer_abundance_comp_neg[j,2]), homo_dimer_abundance_comp_neg[j,3], col = col_chosen[1], pch = 16)
}
legend("bottomright", c("Negative: 283 (r= 0.17)", "Positive: 159 (r= 0.46)"), pch = c(16, 16), col = alpha(apple_colors[c(5,7)], alpha = 0.8),
       bty = "n")
dev.off()

