###########################
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
#FigureS2C
# Check the correaltion between fitness and protein abundance for different groups of PPIs
setwd("~/Dropbox/PPiSeq_02/")
protein_abundance = csvReader_T("Working_data/protein_abundance/table_S4.csv")
PPI_fit = csvReader_T("Working_data/Positive_PPI_environment/Normalzied_fitness_PPI_all.csv")
PPI_env_count = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
#PPI_env_count = PPI_env_count[which(PPI_env_count[,3] == "1"),] # only consider DMSO environment
### (1) Only take the DMSO fitness
#PPI_fit_count = PPI_fit[match(PPI_env_count[,1], PPI_fit[,1]),c(3,13)] # 4786
#mean_fit = rep(0, nrow(PPI_env_count))
#for(i in 1:length(mean_fit)){
#        mean_fit[i] = mean(na.omit(as.numeric(PPI_fit_count[i,])))    
#}
### (2) Take all the mean fitness
#mean_fit = as.numeric(PPI_fit[match(PPI_env_count[,1], PPI_fit[,1]), 2])

### (3) Consider all PPIs
mean_fit = as.numeric(PPI_fit[match(PPI_env_count[,1], PPI_fit[,1]), 2])
PPI_abundance = function(PPI, protein_abundance){
        PP_pair = split_string_vector(PPI)
        protein_1_abun = as.numeric(protein_abundance[match(PP_pair[,1], protein_abundance[,1]),4])
        protein_2_abun = as.numeric(protein_abundance[match(PP_pair[,2], protein_abundance[,1]),4])
        PPI_pair_abun = cbind(protein_1_abun, protein_2_abun)
        return(PPI_pair_abun)
}

PPI_abun = PPI_abundance(PPI_env_count[,1], protein_abundance)
PPI_abun_min = rep(0, nrow(PPI_abun)) # take the smaller abundance value for each PPI
for(i in 1:nrow(PPI_abun)){
        PPI_abun_min[i] = min(PPI_abun[i,])
}

PPI_env_abun = data.frame(PPI_env_count[,1], as.character(PPI_env_count[,2]), mean_fit, log2(PPI_abun_min)) #4736
PPI_env_abun = na.omit(PPI_env_abun) # 13479
colnames(PPI_env_abun) = c("PPI", "Environment", "Fitness", "Abundance")
cor(PPI_env_abun[,3], PPI_env_abun[,4], method = "spearman") # 0.2130625

# Remove the abundance level smaller than 2^8
PPI_env_abun = PPI_env_abun[which(PPI_env_abun[,4] >= 8),] # 13124

pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2D_correlation_fit_abundance/FigureS2D_correlation_fit_abundance_all.pdf", width=20, height =15)
par(mfrow = c(3, 4))
CR= cor(PPI_env_abun[,3], PPI_env_abun[,4], method = "spearman")
CR = round(CR, digits = 4)
plot(PPI_env_abun[,3], PPI_env_abun[,4], xlab = "Mean fitness across environments",
     ylab = "Log2(mean number of minor protein of a PPI per cell)", type = "p",
     col = rgb(0,0,1, alpha = 0.3), pch = 16, cex =1.2, cex.lab = 1.5,
     main = paste("All PPIs (r = ", as.character(CR), ")", sep= ""))

PPI_chosen = PPI_env_abun[which(PPI_env_abun[,2] == 1),]
CR1 = cor(PPI_chosen[,3], PPI_chosen[,4], method = "spearman")
CR1 = round(CR1, digits = 4)
plot(PPI_chosen[,3], PPI_chosen[,4], xlab = "Mean fitness across environments",
          ylab = "Log2(mean number of minor protein of a PPI per cell)", type = "p",
          main = paste("1 environment (r = ", as.character(CR1), ")", sep=""),
          col = rgb(0,0,1, alpha = 0.5), pch = 16, cex = 1.2, cex.lab = 1.5)
     
for(i in 2:9){
        PPI_chosen = PPI_env_abun[which(PPI_env_abun[,2] == i),]
        CR2 = cor(PPI_chosen[,3], PPI_chosen[,4], method = "spearman")
        CR2 = round(CR2, digits = 4)
        plot(PPI_chosen[,3], PPI_chosen[,4], xlab = "Mean fitness across environments",
             ylab = "Log2(mean number of minor protein of a PPI per cell)", type = "p",
             main = paste(as.character(i), " environments (r = ", as.character(CR2), ")",sep=""),
             col = rgb(0,0,1, alpha = 0.5), pch = 16, cex = 1.2, cex.lab = 1.5)
}
dev.off()






