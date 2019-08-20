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
#### First there might be different lineages that have fitness measurement (Use NA  or a specific numuber: -0.5)
PPI_norm = csvReader_T("Working_data/Positive_PPI_environment/Normalized_fitness_PPI_all_primary.csv")
PPI_count = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
PPI_norm_matrix = matrix(0, nrow(PPI_norm), 11)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
        for(j in 3:12){
                mean_fit = mean(as.numeric(na.omit(PPI_norm[i,c(j, j + 11)])))
                if(is.na(mean_fit)){
                        PPI_norm_matrix[i, j-1] = NA
                }
                else{
                        PPI_norm_matrix[i,j-1] = mean_fit
                }
                
        }  
}

sum(PPI_norm_matrix[,1] == PPI_count[,1])

PPI_both_pos = PPI_norm_matrix[which(PPI_count[,3] == "1" & PPI_count[,4] == "1"), c(1,2,3)] # 3417 cor 0.9382537
PPI_both_pos = cbind(PPI_both_pos, rep("SD + SD2 +", nrow(PPI_both_pos))) # 3417
PPI_both_non = na.omit(PPI_both_pos) # 3417
PPI_both_NA = PPI_both_pos[which(!PPI_both_pos[,1] %in% PPI_both_non[,1]),]
cor(as.numeric(PPI_both_non[,2]),as.numeric(PPI_both_non[,3]), method = "spearman")

PPI_SD_pos = PPI_norm_matrix[which(PPI_count[,3] == "1" & PPI_count[,4] == "0"), c(1,2,3)] # 1336
PPI_SD_pos = cbind(PPI_SD_pos, rep("SD + SD2 -", nrow(PPI_SD_pos)))
PPI_SD_pos_non = na.omit(PPI_SD_pos) # 1257 : 1336 

PPI_SD2_pos = PPI_norm_matrix[which(PPI_count[,3] == "0" & PPI_count[,4] == "1"), c(1,2,3)] # 1534
PPI_SD2_pos = cbind(PPI_SD2_pos, rep("SD - SD2 +", nrow(PPI_SD2_pos)))
PPI_SD2_pos_non = na.omit(PPI_SD2_pos) # 1260 : 1543 

### Check the overlap between each group with other environments
PPI_both_non_count = PPI_count[match(PPI_both_non[,1], PPI_count[,1]), c(1,3,5:ncol(PPI_count))]
count_matrix = PPI_both_non_count[,2:ncol(PPI_both_non_count)]
sum_count = rowSums(apply(count_matrix,2, as.numeric))
PPI_both_non_count = cbind(PPI_both_non, sum_count, count_matrix)
colnames(PPI_both_non_count) = c("PPI", "SD_fit", "SD2_fit", "label","Environment_number",
                                 "SD(2X)",colnames(count_matrix)[2:9])
csvWriter(PPI_both_non_count, "Working_data/Positive_PPI_environment/SD_replicate/SD+_SD2+_fit_count_summary.csv")

### SD + SD2 -
PPI_SD_pos_non_count = PPI_count[match(PPI_SD_pos_non[,1], PPI_count[,1]), c(1,3, 5:ncol(PPI_count))]
count_matrix = PPI_SD_pos_non_count[,2:ncol(PPI_SD_pos_non_count)]
sum_count = rowSums(apply(count_matrix,2, as.numeric))
PPI_SD_pos_non_count = cbind(PPI_SD_pos_non, sum_count, count_matrix)
colnames(PPI_SD_pos_non_count) = c("PPI", "SD_fit", "SD2_fit", "label","Environment_number",
                                 "SD",colnames(count_matrix)[2:9])
csvWriter(PPI_SD_pos_non_count, "Working_data/Positive_PPI_environment/SD_replicate/SD+_SD2-_fit_count_summary.csv")

### SD - SD2 +
PPI_SD2_pos_non_count = PPI_count[match(PPI_SD2_pos_non[,1], PPI_count[,1]), c(1,4, 5:ncol(PPI_count))]
count_matrix = PPI_SD2_pos_non_count[,2:ncol(PPI_SD2_pos_non_count)]
sum_count = rowSums(apply(count_matrix,2, as.numeric))
PPI_SD2_pos_non_count = cbind(PPI_SD2_pos_non, sum_count, count_matrix)
colnames(PPI_SD2_pos_non_count) = c("PPI", "SD_fit", "SD2_fit", "label","Environment_number",
                                   "SD2",colnames(count_matrix)[2:9])
csvWriter(PPI_SD2_pos_non_count, "Working_data/Positive_PPI_environment/SD_replicate/SD-_SD2+_fit_count_summary.csv")

PPI_replicate_matrix = rbind(PPI_both_pos, PPI_SD_pos, PPI_SD2_pos)
colnames(PPI_replicate_matrix) = c("PPI", "SD", "SD2", "Label")
csvWriter(PPI_replicate_matrix, "Working_data/Positive_PPI_environment/SD_replicate/SD_fitness_replicates.csv")

##### Scatter plot to show the correlation between SD replicates
SD_rep = dataFrameReader_T("Working_data/Positive_PPI_environment/SD_replicate/SD_fitness_replicates.csv") # 6296
SD_rep = na.omit(SD_rep) # 5934
SD_rep_both_pos = SD_rep[which(SD_rep$Label == "SD + SD2 +"),] # 3417 # 3417
SD_rep_SD_pos = SD_rep[which(SD_rep$Label == "SD + SD2 -"),] # 1336 # 1257  # 79
SD_rep_SD2_pos = SD_rep[which(SD_rep$Label == "SD - SD2 +"),] # 1543 # 1260 # 283

length(which(PPI_count[,3] == "1")) # 4753
length(which(PPI_count[,4] == "1")) # 4960
cor(SD_rep$SD, SD_rep$SD2, method = "spearman") # 0.5259217
cor(SD_rep_both_pos$SD, SD_rep_both_pos$SD2, method = "spearman") # 0.9274757
library(scales)
library(ggplot2)
ggplot() +
        geom_hex(aes(x= SD, y= SD2), SD_rep, bins = 60)+
        scale_fill_gradient(low= apple_colors[10], high = apple_colors[7])+
        #scale_fill_gradientn(colors = apple_colors[c(10,3,7)])+
        # linear regression is heavily afftected by these small fitness values
        #geom_smooth(aes(x = fit01, y = fit02), PPI_fit_final_data, method='lm',se = FALSE, 
        #color = "magenta3", linetype = 2, cex = 0.4)+
        
        #add a line that contain equal fitness values
        #geom_smooth(aes(x = seq(0, 1.2, by = 0.2), y = seq(0, 1.2, by = 0.2)), linetype =2,
                    #method='lm', se= FALSE, col= apple_colors[11], cex = 0.3)+
        annotate("text", x = 0, y = 1.1, label = expression(paste("Spearman's ", italic(r), " = 0.53")),  parse = TRUE, col = apple_colors[11]) +
        
        scale_color_manual('', breaks = c("Positive PPI"),
                           values = apple_colors[8]) +
        
        scale_y_continuous(name = "Fitness of SD2",
                           limits=c(-0.4, 1.2),
                           breaks=seq(-0.4,1.2, by =0.4),
                           labels = seq(-0.4,1.2, by= 0.4)) +
        scale_x_continuous(name = "Fitness of SD", 
                           limits=c(-0.4, 1.2),
                           breaks=seq(-0.4,1.2, by =0.4),
                           labels = seq(-0.4,1.2, by= 0.4))+
        theme(legend.position =c(0.9,0.2), legend.key=element_blank(), legend.text=element_text(size=10)) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))

ggsave("~/Dropbox/PPiSeq_02/working_figure/Figure1/Figure1E_correlation_SD_replicates_hexagonlot_SD.pdf", height =5, width =5)

