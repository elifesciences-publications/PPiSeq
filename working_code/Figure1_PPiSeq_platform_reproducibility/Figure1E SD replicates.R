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
                        PPI_norm_matrix[i, j-1] = -0.6
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
PPI_SD_pos_non = na.omit(PPI_SD_pos) # 1257 : 1336 for -0.5

PPI_SD2_pos = PPI_norm_matrix[which(PPI_count[,3] == "0" & PPI_count[,4] == "1"), c(1,2,3)] # 1534
PPI_SD2_pos = cbind(PPI_SD2_pos, rep("SD - SD2 +", nrow(PPI_SD2_pos)))

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

library(scales)
#pdf("~/Desktop/Positive_PPI_fitness_SD_replicate.pdf", width =5, height = 5)
pdf("working_figure/Figure1/Figure1E_Positive_PPI_fitness_SD_replicate.pdf", width =5, height = 5)
plot(SD_rep_SD_pos[,2], SD_rep_SD_pos[,3], xlim = c(-0.6, 1.5), ylim = c(-0.6, 1.5),
     xlab = "Fitness in SD (SD1)", ylab = "Fitness in SD replicate (SD2)", type = "p", pch = 16,
     col = alpha("#984ea3", 0.3), bty = "n")
points(SD_rep_SD2_pos[,2], SD_rep_SD2_pos[,3], pch = 16, col = alpha("#ff7f00", 0.3))
points(SD_rep_both_pos[,2], SD_rep_both_pos[,3], pch = 16, col = alpha(apple_colors[7], 0.3))
legend(x = -0.5, y = 1.7, legend = c("SD2 + SD1 - (1260)", "SD2 + SD1 + (3417, r = 0.93)", "SD2 - SD1 + (1257)"),
       pch = c(16,16,16), col = c("#ff7f00", apple_colors[7], "#984ea3"), bty = "n", xpd = T)
text(x = -0.4, y = 0, "Fitnesses \n in SD2(283)")
text(x = -0.1, y = -0.5, "Fitneses \nin SD1 (79)")
dev.off()
#### 







SD_rep$Label = factor(SD_rep$Label)
library(ggplot2)
ggplot() + 
        geom_point(aes(x = SD, y = SD2, col = Label), SD_rep, alpha = 0.3)+
        scale_color_manual('', values = c("#ff7f00", "#984ea3", apple_colors[7]))+
        
        scale_y_continuous(name = "Fitness in SD2",
                           limits=c(-0.5, 1.5),
                           breaks=seq(-0.5,1.5, by =0.3),
                           labels = seq(-0.5, 1.5, by= 0.3)) +
        scale_x_continuous(name = "Fitness in SD", 
                           limits=c(-0.5, 1.5),
                           breaks=seq(-0.5,1.5, by =0.3),
                           labels = seq(-0.5,1.5, by= 0.3))+
        theme(legend.position = "top", legend.key=element_blank(), legend.text=element_text(size=10)) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("~/Desktop/Positive_PPI_fitness_SD_replicate.pdf", width = 5, height = 5)
##### The fitness of each lineage and mean fitness for each protein pair
SD = csvReader_T("Paper_data/SD2_PPI_barcodes_fitness_counts.csv") # 5331589
SD2 = csvReader_T("Paper_data/SD_PPI_barcodes_fitness_counts.csv") # 5013513
fit_SD = as.numeric(SD[,4])
fit_SD2 = as.numeric(SD2[,4])
fit_SD_correct = (fit_SD - min(fit_SD))/(max(fit_SD) - min(fit_SD))
fit_SD2_correct = (fit_SD2 - min(fit_SD2))/(max(fit_SD2) - min(fit_SD2))
SD[,4] = fit_SD_correct
SD2[,4] = fit_SD2_correct

BC_overlap = intersect(SD[,3], SD2[,3]) # 4896623
SD_overlap = as.numeric(SD[match(BC_overlap, SD[,3]),4])
SD2_overlap = as.numeric(SD2[match(BC_overlap, SD2[,3]),4])

cor(SD_overlap, SD2_overlap)# 0.2715489

PPI_SD = csvReader_T("Paper_data/SD_mean_fitness_positive.csv") # 1445535
PPI_SD2 = csvReader_T("Paper_data/SD2_mean_fitness_positive.csv") # 1520112
PPI_overlap = intersect(PPI_SD[,1], PPI_SD2[,1]) # 1427398

PPI_SD_overlap = as.numeric(PPI_SD[match(PPI_overlap, PPI_SD[,1]),3])
PPI_SD2_overlap = as.numeric(PPI_SD2[match(PPI_overlap, PPI_SD[,1]),3])

cor(PPI_SD_overlap, PPI_SD2_overlap) # 0.0001269587

PPI_SD_pos = PPI_SD[which(PPI_SD[,7] == "1"),] # 5178
PPI_SD2_pos = PPI_SD2[which(PPI_SD2[,7] == "1"),] # 5348

PPI_pos_overlap = intersect(PPI_SD_pos[,1], PPI_SD2_pos[,1]) # 3745

PPI_SD_pos_overlap = as.numeric(PPI_SD_pos[match(PPI_pos_overlap, PPI_SD_pos[,1]),3])
PPI_SD2_pos_overlap = as.numeric(PPI_SD2_pos[match(PPI_pos_overlap, PPI_SD2_pos[,1]),3])

cor(PPI_SD_pos_overlap, PPI_SD2_pos_overlap) # 0.9393566

PPI_SD_pos_non_overlap = as.numeric(PPI_SD_pos[which(!PPI_SD_pos[,1] %in% PPI_pos_overlap),3]) # 1433
min(PPI_SD_pos_non_overlap) # 0.25
max(PPI_SD_pos_non_overlap) # 1.05418
median(PPI_SD_pos_non_overlap) # 0.34

PPI_SD_pos[which(as.numeric(PPI_SD_pos[,3]) == max(PPI_SD_pos_non_overlap)),]

PPI_SD2_pos_non_overlap = as.numeric(PPI_SD2_pos[which(!PPI_SD2_pos[,1] %in% PPI_pos_overlap),3]) # 1603
min(PPI_SD2_pos_non_overlap) # 0.1028208
max(PPI_SD2_pos_non_overlap) # 0.92554
median(PPI_SD2_pos_non_overlap) # 0.34

PPI_SD2_pos[which(as.numeric(PPI_SD2_pos[,3]) == max(PPI_SD2_pos_non_overlap)),]

hist(PPI_SD_pos_non_overlap, breaks =seq(0.2, 1.1, by= 0.05) )


nrow(na.omit(DMSO_2))
nrow(na.omit(DMSO))
nrow(na.omit(Forskolin))


nrow(DMSO_2) # 5322
nrow(DMSO) # 5036
nrow(Forskolin) # 4232
DMSO_2_dup = mark_duplicates_fast(DMSO_2[,1]) # 4965
DMSO_dup = mark_duplicates_fast(DMSO[,1]) # 4645
Forskolin_dup = mark_duplicates_fast(Forskolin[,1]) # 3877

DMSO_2_DMSO_overlap = match_both_direction(DMSO_2_dup, DMSO_dup[,1]) # 3414
DMSO_2_DMSO_Forskolin_overlap = match_both_direction(DMSO_2_DMSO_overlap, Forskolin_dup[,1]) # 2933
DMSO_2_Forskolin_overlap = match_both_direction(DMSO_2_dup, Forskolin_dup[,1]) # 3131
DMSO_Forskolin_overlap = match_both_direction(DMSO_dup, Forskolin_dup[,1]) # 3235


area1 = nrow(DMSO_dup) #4645
area2 = nrow(DMSO_2_dup) #4965
area3 = nrow(Forskolin_dup) # 3877
n12 = nrow(DMSO_2_DMSO_overlap) # 3414
n23 = nrow(DMSO_2_Forskolin_overlap) # 3131
n13 = nrow(DMSO_Forskolin_overlap) # 3235
n123 = nrow(DMSO_2_DMSO_Forskolin_overlap) # 2933
library(VennDiagram)
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure1/Figure1H_SD_replicates/DMSO_DMSO2_Forskolin_venn_diagram.pdf", width = 5, height =5)
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c("SD", "SD_2", "Forskolin"), euler.d = TRUE,
                 scaled = TRUE, col = apple_colors[c(5,3,6)], fill = apple_colors[c(5,3,6)], cex = 1.5)
dev.off()

### Check the correlation between overlaps DMSO_2_DMSO versus Forskolin_DMSO
DMSO_2_overlap_fit = as.numeric(DMSO_2[match(DMSO_2_DMSO_overlap[,1], DMSO_2[,1]), 3]) # 
DMSO_overlap_fit_1 = as.numeric(DMSO[match(DMSO_2_DMSO_overlap[,1], DMSO[,1]),3])
SD_overlap = na.omit(cbind(DMSO_2_overlap_fit, DMSO_overlap_fit_1)) # 3388
cor(as.numeric(SD_overlap[,1]), as.numeric(SD_overlap[,2]), method = "spearman") # 0.9352761

DMSO_overlap_fit_2 = as.numeric(DMSO[match(DMSO_Forskolin_overlap[,1], DMSO[,1]),3]) # 3235
Forskolin_overlap_fit = as.numeric(Forskolin[match(DMSO_Forskolin_overlap[,1], Forskolin[,1]),3])
SD_forskolin_overlap = na.omit(cbind(DMSO_overlap_fit_2, Forskolin_overlap_fit)) # 3190
cor(as.numeric(SD_forskolin_overlap[,1]), as.numeric(SD_forskolin_overlap[,2]), method = "spearman") # 0.8811259

library(scales)
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure1/Figure1H_SD_replicates/Correlation between SD and SD_2 and correlation between SD and Forskolin.pdf", width = 10, height =5)
par(mfrow=c(1,2))
plot(as.numeric(SD_overlap[,1]), as.numeric(SD_overlap[,2]), pch = 16, 
     col = alpha(apple_colors[7], 0.1), type= "p", xlab= "Mean fitness in SD",
     ylab = "Mean fitness in SD replicate")
text(x = 0.3, y = 1.0, "Spearman's r = 0.94")
plot(as.numeric(SD_forskolin_overlap[,1]), as.numeric(SD_forskolin_overlap[,2]), 
     pch = 16, col = alpha(apple_colors[7], 0.1), type = "p", xlab = "Mean fitness in SD",
     ylab = "Mean fitness in Forskolin environment")
text(x= 0.5, y = 1.2, "Spearman's r = 0.88")
dev.off()
