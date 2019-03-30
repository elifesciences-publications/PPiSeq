##########################
#### This script can be used to generate a plot of lineage trajectories 
#### and use gradient colors (fitness values) to label each lineage (Figure 1B).
#### Meanwhile generate 1C, a barplot to show how we call positive PPIs.
##########################

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

#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

PPI_lineages = dataFrameReader_T("Paper_data/DMSO_PPI_barcodes_fitness_counts.csv")
PPI_lineages_select = PPI_lineages[, c(1, 3, 4, 6:10)] 
for (i in 4:8){
        PPI_lineages_select[,i] = frequency(as.numeric(PPI_lineages_select[,i]))
}
PPI_lineages_select[PPI_lineages_select == "0"] = 1e-9 # 5081689

# Extract Positive and Negative control
PPI_pos_DHFR = PPI_lineages_select[which(PPI_lineages_select[,1] == "positive_DHFR"),]
PPI_pos_DHFR[,1] = "DHFR(+)"
PPI_neg_DHFR = PPI_lineages_select[which(PPI_lineages_select[,1] == "negative_non_DHFR"),]
PPI_neg_DHFR[,1] = "DHFR(-)"

# Get PPIs containing IMD3 (YLR432W) as bait, EPR3 as prey, 
# YIL143C (promiscuous protein) as bit, and YPL139C as prey. 
PPI_ERP3 = PPI_lineages_select[grep("YDL018C_", PPI_lineages_select[,1]),] # 5258
PPI_IMD3 = PPI_lineages_select[grep("_YLR432W", PPI_lineages_select[,1]),] # 3545
PPI_promiscuous_YPL139C = PPI_lineages_select[grep("YPL139C_", PPI_lineages_select[,1]),] # 6410
PPI_promiscuous_YIL143C = PPI_lineages_select[grep("_YIL143C", PPI_lineages_select[,1]), ] # 3425

# Extract PRS and RRS strains
PRS = PPI_lineages_select[grep("Pos", PPI_lineages_select[,1]),]
RRS = PPI_lineages_select[grep("Neg", PPI_lineages_select[,1]),]
PRS[,1] = "PRS"
RRS[,1] = "RRS"

# Extract strains that contain ORF X Fragment only strains
MATa_DHFR12 = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/PPI_barcodes/MATa_genome_combine.csv")
MATalpha_DHFR3 = csvReader_T("~/Dropbox/PPiSeq_02/Working_data//PPI_barcodes/MATalpha_genome_combine.csv")
PPI_DHFR12 = PPI_lineages_select[which(PPI_lineages_select[,2] %in% MATa_DHFR12[,3]),]
PPI_DHFR12[,1] = "ORF X DHFR[1,2]"
PPI_DHFR3 = PPI_lineages_select[which(PPI_lineages_select[,2] %in% MATalpha_DHFR3[,3]),]
PPI_DHFR3[,1] = "DHFR[3] X ORF"

################### Figure1B
##### A function to transform lineage matrix to data frame which can be plotted in ggplot2
transform_lineages <- function(x, PPI_pos, barcode_pos, fitness_pos, G0_pos){
        time_points = rep(c(0,3,6,9,12), nrow(x))
        lineages = as.numeric(c(t(x[, G0_pos:(G0_pos + 4)])))
        PPI = rep("0", length(time_points))
        barcode = rep("0", length(time_points))
        Fitness = rep(0, length(time_points))
        for(i in 1:nrow(x)){
                PPI[(i*5 -4): (i*5)] = as.character(rep(x[i,PPI_pos], 5))
                barcode[(i*5 -4): (i*5)] = as.character(rep(x[i,barcode_pos], 5))
                Fitness[(i*5 -4): (i*5)] = as.numeric(rep(x[i,fitness_pos], 5))
        }
        matrix_lineages = data.frame(PPI, barcode, Fitness, time_points, lineages)
        return (matrix_lineages)
}
# Trasnform the lineages into a format for plot
PPI_pos_DHFR_transform = transform_lineages(PPI_pos_DHFR, 1,2,3,4)
PPI_neg_DHFR_transform = transform_lineages(PPI_neg_DHFR, 1,2,3,4)
PRS_transform = transform_lineages(PRS, 1,2,3,4)
RRS_transform = transform_lineages(RRS, 1,2,3,4)
PPI_DHFR12_transform = transform_lineages(PPI_DHFR12, 1,2,3,4)
PPI_DHFR3_transform = transform_lineages(PPI_DHFR3, 1,2,3,4)
PPI_ERP3_transform = transform_lineages(PPI_ERP3, 1,2,3,4)
PPI_IMD3_transform = transform_lineages(PPI_IMD3, 1,2,3,4)
PPI_promiscuous_YPL139C_transform = transform_lineages(PPI_promiscuous_YPL139C, 1,2,3,4)
PPI_promiscuous_YIL143C_transform = transform_lineages(PPI_promiscuous_YIL143C, 1,2,3,4)

library(ggplot2)
library(scales)

Lineage_plot = function(PPI_pos_DHFR_transform, output){
        PPI_pos_DHFR_transform = PPI_pos_DHFR_transform
        ggplot(data = PPI_pos_DHFR_transform, aes(x = time_points, y= lineages, group = barcode, color = Fitness)) +
                geom_line(alpha =0.7, size = 0.3) +
                scale_color_gradientn(colors = apple_colors[c(5,10,7)], limits = c(-0.5, 1.3)) +
                scale_y_continuous(name = "Frequency",
                                   limits = c(1e-9, 1e-1),
                                   trans = log10_trans(),
                                   breaks = trans_breaks("log10", function(x) 10^x),
                                   labels = trans_format("log10", math_format(10^.x)))+
                scale_x_continuous(name = "Generation",
                                   breaks = seq(0,12, by =3),
                                   labels = seq(0,12, by= 3)) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                theme(axis.text.x = element_text(size = 10, color = "black"), 
                      axis.text.y.left = element_text(size = 10, color = "black")) + 
                theme(text = element_text(size=12))
        
        ggsave(output, width = 5, height = 4)
        
}

Lineage_plot(PPI_pos_DHFR_transform, "Working_figure/DMSO_lineage_plot/Positive DHFR lineages.pdf")
Lineage_plot(PPI_neg_DHFR_transform, "Working_figure/DMSO_lineage_plot/Negative DHFR lineages.pdf")
Lineage_plot(PRS_transform, "Working_figure/DMSO_lineage_plot/PRS lineages.pdf")
Lineage_plot(RRS_transform, "Working_figure/DMSO_lineage_plot/RRS lineages.pdf")
Lineage_plot(PPI_DHFR12_transform, "Working_figure/DMSO_lineage_plot/ORF X DHFR[1,2] lineages.pdf")
Lineage_plot(PPI_DHFR3_transform, "Working_figure/DMSO_lineage_plot/DHFR[3] X ORF lineages.pdf")
Lineage_plot(PPI_ERP3_transform, "Working_figure/DMSO_lineage_plot/ERP3 X ORF lineages.pdf")
Lineage_plot(PPI_IMD3_transform, "Working_figure/DMSO_lineage_plot/ORF X IMD3 lineages.pdf")
Lineage_plot(PPI_promiscuous_YIL143C_transform, "Working_figure/DMSO_lineage_plot/ORF X YIL143C.pdf")
Lineage_plot(PPI_promiscuous_YPL139C_transform, "Working_figure/DMSO_lineage_plot/YPL139C X ORF.pdf")

#Randomly choose some lineages after removing these controls strains and make a plot
all_control = do.call("rbind", list(PPI_pos_DHFR, PPI_neg_DHFR, PRS, RRS, PPI_DHFR12, PPI_DHFR3))# Combine all the controls
PPI_lineages_remaining = PPI_lineages_select[which(!PPI_lineages_select[,2] %in% all_control[,2]),]
PPI_lineages_remaining = PPI_lineages_remaining[order(as.numeric(PPI_lineages_remaining[,3]), decreasing = T),] # order the lineage based on their fitness

# Different groups
PPI_lineages_high = PPI_lineages_remaining[which(as.numeric(PPI_lineages_remaining[,3]) > 0.4 & as.numeric(PPI_lineages_remaining[,4]) > 1e-7),] # 5568
PPI_lineages_medium = PPI_lineages_remaining[which(as.numeric(PPI_lineages_remaining[,3]) > 0.1 & as.numeric(PPI_lineages_remaining[,3]) < 0.4
                                                   & as.numeric(PPI_lineages_remaining[,4]) > 1e-8 ),] # 2417400
PPI_lineages_neutral = PPI_lineages_remaining[which(as.numeric(PPI_lineages_remaining[,3]) > -0.2 & as.numeric(PPI_lineages_remaining [,3]) < 0.2
                                                    & as.numeric(PPI_lineages_remaining[,4]) > 1e-8),] # 4460803
PPI_lineages_low = PPI_lineages_remaining[which(as.numeric(PPI_lineages_remaining[,3]) > -0.5 & as.numeric(PPI_lineages_remaining [,3]) < -0.1
                                                & as.numeric(PPI_lineages_remaining[,4]) > 1e-8),] # 97172
# 1000 lineages with fitness > 0.4; 2000 lineages with fitness > 0.1 & < 0.4; 2000 lineages > -0.2 & < 0.2; 1000 lineages < -0.1
PPI_lineages_final = do.call("rbind", list(PPI_lineages_high[sample(1:nrow(PPI_lineages_high), 1000),],
                                           PPI_lineages_medium[sample(1:nrow(PPI_lineages_medium), 2000),],
                                           PPI_lineages_neutral[sample(1:nrow(PPI_lineages_neutral), 2000),],
                                           PPI_lineages_low[sample(1:nrow(PPI_lineages_low), 1000),]))
PPI_lineages_final_transform = transform_lineages(PPI_lineages_final, 1,2,3,4)

Lineage_plot(PPI_lineages_final_transform, "Working_figure/DMSO_lineage_plot/Figure1B_6000_representative_lineages.pdf")

##################### Figure 1C
### A plot to show the fitness distribution of positive control, negative control, 
### PPIs that contain fragments, several positive PPIs with different fitness values or Q-values,
### and negative PPIs with low fitness or large Q-values

#Input data of positive PPIs in DMSO
ORF_fragments = rbind(PPI_DHFR12, PPI_DHFR3) # 17558
ORF_fragments[,1] = "ORF X DHFR fragment"
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv") # 1459163
DMSO_pos = DMSO_mean[which(DMSO_mean[,7] == 1),] # 5211
#Here I only choose these reported PPIs in this figure
PPI_reported = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/multiple_validated_PPI.csv") # summary of BIOGRID data
DMSO_pos = DMSO_pos[which(DMSO_pos[,1] %in% PPI_reported[,1]),] # 862
fitness_pos = as.numeric(DMSO_pos[,3])
DMSO_pos_high = DMSO_pos[which(fitness_pos > 0.6 & fitness_pos < 0.8),]
DMSO_pos_medium = DMSO_pos[which(fitness_pos > 0.5 & fitness_pos < 0.6),]
DMSO_pos_medium_low = DMSO_pos[which(fitness_pos > 0.4 & fitness_pos < 0.5),]
DMSO_pos_low = DMSO_pos[which(fitness_pos > 0.2 & fitness_pos < 0.4),]
DMSO_neg = DMSO_mean[which(DMSO_mean[,7] == 0),] # 1453952
fitness_neg = as.numeric(DMSO_neg[,3])
DMSO_neg_medium_low = DMSO_neg[which(fitness_neg > 0.4 & as.numeric(DMSO_neg[,5]) > 0.1),]
DMSO_neg_low = DMSO_neg[which(fitness_neg > 0 & fitness_neg < 0.4),]

# Random sample one protein, return the PPI name
random_sample_one = function(DMSO_pos_high){
        a = sample(1:nrow(DMSO_pos_high), 1)
        return (DMSO_pos_high[a,1])
}
# get all the replciates from PPI_lineages data
#pos_PPI_high = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_pos_high)),] # YGR106C_YOR270C
pos_PPI_high = PPI_lineages_select[which(PPI_lineages_select[,1] == "YGR106C_YOR270C"),]
pos_PPI_high[,1] = "VOA1_VPH1"
#pos_PPI_medium = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_pos_medium)),] # YDR508C_YBR106W
pos_PPI_medium = PPI_lineages_select[which(PPI_lineages_select[,1] == "YDR508C_YBR106W"),]
pos_PPI_medium[,1] = "GNP1_SND3"
#pos_PPI_medium_low = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_pos_medium_low)),] # YGL077C_YGL203C
pos_PPI_medium_low = PPI_lineages_select[which(PPI_lineages_select[,1] == "YGL077C_YGL203C"),] 
pos_PPI_medium_low[,1] = "HNM1_KEX1"
#pos_PPI_low = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_pos_low)),] # YIL035C_YGL019W
pos_PPI_low = PPI_lineages_select[which(PPI_lineages_select[,1] == "YIL035C_YGL019W"),] # YIL038C_YNL091W
pos_PPI_low[,1] = "CKA1_CKB1"
#neg_PPI_medium_low = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_neg_medium_low)),] # YDR086C_YPR028W
neg_PPI_medium_low = PPI_lineages_select[which(PPI_lineages_select[,1] == "YDR086C_YPR028W"),]
neg_PPI_medium_low[,1] = "SSS1_YOP1" # not reported in all other environments 
#neg_PPI_low = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_neg_low)),] # YGL085W_YIL030C
neg_PPI_low = PPI_lineages_select[which(PPI_lineages_select[,1] == "YGL085W_YIL030C"),]
neg_PPI_low[,1] = "LCL3_SSM4" # not reported in all other environments

all_data = rbind(ORF_fragments, PPI_pos_DHFR, pos_PPI_high, pos_PPI_medium, pos_PPI_medium_low, 
                 pos_PPI_low, neg_PPI_medium_low, neg_PPI_low, PPI_neg_DHFR)
color_label = rep("Negative PPI", nrow(all_data))
pos_PPI_name = c(PPI_pos_DHFR[1,1], as.character(pos_PPI_high[1,1]), as.character(pos_PPI_medium[1,1]),
                   as.character(pos_PPI_medium_low[1,1]), as.character(pos_PPI_low[1,1]))
#Control_PPI_name = c(unique(PPI_neg_DHFR[1,1], neg_PPI_medium_low[1,1], neg_PPI_low[1,1]))
color_label[which(all_data[,1] %in% pos_PPI_name)] = "Positive PPI"
color_label[which(all_data[,1] == "ORF X DHFR fragment")] = "Negative control"

PPI = all_data[,1]
fitness = all_data[,3]
color = color_label
bar_plot_data = data.frame(PPI, fitness, color)
bar_plot_data$PPI = factor(bar_plot_data$PPI, 
                           levels = c("ORF X DHFR fragment", "DHFR(+)", "VOA1_VPH1",
                                      "GNP1_SND3", "HNM1_KEX1", "CKA1_CKB1",
                                      "SSS1_YOP1", "LCL3_SSM4", "DHFR(-)"))
library(ggplot2)

ggplot(bar_plot_data, aes(x = PPI, y = fitness, group = PPI, col = color))+
        #geom_boxplot()+
        geom_point(alpha =0.3) +
        #geom_boxplot(col = apple_colors[8])+
        #geom_jitter(position=position_jitter(0.1), size = 0.6, alpha = 0.3) + 
        stat_summary(fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
        scale_color_manual(name = "", breaks = c('Positive PPI', "Negative PPI", "Negative control"),
                           values  = apple_colors[c(6, 5,7)])  +
        scale_y_continuous(name = "Fitness", 
                           limits=c(-1, 1.2),
                           breaks = seq(-1,1.2, by =0.2),
                           labels = seq(-1,1.2, by= 0.2))+
  theme(legend.key=element_blank(), legend.position = "top")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = 60, hjust =1),
        axis.title.x = element_blank(),axis.text.y.left = element_text(size = 10, color = "black")) + 
  theme(text = element_text(size=12))+ theme(plot.margin = unit(c(0.2,0.1,0.4,0.5), "cm"))
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure1C_Calling_PPIs.pdf", width= 5, height = 5)

########################## Figure 1D and Figure 1E 
### Figure 1D: scatter plot to show the correlation between each pair of barcodes for the same PPI
### Figure 1E: scatter plot to show the correlation between mean fitness values of the same PPI with two directions

# Figure 1D: only consider positive PPIs (removing control strains)
setwd("~/Dropbox/PPiSeq_02/")
PPI_lineages = dataFrameReader_T("Paper_data/DMSO_PPI_barcodes_fitness_counts.csv")
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv") # 1459163
# First remove control strains in the data. These strains have larger number of replciates make the analysis more difficult.
PPI_RRS = DMSO_mean[grep("Neg_PPI", DMSO_mean[,1]),1] #97
PPI_PRS = DMSO_mean[grep("Pos_PPI", DMSO_mean[,1]),1] #108
PPI_pos = DMSO_mean[grep("positive_DHFR", DMSO_mean[,1]),1] # 1
PPI_neg = DMSO_mean[grep("negative_non_DHFR", DMSO_mean[,1]),1] # 1
PPI_control = c(PPI_PRS, PPI_RRS, PPI_pos, PPI_neg)
DMSO_pos = DMSO_mean[which(DMSO_mean[,7] == "1"),] # 5211
DMSO_pos_select = DMSO_pos[which(!DMSO_pos[,1] %in% PPI_control),] #5178
PPI_lineages_select= PPI_lineages[which(PPI_lineages[,1] %in% DMSO_pos_select[,1]),] #18411

# put the fitness values of replicates onto the same row
PPI_unique= unique(DMSO_pos_select[,1])
PPI_indiv_matrix= matrix(0, length(PPI_unique), 7)
PPI_indiv_matrix[,1]= PPI_unique
PPI_indiv_matrix[,2] = as.numeric(DMSO_pos_select[match(PPI_indiv_matrix[,1], DMSO_pos_select[,1]),2])
PPI_indiv_matrix[,3] = as.numeric(DMSO_pos_select[match(PPI_indiv_matrix[,1], DMSO_pos_select[,1]),3])
index = 0
for (i in 1:length(PPI_unique)){
        if(as.numeric(PPI_indiv_matrix[i,2]) == 4){
                PPI_indiv_matrix[i, 4]= PPI_lineages_select[index + 4 -3, 4]
                PPI_indiv_matrix[i, 5]= PPI_lineages_select[index + 4 -2, 4]
                PPI_indiv_matrix[i, 6]= PPI_lineages_select[index + 4 -1, 4]
                PPI_indiv_matrix[i, 7]= PPI_lineages_select[index + 4, 4]
                index = index + 4
        }else if(as.numeric(PPI_indiv_matrix[i,2] == 3)){
                PPI_indiv_matrix[i, 4]= PPI_lineages_select[index + 3 -2, 4]
                PPI_indiv_matrix[i, 5]= PPI_lineages_select[index + 3 -1, 4]
                PPI_indiv_matrix[i, 6]= PPI_lineages_select[index + 3, 4]
                index = index + 3
        }else if(as.numeric(PPI_indiv_matrix[i,2] == 2)){
                PPI_indiv_matrix[i, 4]= PPI_lineages_select[index + 2 -1, 4]
                PPI_indiv_matrix[i, 5]= PPI_lineages_select[index + 2, 4]
                index = index + 2
        }
}
colnames(PPI_indiv_matrix)= c("PPI", "Barcodes", "Mean_fitness","fit01", "fit02", "fit03", "fit04")

# Transfer this matrix into a matrix containing 3 colums: PPI, rep_01, rep_02
PPI_fit_matrix_01 = PPI_indiv_matrix[,c(1, 4,5)]
PPI_fit_matrix_02 = PPI_indiv_matrix[,c(1, 4,6)]
PPI_fit_matrix_03 = PPI_indiv_matrix[,c(1, 4,7)]
PPI_fit_matrix_04 = PPI_indiv_matrix[,c(1, 5,6)]
PPI_fit_matrix_05 = PPI_indiv_matrix[,c(1, 5,7)]
PPI_fit_matrix_06 = PPI_indiv_matrix[,c(1, 6,7)]
PPI_fit_all = rbind(PPI_fit_matrix_01, PPI_fit_matrix_02, PPI_fit_matrix_03,
                    PPI_fit_matrix_04, PPI_fit_matrix_05, PPI_fit_matrix_06) # 31068
# Remove any pair with at least one value >= 0
PPI_fit_final = PPI_fit_all[which(as.numeric(PPI_fit_all[,2]) != 0 & as.numeric(PPI_fit_all[,3]) != 0),] # 24887
cor(as.numeric(PPI_fit_final[,2]), as.numeric(PPI_fit_final[,3])) # 0.7000878

'''
######### Traditional method to plot correlation
lm_PPI = lm(as.numeric(PPI_fit_final[,3]) ~ as.numeric(PPI_fit_final[,2]))
predict_PPI = predict(lm_PPI)
lm_matrix = cbind(as.numeric(PPI_fit_final[,2]), predict_PPI)
lm_matrix = lm_matrix[order(lm_matrix[,2]),]
### Scatter plot
pdf("~/Desktop/Figure1D_correlation_two_replicates.pdf", width=5, height =5 )
plot(as.numeric(PPI_fit_final[,2]), as.numeric(PPI_fit_final[,3]), xlim = c(-0.4, 1.6),
     ylim = c(-0.4, 1.6), col = rgb(1,0,0,0.3), xlab = "Fitness of replicate 1", 
     ylab = "Fitness of replicate 2", pch = 16)
lines(lm_matrix[,1], lm_matrix[,2], col= apple_colors[11])
text(1.3,-0.3, "r = 0.7", col= apple_colors[11])
dev.off()
'''
####### Use ggplot to make scatter plots and hexagon plot
PPI = as.character(PPI_fit_final[,1])
fit01 = as.numeric(PPI_fit_final[,2])
fit02 = as.numeric(PPI_fit_final[,3])
PPI_fit_final_data = data.frame(PPI, fit01, fit02) # Transform the matrix into data.frame

library(ggplot2)
# scatter plot
ggplot() +
        geom_point(aes(x= fit01, y= fit02),PPI_fit_final_data, col = apple_colors[7], alpha = 0.3) +
        #geom_smooth(aes(x = fit01, y = fit02), method='lm',se = FALSE, color = apple_colors[11], cex = 0.6) +
        geom_smooth(aes(x = seq(0, 1.2, by = 0.2), y = seq(0, 1.2, by = 0.2)), method='lm',
                    se = FALSE, color = apple_colors[11], cex = 0.3, linetype = 2) +
        annotate("text", x = 1, y = 1.2, label = "italic(r) == 0.7", parse = TRUE, size = 5) + 
        #theme(legend.position =c(0.7, 0.12), legend.key=element_blank(), legend.title = element_blank())+
        #theme(legend.key=element_blank(), legend.position = "top")+
        scale_y_continuous(name = "Fitness of replicate1",
                           limits=c(0, 1.2),
                           breaks=seq(0,1.2, by =0.2),
                           labels = seq(0,1.2, by= 0.2)) +
        scale_x_continuous(name = "Fitness of replicate2", 
                           limits=c(0, 1.2),
                           breaks=seq(0,1.2, by =0.2),
                           labels = seq(0,1.2, by= 0.2)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black",hjust =1),
              axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("~/Dropbox/PPiSeq_02/working_figure/Figure1D_correlation_two_replicates_ggplot.pdf", width = 5, height = 5)

### Hexagon plot I think Hexagon plot is better than scatter plot
ggplot() +
        geom_hex(aes(x= fit01, y= fit02, fill = log10(..count..)), PPI_fit_final_data, bins = 60)+
        scale_fill_gradient(low= "white", high = apple_colors[7])+
        # linear regression is heavily afftected by these small fitness values
        #geom_smooth(aes(x = fit01, y = fit02), PPI_fit_final_data, method='lm',se = FALSE, 
        #color = "magenta3", linetype = 2, cex = 0.4)+
        
        #add a line that contain equal fitness values
        geom_smooth(aes(x = seq(0, 1.2, by = 0.2), y = seq(0, 1.2, by = 0.2)), linetype =2,
                    method='lm', se= FALSE, col= apple_colors[11], cex = 0.3)+
        annotate("text", x = 0.2, y = 1.1, label = "italic(r) == 0.7",  parse = TRUE, col = apple_colors[11]) +
        
        scale_color_manual('', breaks = c("Positive PPI"),
                           values = apple_colors[8]) +
        
        scale_y_continuous(name = "Fitness of replicate 2",
                           limits=c(0, 1.2),
                           breaks=seq(0,1.2, by =0.2),
                           labels = seq(0,1.2, by= 0.2)) +
        scale_x_continuous(name = "Fitness of replicate 1", 
                           limits=c(0, 1.2),
                           breaks=seq(0,1.2, by =0.2),
                           labels = seq(0,1.2, by= 0.2))+
        theme(legend.position ="right", legend.key=element_blank(), legend.text=element_text(size=10)) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black", hjust =1),
              axis.text.y.left = element_text(size = 10, color = "black"))

ggsave("~/Dropbox/PPiSeq_02/working_figure/Figure1D_correlation_two_replicates_hexagonlot.pdf", height =4, width =5)

######### Figure 1E: scatter plot to show the correlation between mean fitness values of the same PPI with two directions
setwd("~/Dropbox/PPiSeq_02/")
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv") # 1459163
DMSO_pos = DMSO_mean[which(DMSO_mean[,7] == "1"),] # 5211
# first find these PPIs that have two orientations in our data (check function of mark_duplicates_fast)
DMSO_pos_two = mark_duplicates_fast(DMSO_pos[,1]) # 4819  The second column are the opposite orientation (0 if none)
pos_PPI_two_direction = DMSO_pos_two[which(DMSO_pos_two[,2] != 0),] # 392
fitness_two= matrix(0, nrow(pos_PPI_two_direction), 2)
for(i in 1:nrow(fitness_two)){
        fitness_two[i,1] = DMSO_pos[which(DMSO_pos[,1] == pos_PPI_two_direction[i,1]),3]
        fitness_two[i,2] = DMSO_pos[which(DMSO_pos[,1] == pos_PPI_two_direction[i,2]),3]
}

#Transform the data into a data.frame
fit_01 = as.numeric(fitness_two[,1])
fit_02 = as.numeric(fitness_two[,2])
PPI = pos_PPI_two_direction[,1]
fitness_two_data = data.frame(PPI, fit_01, fit_02)
cor(fitness_two_data[,2], fitness_two_data[,3]) #0.58

library(ggplot2)
ggplot() +
        geom_point(aes(x= fit_01, y= fit_02), col = apple_colors[7], alpha = 0.8) +
        geom_smooth(aes(x = seq(0.2, 1, by= 0.2), y = seq(0.2, 1, by = 0.2)), 
                    method='lm',se = FALSE, color = apple_colors[11], cex = 0.4, linetype = 2) +
        scale_y_continuous(name = "Fitness of ORF1-DHFR[3] X ORF2-DHFR[1,2]",
                           limits=c(0.2, 1),
                           breaks=seq(0.2,1, by =0.2),
                           labels = seq(0.2,1, by= 0.2)) +
        scale_x_continuous(name = "Fitness of ORF1-DHFR[1,2] X ORF2-DHFR[3]", 
                           limits=c(0.2, 1),
                           breaks=seq(0.2,1, by =0.2),
                           labels = seq(0.2,1, by= 0.2)) +
        
        annotate("text", x = 0.3, y = 0.95, label = "italic(r) == 0.58", parse = TRUE, size = 5) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black", hjust =1),
              axis.text.y.left = element_text(size = 10, color = "black"))

ggsave("~/Dropbox/PPiSeq_02/working_figure/Figure1E_correlation_fitness_two_directions.pdf", width = 5, height = 5)

################# Figure 1F The performance of PRS and RRS. That is comparable with Y2H an PCA.
setwd("~/Dropbox/PPiSeq_02/")
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv") # 1459163
DMSO_pos = DMSO_mean[which(DMSO_mean[,7] == "1"),] # 5211
DMSO_RRS = DMSO_mean[grep("Neg_PPI", DMSO_mean[,1]),1] #97
DMSO_PRS = DMSO_mean[grep("Pos_PPI", DMSO_mean[,1]),1] #108
DMSO_pos_RRS = DMSO_pos[grep("Neg_PPI", DMSO_pos[,1]),1] # 29
DMSO_pos_PRS = DMSO_pos[grep("Pos_PPI", DMSO_pos[,1]),1] # 3
DMSO_pos_RRS_unique = unique(split_string_vector_line(DMSO_pos_RRS)[,1]) # 2
DMSO_pos_PRS_unique = unique(split_string_vector_line(DMSO_pos_PRS)[,1]) # 20
DMSO_RRS_unique = unique(split_string_vector_line(DMSO_RRS)[,1]) # 67
DMSO_PRS_unique = unique(split_string_vector_line(DMSO_PRS)[,1]) # 71
PRS_percent = length(DMSO_pos_PRS_unique)/length(DMSO_PRS_unique)* 100
RRS_percent = length(DMSO_pos_RRS_unique)/length(DMSO_RRS_unique) * 100
# make a barplot to show the percentage 
PRS_percent #28.17
RRS_percent #2.99

pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure1F_PRS_RRS_barplot.pdf", width = 3, height =4)
par(mar= c(3,4.5,1,1))
name = c("PRS", "RRS")
barCenter = barplot(c(PRS_percent, RRS_percent), horiz=F, beside=F, ylim = c(0,40),
                    ylab= "Fraction positive by PPiSeq (%)", space = c(0.2),
                    col = c(apple_colors[7], apple_colors[5]), axisnames = F,
                    border = NA)
text(x = barCenter+0.2, y = -3, adj = 1, labels = name, xpd = TRUE)
text(x = barCenter, y = c(PRS_percent + 2, RRS_percent + 2), labels = c("20/71", "2/67"))
dev.off()

################# Figure 1G Comparison among PPiSeq, DHFR_PCA, and BioGRID. 
### Make a venn diagram to show the overlap the current data with PCA and BIOGRID and fitness distributions for different parts
# Generate BIOGRID databse that not contain PCA PPIs
setwd("~/Dropbox/PPiSeq_02/")
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv") # 1459163
yeast_PPI = read.delim("Working_data/BIOGRID-ORGANISM-3.5.165.tab2/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.165.tab2.txt", header = T)
pos_PPI= yeast_PPI[which(yeast_PPI$Experimental.System.Type == "physical"), c(6:7, 12:14)] # 164992
nrow(pos_PPI[which(pos_PPI[,5] == 'Tarassov K (2008)'),]) # 2616
BIOGRID_PPI = pos_PPI[which(pos_PPI[,5] != 'Tarassov K (2008)'),] #remove mDHFR-PCA, 162376
BIOGRID_PPI_unique = unique(paste(BIOGRID_PPI[,1], BIOGRID_PPI[,2], sep = "_")) # 116665

## Obtain promiscuous proteins from PPiSeq data
fragment_select = dataFrameReader_T("Working_data/Promiscuous_PPIs/DMSO_Promiscuous_proteins.csv")
fragment_protein = unique(as.vector(split_string_vector(fragment_select[,1])))

## Remove PPIs that contain promiscuous proteins from BioGRID
promiscuous_BIOGRID_select = rep(0, length(BIOGRID_PPI_unique))
for(i in 1: length(promiscuous_BIOGRID_select)){
        PPI = split_string(BIOGRID_PPI_unique[i])
        if(PPI[1] %in% fragment_protein | PPI[2] %in% fragment_protein){
                promiscuous_BIOGRID_select[i] = 1
        }
}
length(which(promiscuous_BIOGRID_select == 1)) #113
BIOGRID_PPI_unique_filter = BIOGRID_PPI_unique[which(promiscuous_BIOGRID_select != 1)] #116552

# Here we only consider comparison of one direction PPIs among three datasets
BIOGRID_PPI_unique_deduplicate = mark_duplicates_fast(BIOGRID_PPI_unique_filter) # 106837

## Check which PPIs in BIOGRID after removing promisicuous PPIs that were covered by PPiseq
# match_both_direction: output the PPIs in the first input matrix that were matched two orientations of the second input
BIOGRID_PPiseq = match_both_direction(BIOGRID_PPI_unique_deduplicate, as.character(DMSO_mean[,1]))# 11480
BIOGRID_PPI_uncover = BIOGRID_PPI_unique_deduplicate[which(!BIOGRID_PPI_unique_deduplicate[,1] %in% BIOGRID_PPiseq[,1]),] # 95357

# Check which PPI in PCA database after removing promiscuous PPIS that were covered by PPiseq 
PPI_paper= read.delim(file = "Working_data/PPI_set_PCA_science.txt", sep=" ")
PPI_PCA= paste(PPI_paper[,1], PPI_paper[,3], sep="_") # 2770
# Remove PPIs that contain promiscuous proteins
promiscuous_PCA_select = rep(0, length(PPI_PCA))
for(i in 1: length(promiscuous_PCA_select)){
        PPI = split_string(PPI_PCA[i])
        if(PPI[1] %in% fragment_protein | PPI[2] %in% fragment_protein){
                promiscuous_PCA_select[i] = 1
        }
}
length(which(promiscuous_PCA_select == 1)) #35
PCA_PPI_unique_filter = PPI_PCA[which(promiscuous_PCA_select != 1)] #2735
PCA_PPI_unique_deduplicate = mark_duplicates_fast(PCA_PPI_unique_filter) # 2735 no duplicated PPIs (two orientations)
PCA_PPI_matrix = cbind(PCA_PPI_unique_filter, 1:length(PCA_PPI_unique_filter)) # add line number to make it can be input of match_both_direction
# Get PCA data covered by PPiseq
PCA_PPiseq = match_both_direction(PCA_PPI_matrix, as.character(DMSO_mean[,1])) #2310
PCA_PPiseq_uncover = PCA_PPI_unique_filter[which(!PCA_PPI_unique_filter %in% PCA_PPiseq[,1])] # 425

#Get one direction of positive PPIs called by PPiSeq
pos_PPI = DMSO_mean[which(DMSO_mean[,7] == "1"),]
#Remove all the control strains
pos_RRS =pos_PPI[grep("Neg_PPI", pos_PPI[,1]),1] #29
pos_PRS =pos_PPI[grep("Pos_PPI", pos_PPI[,1]),1] #3
pos_DHFR =pos_PPI[grep("positive_DHFR", pos_PPI[,1]),1] # 1
pos_PPI_real = pos_PPI[which(!pos_PPI[,1] %in% c(pos_RRS, pos_PRS, pos_DHFR)),1] # 5178
pos_PPI_unique_filter = mark_duplicates_fast(pos_PPI_real) # 4786

BIOGRID_PPiseq_overlap = match_both_direction(BIOGRID_PPiseq, pos_PPI_unique_filter[,1])#656
BIOGRID_PPiseq_PCA_overlap = match_both_direction(BIOGRID_PPiseq_overlap, as.character(PCA_PPiseq[,1])) #456
PCA_BIOGRID_overlap = match_both_direction(PCA_PPiseq, BIOGRID_PPiseq[,1]) #652
PCA_PPiseq_overlap = match_both_direction(PCA_PPiseq, pos_PPI_unique_filter[,1]) #1198

area1 = nrow(pos_PPI_unique_filter) #4786
area2 = nrow(PCA_PPiseq) # 2310
area3 = nrow(BIOGRID_PPiseq) # 11480
n12 = nrow(PCA_PPiseq_overlap) # 1198
n23 = nrow(PCA_BIOGRID_overlap) # 652
n13 = nrow(BIOGRID_PPiseq_overlap) # 656
n123 = nrow(BIOGRID_PPiseq_PCA_overlap) #456
library(VennDiagram)
pdf("~//Dropbox/PPiSeq_02/Working_figure/Figure1G_draft_Venn diagram of overlapped PPIs among PPiseq, PCA, BIOGRID.pdf", width = 5, height =5)
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c("PPiSeq", "mDHFR-PCA", "BIOGRID data excluding mDHFR-PCA"), euler.d = TRUE,
                 scaled = TRUE, col = apple_colors[c(5,3,6)], fill = apple_colors[c(5,3,6)], cex = 1.5)
dev.off()

# Also include BIOGRID Uncovered, PCA Uncovered, and the intersection between them
a = nrow(BIOGRID_PPI_uncover) # 95357
b = length(PCA_PPiseq_uncover) # 425
c = length(intersect(BIOGRID_PPI_uncover[,1], PCA_PPiseq_uncover)) # 97
a-c # 95260
b-c # 328
nrow(DMSO_mean) # 1459163
# Manually draw the figure in keynote
