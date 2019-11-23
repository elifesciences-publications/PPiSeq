
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

##################### Figure 1C
### A plot to show the fitness distribution of positive control, negative control, 
### PPIs that contain fragments, several positive PPIs with different fitness values or Q-values,
### and negative PPIs with low fitness or large Q-values
setwd("~/Dropbox/PPiSeq_02/")
PPI_lineages = dataFrameReader_T("Paper_data/Lineage_barcode_fitness_files/SD_PPI_barcodes_fitness_counts.csv")
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
bait_ORF_fragment = csvReader_T("Paper_data/PPiseq_barcodes/Bait_ORF_DHFR3_control.csv")
prey_ORF_fragment = csvReader_T("Paper_data/PPiseq_barcodes/DHFR12_prey_ORF_control.csv")

PPI_DHFR12 = PPI_lineages_select[which(PPI_lineages_select[,2] %in% prey_ORF_fragment[,3]),]
PPI_DHFR12[,1] = "ORF X DHFR[1,2]"
PPI_DHFR3 = PPI_lineages_select[which(PPI_lineages_select[,2] %in% bait_ORF_fragment[,3]),]
PPI_DHFR3[,1] = "DHFR[3] X ORF"

#Input data of positive PPIs in DMSO
ORF_fragments = rbind(PPI_DHFR12, PPI_DHFR3) # 17558
ORF_fragments[,1] = "ORF X DHFR fragment"
DMSO_mean = csvReader_T("Paper_data/PPI_mean_fitness_calling_files/SD_mean_fitness_positive.csv") # 1448699
DMSO_pos = DMSO_mean[which(DMSO_mean[,8] == 1),] # 5286
#Here I only choose these reported PPIs in this figure
PPI_reported = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/Useful_datasets/multiple_validated_PPI.csv") # summary of BIOGRID data
DMSO_pos = DMSO_pos[which(DMSO_pos[,1] %in% PPI_reported[,1]),] # 880
fitness_pos = as.numeric(DMSO_pos[,3])
DMSO_pos_high = DMSO_pos[which(fitness_pos > 0.6 & fitness_pos < 0.8),]
DMSO_pos_medium = DMSO_pos[which(fitness_pos > 0.5 & fitness_pos < 0.6),]
DMSO_pos_medium_low = DMSO_pos[which(fitness_pos > 0.4 & fitness_pos < 0.5),]
DMSO_pos_low = DMSO_pos[which(fitness_pos > 0.2 & fitness_pos < 0.4),]
DMSO_neg = DMSO_mean[which(DMSO_mean[,8] == 0),] 
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
pos_PPI_high[,1] = "VOA1 x VPH1"
#pos_PPI_medium = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_pos_medium)),] # YDR508C_YBR106W
pos_PPI_medium = PPI_lineages_select[which(PPI_lineages_select[,1] == "YDR508C_YBR106W"),]
pos_PPI_medium[,1] = "GNP1 x SND3"
#pos_PPI_medium_low = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_pos_medium_low)),] # YGL077C_YGL203C
pos_PPI_medium_low = PPI_lineages_select[which(PPI_lineages_select[,1] == "YGL077C_YGL203C"),] 
pos_PPI_medium_low[,1] = "HNM1 x KEX1"

#pos_PPI_low = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_pos_low)),] # YIL035C_YGL019W
pos_PPI_low = PPI_lineages_select[which(PPI_lineages_select[,1] == "YJL061W_YIL115C"),] # YJL061W_YIL115C
pos_PPI_low[,1] = "NUP82 x NUP159"

#neg_PPI_medium_low = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_neg_medium_low)),] # YDR086C_YPR028W
neg_PPI_medium_low = PPI_lineages_select[which(PPI_lineages_select[,1] == "YCR005C_YMR120C"),]
neg_PPI_medium_low[,1] = "CIT2 x ADE17" # not reported in all other environments 

#neg_PPI_low = PPI_lineages_select[which(PPI_lineages_select[,1] == random_sample_one(DMSO_neg_low)),] # YGL085W_YIL030C
neg_PPI_low = PPI_lineages_select[which(PPI_lineages_select[,1] == "YGL085W_YIL030C"),]
neg_PPI_low[,1] = "LCL3 x SSM4" # not reported in all other environments

ORF_fragments[,1] = "ORF x Null"

all_data = rbind(ORF_fragments, PPI_pos_DHFR, pos_PPI_high, pos_PPI_medium, pos_PPI_medium_low, 
                 pos_PPI_low, neg_PPI_medium_low, neg_PPI_low, PPI_neg_DHFR)
color_label = rep("Negative PPI", nrow(all_data))
pos_PPI_name = c(PPI_pos_DHFR[1,1], as.character(pos_PPI_high[1,1]), as.character(pos_PPI_medium[1,1]),
                 as.character(pos_PPI_medium_low[1,1]), as.character(pos_PPI_low[1,1]))
#Control_PPI_name = c(unique(PPI_neg_DHFR[1,1], neg_PPI_medium_low[1,1], neg_PPI_low[1,1]))
color_label[which(all_data[,1] %in% pos_PPI_name)] = "Positive PPI"
#color_label[which(all_data[,1] == "ORF X DHFR fragment")] = "Negative control"

PPI = all_data[,1]
fitness = all_data[,3]
color = color_label
bar_plot_data = data.frame(PPI, fitness, color)
bar_plot_data$PPI = factor(bar_plot_data$PPI, 
                           levels = c("ORF x Null","DHFR(-)", "LCL3 x SSM4","CIT2 x ADE17", 
                                      "NUP82 x NUP159", "HNM1 x KEX1", "GNP1 x SND3", "VOA1 x VPH1", 
                                      "DHFR(+)"))
bar_plot_data_control = bar_plot_data[which(bar_plot_data[,1] == "ORF x Null"),]
bar_plot_data_neg = bar_plot_data[which(bar_plot_data[,1] %in% c("DHFR(-)", "DHFR(+)")),]
bar_plot_data_real = bar_plot_data[which(!bar_plot_data[,1] %in% unique(c(as.character(bar_plot_data_control[,1]),
                                                                          as.character(bar_plot_data_neg[,1])))),]
bar_plot_data_real$PPI = factor(bar_plot_data_real$PPI, 
                                levels = c("LCL3 x SSM4","CIT2 x ADE17", 
                                           "NUP82 x NUP159", "HNM1 x KEX1", "GNP1 x SND3", "VOA1 x VPH1"))

library(ggplot2)

ggplot()+
  geom_dotplot(aes(x = PPI, y = fitness, group = PPI, fill = color, col = color), bar_plot_data_real, 
               binaxis="y", stackdir="center",  binwidth = 0.04, alpha =0.5,show.legend = FALSE)+
  geom_jitter(aes(x = PPI, y = fitness, group = PPI, col = color), bar_plot_data_neg, width = 0.1, alpha = 0.5) +
  geom_violin(aes(x = PPI, y = fitness, group = PPI, col = color), bar_plot_data_control, 
              draw_quantiles = c(0.25, 0.5, 0.75), show.legend = FALSE)+
  #geom_boxplot(aes(x = PPI, y = fitness, group = PPI, col = color), bar_plot_data_control, 
  #show.legend = FALSE)+
  
  scale_x_discrete(limits = c("ORF x Null","DHFR(-)", "LCL3 x SSM4","CIT2 x ADE17", 
                              "NUP82 x NUP159", "HNM1 x KEX1", "GNP1 x SND3", "VOA1 x VPH1", 
                              "DHFR(+)")) +
  
  stat_summary(aes(x = PPI, y = fitness, group = PPI, col = color),bar_plot_data,
               fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
  scale_color_manual(name = "", breaks = c("Negative PPI", 'Positive PPI'),
                     values  = apple_colors[c(5,7)]) +
  scale_fill_manual(name = "", breaks = c( "Negative PPI", 'Positive PPI'),
                    values  = apple_colors[c(5,7)])+
  scale_y_continuous(name = "Fitness", 
                     limits=c(-1, 1.2),
                     breaks = seq(-1,1.2, by =0.2),
                     labels = seq(-1,1.2, by= 0.2))+
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 0.5)))+
  theme(legend.key = element_blank(), legend.position = c(0.8,0.2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black", angle = 60, hjust =1),
        axis.title.x = element_blank(),axis.text.y.left = element_text(size = 10, color = "black")) + 
  theme(text = element_text(size=10))
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure1/Figure1D_Calling_PPIs_violin_jitter_dot.pdf", width= 4.5, height = 4.5)

