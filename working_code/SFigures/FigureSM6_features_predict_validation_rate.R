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

all_Tecan = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/Tecan_DMSO_positive_combine.csv")

vScore = csvReader_T("~/Dropbox/PPiseq_02/Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
count_summary = csvReader_T("~/Dropbox/PPiseq_02/Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")
SD_merge = csvReader_T("~/Dropbox/PPiseq_02/Paper_data/PPI_mean_fitness_calling_files/SD_merge_mean_fitness_positive.csv")
SD_select = count_summary[which(count_summary[,3] == "1"),]
vScore_select = vScore[which(vScore[,1] %in% SD_select[,1]),]
Tecan_select = all_Tecan[which(all_Tecan[,1] %in% vScore_select[,1]),]
vScore_Tecan_select = vScore[match(Tecan_select[,1], vScore[,1]),]
p_value_select = SD_merge[match(Tecan_select[,1], SD_merge[,1]), 6]
vScore_Tecan_combine = cbind(vScore_Tecan_select, p_value_select,Tecan_select[,ncol(Tecan_select)])
vScore_Tecan_combine[,ncol(vScore_Tecan_combine)] = (as.numeric(vScore_Tecan_combine[,ncol(vScore_Tecan_combine)]) <= 0.05)
colnames(vScore_Tecan_combine) = c(colnames(vScore_Tecan_select), "p-value", "Tecan_q")
csvWriter(vScore_Tecan_combine, "~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/Method/SFigureSM6_Modeling_validation_rate/Tecan_PPiseq_training_data.csv")


training_data = dataFrameReader_T("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/Method/SFigureSM6_Modeling_validation_rate/Tecan_PPiseq_training_data.csv")
library(ggplot2)
col_chosen = apple_colors[c(5,7)]
ggplot(training_data, aes(x = Environment_number, fill = Tecan_q, col = Tecan_q))+
  geom_density(alpha = 0.3)+
  scale_color_manual(name = "Validated", values = col_chosen)+
  scale_fill_manual(name = "Validated", values = col_chosen)+
  scale_x_continuous(name = expression(italic(n)), 
                     limits=c(1,9),
                     breaks = seq(1,9, by =1),
                     labels = seq(1,9, by= 1)) +
  ylab("Density") +
  guides(fill=guide_legend(ncol=3), col = guide_legend(ncol= 3))+
  theme(legend.key = element_blank(),legend.position = c(0.3, 0.9),
        legend.text=element_text(size=10),legend.title=element_text(size=10),
        legend.key.size = unit(0.4, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"),
        axis.title.y=element_text(size=10)) + 
  theme(text = element_text(size=10))

ggsave("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/Method/SFigureSM6_Modeling_validation_rate/Environment_number.pdf", width =4, height =4)


ggplot(training_data, aes(x = SD, fill = Tecan_q, col = Tecan_q))+
  geom_density(alpha = 0.3)+
  scale_color_manual(name = "Validated", values = col_chosen)+
  scale_fill_manual(name = "Validated", values = col_chosen)+
  scale_x_continuous(name = expression(italic(f)), 
                     limits=c(0.2, 1),
                     breaks = seq(0.2,1, by =0.1),
                     labels = seq(0.2,1, by= 0.1)) +
  ylab("Density") +
  guides(fill=guide_legend(ncol=3), col = guide_legend(ncol= 3))+
  theme(legend.key = element_blank(),legend.position = c(0.8, 0.9),
        legend.text=element_text(size=10),legend.title=element_text(size=10),
        legend.key.size = unit(0.4, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"),
        axis.title.y=element_text(size=10)) + 
  theme(text = element_text(size=10))

ggsave("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/Method/SFigureSM6_Modeling_validation_rate/SD_fitness.pdf", width =4, height =4)

training_data = dataFrameReader_T("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/Method/SFigureSM6_Modeling_validation_rate/Tecan_PPiseq_training_data.csv")
training_data$p.value = log10(training_data$p.value)
ggplot(training_data, aes(x = p.value, fill = Tecan_q, col = Tecan_q))+
  geom_density(alpha = 0.3)+
  scale_color_manual(name = "Validated", values = col_chosen)+
  scale_fill_manual(name = "Validated", values = col_chosen)+
  scale_x_continuous(name = expression(italic(p)), 
                     limits=c(-12, 0),
                     breaks = seq(-12,0, by =2),
                     labels = seq(-12,0, by= 2)) +
  ylab("Density") +
  guides(fill=guide_legend(ncol=3), col = guide_legend(ncol= 3))+
  theme(legend.key = element_blank(),legend.position = c(0.8, 0.95),
        legend.text=element_text(size=10),legend.title=element_text(size=10),
        legend.key.size = unit(0.4, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black"),
        axis.title.y=element_text(size=10)) + 
  theme(text = element_text(size=10))

ggsave("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/Method/SFigureSM6_Modeling_validation_rate/p_value.pdf", width =4, height =4)


cor(training_data$SD, training_data$p.value)
cor(training_data$SD, training_data$Environment_number)
cor(training_data$p.value, training_data$Environment_number)
plot(training_data$SD, training_data$p.value)
plot(training_data$SD, training_data$Environment_number)
