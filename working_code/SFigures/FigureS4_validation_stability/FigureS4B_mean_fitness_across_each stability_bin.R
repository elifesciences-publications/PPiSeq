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

### Check the normalized fitness values for these PPIs
setwd("~/Dropbox/PPiseq_02/")

################################
# Figure S4A mean fitness for each stability bin
vScore = dataFrameReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
count_summary = csvReader_T("Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")
vScore_fit = vScore[,4:ncol(vScore)]
### Take the mean fitness for all the positive PPIs
mean_fitness_pos = rep(0, nrow(vScore))
for(i in 1:length(mean_fitness_pos)){
        pos_env_index = which(count_summary[i,3:ncol(count_summary)] == "1")
        pos_env_index_fit = pos_env_index + 3
        mean_fitness_pos[i] = mean(as.numeric(vScore[i, pos_env_index_fit]))
}
min(mean_fitness_pos) # 0.1121392
fitness_count = data.frame(mean_fitness_pos, count_summary[,2])
colnames(fitness_count) = c("Mean_fitness", "Env_count")
ggplot() +
        
        geom_boxplot(aes(x = Env_count, y = Mean_fitness), fitness_count, outlier.shape=NA) +
        geom_dotplot(aes(x = Env_count, y = Mean_fitness), fitness_count, 
                     binaxis="y",stackdir="center",binwidth=0.002, alpha=0.2, col = apple_colors[8]) +
        xlab("Number of environments in which a PPI is identified") +
        ylab("Mean fitness of a PPI across different environments") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black")) 
#theme(plot.margin = unit(c(1,1,2,1), "cm"))
ggsave("Working_figure/Sfigures/paper/FigureS4_validation_each_stability_bin/FigureS4A_Mean_fitness_PPI_each_stability_bin.pdf", width=5, height =5)

################# Figure S4B can be directly Figure 3B
### Make barplot to show the percentage
### put the reported and unreported on to the same figure
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/")
rep_PPI_matrix = dataFrameReader_T("Reported_validation_matrix_SD_merge.csv")
unrep_PPI_matrix = dataFrameReader_T("Unreported_validation_matrix_SD_merge.csv")
ratio_rep = rep_PPI_matrix[4,-1]
ratio_unrep = unrep_PPI_matrix[4,-1]
ratio_all = as.numeric(c(ratio_rep[1], ratio_unrep[1], ratio_rep[2], ratio_unrep[2], 
                         ratio_rep[3], ratio_unrep[3], ratio_rep[4], ratio_unrep[4],
                         ratio_rep[5], ratio_unrep[5], ratio_rep[6], ratio_unrep[6],
                         ratio_rep[7], ratio_unrep[7], ratio_rep[8], ratio_unrep[8],
                         ratio_rep[9], ratio_unrep[9]))


rep_PPI_matrix[1,] #   0       0       5       7      14      15      23      43      16
rep_PPI_matrix[3,] #   1       1       6       9      19      18      27      44      16
unrep_PPI_matrix[1,]# 9      22      22      13      28      38      30      37      26
unrep_PPI_matrix[3,]#13      32      31      20      33      45      32      38      27
counts_label = c("0/1", "9/13", "0/1", "22/32", "5/6", "22/31",
                 "7/9", "13/20", "14/19", "28/33", "15/18", "38/45",
                 "23/27", "30/32", "43/44", "37/38", "16/16", "26/27")
library(RColorBrewer)
#col_chosen = brewer.pal(3,"Dark2")[1:2]
col_chosen = apple_colors[c(1,4)]
pdf("~/Dropbox/PPiSeq_02/Working_figure/Sfigures/paper/FigureS4_validation_each_stability_bin/FigureS4B_Validation_bar_plot_merge_reported_unreported.pdf", 
    width= 6, height=5)
barCenter = barplot(ratio_all*100, horiz=F, beside=F, ylim=c(0,100), ylab="Validation rate (%)",
                    space= c(0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15,
                             0.4, 0.15, 0.4, 0.15, 0.4, 0.15, 0.4, 0.15),
                    col= col_chosen , axisnames=F, border=NA, cex.axis=0.8)
legend(-0.5,120, legend=c("Previously reported", "Previously unreported"),fill=col_chosen, cex=0.8, bty="n",
       border=FALSE, xpd = TRUE)
text(x= barCenter, y = ratio_all*100 + 2, labels = counts_label, cex=0.5, xpd = TRUE)
env_num_loc = rep(0, 9)
for(i in 1:9){
        env_num_loc[i] = mean(barCenter[(2*i-1):(2*i)])
}
text(x = env_num_loc, y = -8, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -16, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()


#### Figure S4C False negative rate in different environmentsxs
