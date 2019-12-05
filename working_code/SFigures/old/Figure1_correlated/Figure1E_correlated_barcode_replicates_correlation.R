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

### Figure 1D: scatter plot to show the correlation between each pair of barcodes for the same PPI
# Figure 1D: only consider positive PPIs (removing control strains)
setwd("~/Dropbox/PPiSeq_02/")
PPI_lineages = dataFrameReader_T("Paper_data/FK506_PPI_barcodes_fitness_counts.csv")
DMSO_mean = csvReader_T("Paper_data/FK506_mean_fitness_positive.csv") # 1459163
# First remove control strains in the data. These strains have larger number of replciates make the analysis more difficult.
PPI_RRS = DMSO_mean[grep("Neg_PPI", DMSO_mean[,1]),1] #97
PPI_PRS = DMSO_mean[grep("Pos_PPI", DMSO_mean[,1]),1] #108
PPI_pos = DMSO_mean[grep("positive_DHFR", DMSO_mean[,1]),1] # 1
PPI_neg = DMSO_mean[grep("negative_non_DHFR", DMSO_mean[,1]),1] # 1
PPI_control = c(PPI_PRS, PPI_RRS, PPI_pos, PPI_neg)
DMSO_pos = DMSO_mean[which(DMSO_mean[,7] == "1"),] # 5069
DMSO_pos_select = DMSO_pos[which(!DMSO_pos[,1] %in% PPI_control),] #5036
PPI_lineages_select= PPI_lineages[which(PPI_lineages[,1] %in% DMSO_pos_select[,1]),] #18016

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
PPI_fit_final = PPI_fit_all[which(as.numeric(PPI_fit_all[,2]) != 0 & as.numeric(PPI_fit_all[,3]) != 0),] # 23724
cor(as.numeric(PPI_fit_final[,2]), as.numeric(PPI_fit_final[,3]), method = "spearman") # 0.7925929

####### Use ggplot to make scatter plots and hexagon plot
PPI = as.character(PPI_fit_final[,1])
fit01 = as.numeric(PPI_fit_final[,2])
fit02 = as.numeric(PPI_fit_final[,3])
PPI_fit_final_data = data.frame(PPI, fit01, fit02) # Transform the matrix into data.frame

library(ggplot2)
### Hexagon plot I think Hexagon plot is better than scatter plot
ggplot() +
        geom_hex(aes(x= fit01, y= fit02, fill = log10(..count..)), PPI_fit_final_data, bins = 60)+
        scale_fill_gradient(low= "white", high = apple_colors[7])+
        # linear regression is heavily afftected by these small fitness values
        #geom_smooth(aes(x = fit01, y = fit02), PPI_fit_final_data, method='lm',se = FALSE, 
        #color = "magenta3", linetype = 2, cex = 0.4)+
        
        #add a line that contain equal fitness values
        geom_smooth(aes(x = seq(0.3, 2.1, by = 0.3), y = seq(0.3, 2.1, by = 0.3)), linetype =2,
                    method='lm', se= FALSE, col= apple_colors[11], cex = 0.3)+
        annotate("text", x = 0.6, y = 2.3, label = expression(paste("Spearman's ", italic(r), " = 0.79")),  parse = TRUE, col = apple_colors[11]) +
        
        scale_color_manual('', breaks = c("Positive PPI"),
                           values = apple_colors[8]) +
        
        scale_y_continuous(name = "Fitness of replicate 2",
                           limits=c(0.3, 2.4),
                           breaks=seq(0.3,2.4, by =0.3),
                           labels = seq(0.3,2.4, by= 0.3)) +
        scale_x_continuous(name = "Fitness of replicate 1", 
                           limits=c(0.3, 2.4),
                           breaks=seq(0.3,2.4, by =0.3),
                           labels = seq(0.3,2.4, by= 0.3))+
        theme(legend.position =c(0.9,0.2), legend.key=element_blank(), legend.text=element_text(size=10)) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black", hjust =1),
              axis.text.y.left = element_text(size = 10, color = "black"))

ggsave("~/Dropbox/PPiSeq_02/working_figure/SFigures/Figure1_related/FigureSX_Replicate_pos_all_environments_related_Figure1D/Figure1D_other_environments/Figure1D_correlation_two_replicates_hexagonlot_FK506.pdf", height =5, width =5)
