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
mean_fit = as.numeric(PPI_fit[match(PPI_env_count[,1], PPI_fit[,1]),2])

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

PPI_env_abun = data.frame(PPI_env_count[,1], as.character(PPI_env_count[,2]), log2(PPI_abun_min))
PPI_env_abun = na.omit(PPI_env_abun) # 13504
colnames(PPI_env_abun) = c("PPI", "Environment", "Abundance")
PPI_env_abun$Environment = factor(PPI_env_abun$Environment, levels = as.character(1:9))
min(PPI_env_abun$Abundance) # 1.5
max(PPI_env_abun$Abundance) # 17.3
## Get the number of PPIs for each environment
table(PPI_env_abun$Environment)
# 1:8085; 2:1336; 3:650; 4:531; 5:552; 6:657; 7:567; 8:613; 9:513
### All PPIs
count = c("8085", "1336", "650", "531", "552", "657", "567", "613", "513")


ggplot(PPI_env_abun, aes(x = Environment, y = Abundance, group = Environment))+
        #geom_boxplot(col = apple_colors[5])+
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), col= apple_colors[5])+
        #geom_point(col = apple_colors[3])
        
        stat_summary(aes(x = Environment, y = Abundance, group = Environment), PPI_env_abun,
                     fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
        scale_x_discrete(limits = c("1","2", "3","4", "5", "6", "7", "8", "9"))+
        scale_y_continuous(name = "Log2(mean number of minor protein of a PPI per cell)", 
                           limits=c(0, 18),
                           breaks = seq(0,18, by =2),
                           labels = seq(0,18, by= 2))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))+
        xlab("Number of environments in which a PPI is observed") + 
        annotate("text", x = 1:9,  y = rep(1, 9), label = count)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2C_minor_protein_abundance/SFigure2C_all_protein_abundance.pdf", 
       width= 5, height = 5)

### Tease out the reported PPIs
reported_PPI = csvReader_T("Working_data/multiple_validated_PPI.csv")
PPI_env_count_reported = match_both_direction(PPI_env_count, reported_PPI[,1]) # 1853
mean_fit = as.numeric(PPI_fit[match(PPI_env_count_reported[,1], PPI_fit[,1]),2]) # 1853

PPI_abun = PPI_abundance(PPI_env_count_reported[,1], protein_abundance)
PPI_abun_min = rep(0, nrow(PPI_abun)) # take the smaller abundance value for each PPI
for(i in 1:nrow(PPI_abun)){
        PPI_abun_min[i] = min(PPI_abun[i,])
}

PPI_env_abun = data.frame(PPI_env_count_reported[,1], as.character(PPI_env_count_reported[,2]), log2(PPI_abun_min))
PPI_env_abun = na.omit(PPI_env_abun) # 13504
colnames(PPI_env_abun) = c("PPI", "Environment", "Abundance")
PPI_env_abun$Environment = factor(PPI_env_abun$Environment, levels = as.character(1:9))
min(PPI_env_abun$Abundance) # 7.6
max(PPI_env_abun$Abundance) # 16.53
## Get the number of PPIs for each environment
table(PPI_env_abun$Environment)
# 1:306; 2:144; 3:112; 4:115; 5:148; 6:212; 7:254; 8:338; 9:222
### All PPIs
count = c("306", "144", "112", "115", "148", "212", "254", "338", "222")

library(ggplot2)

ggplot(PPI_env_abun, aes(x = Environment, y = Abundance, group = Environment))+
        #geom_boxplot(col = apple_colors[5])+
        geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), col= apple_colors[5])+
        #geom_point(col = apple_colors[3])
        
        stat_summary(aes(x = Environment, y = Abundance, group = Environment), PPI_env_abun,
                     fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
        scale_x_discrete(limits = c("1","2", "3","4", "5", "6", "7", "8", "9"))+
        scale_y_continuous(name = "Log2(mean number of minor protein of a PPI per cell)", 
                           limits=c(6, 18),
                           breaks = seq(6, 18, by =2),
                           labels = seq(6,18, by= 2))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))+
        xlab("Number of environments in which a PPI is observed") + 
        annotate("text", x = 1:9,  y = rep(7, 9), label = count)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/SFigure2/SFigure2C_minor_protein_abundance/SFigure2C_all_protein_abundance_reported.pdf", 
       width= 5, height = 5)