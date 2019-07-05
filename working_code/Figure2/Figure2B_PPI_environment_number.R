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

##################################################
# Figure2B make a pie plot to show the number of detected environments for each PPI 
# The idea is attach positive environments to each PPI 
setwd("~/Dropbox/PPiSeq_02/")
DMSO_pos = csvReader_T("Working_data/Positive_PPI_environment/DMSO_Pos_PPI_real.csv")
H2O2_pos = csvReader_T("Working_data/Positive_PPI_environment/H2O2_Pos_PPI_real.csv")
HU_pos = csvReader_T("Working_data/Positive_PPI_environment/Hydroxyurea_Pos_PPI_real.csv")
Dox_pos = csvReader_T("Working_data/Positive_PPI_environment/Doxorubicin_Pos_PPI_real.csv")
Forskolin_pos = csvReader_T("Working_data/Positive_PPI_environment/Forskolin_Pos_PPI_real.csv")
Raffinose_pos = csvReader_T("Working_data/Positive_PPI_environment/Raffinose_Pos_PPI_real.csv")
NaCl_pos = csvReader_T("Working_data/Positive_PPI_environment/NaCl_Pos_PPI_real.csv")
cold_pos = csvReader_T("Working_data/Positive_PPI_environment/Cold_16C_Pos_PPI_real.csv")
FK506_pos = csvReader_T("Working_data/Positive_PPI_environment/FK506_Pos_PPI_real.csv")
PPI_list = list(DMSO_pos[,1], H2O2_pos[,1], HU_pos[,1], Dox_pos[,1], Forskolin_pos[,1],
                Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1]) # store PPIs of one environment as an element in a list
all_PPI = unique(c(DMSO_pos[,1], H2O2_pos[,1], HU_pos[,1], Dox_pos[,1], Forskolin_pos[,1],
                   Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1])) #13050
all_PPI_unique = mark_duplicates_fast(all_PPI) # 12333
all_PPI_matrix = matrix(0, nrow(all_PPI_unique), 10)
all_PPI_matrix[,1] = all_PPI_unique[,1]
colnames(all_PPI_matrix)= c("PPI","DMSO", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
# Check each PPI in each environment, if reported, put a value of 1 into the space.
for(i in 1:nrow(all_PPI_matrix)){
  if(all_PPI_unique[i,2] != "0"){
    for(j in 1:length(PPI_list)){
      if(all_PPI_unique[i,1] %in% PPI_list[[j]] | all_PPI_unique[i,2] %in% PPI_list[[j]]){
        all_PPI_matrix[i,j +1] = 1
      }
    }    
  }else {
    for(j in 1:length(PPI_list)){
      if(all_PPI_unique[i,1] %in% PPI_list[[j]]){
        all_PPI_matrix[i,j +1] = 1
      }
    }    
  }
  
}

environment_number = rep(0, nrow(all_PPI_matrix)) # adding the column of 2:9 in each row will give the number of positive environment
for(i in 1:length(environment_number)){
  environment_number[i] = sum(as.numeric(all_PPI_matrix[i,2:10]))
}
all_PPI_matrix_final = cbind(all_PPI_matrix[,1], environment_number, all_PPI_matrix[,2:10])
csvWriter(all_PPI_matrix_final, "Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")

# Or make a barplot to show how many of them have been reproted
setwd("~/Dropbox/PPiSeq_02/")
all_PPI_matrix_final = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
reported_PPI = csvReader_T("Working_data/multiple_validated_PPI.csv")
matrix_PPI_env_rep = matrix(0, 2, 9)
for(i in 1:9){
  all = all_PPI_matrix_final[which(all_PPI_matrix_final[,2] == i),]
  all_reported = match_both_direction(all,reported_PPI[,1])
  all_unreported = all[which(!all[,1] %in% all_reported[,1]),]
  matrix_PPI_env_rep[1,i] = nrow(all_reported)
  matrix_PPI_env_rep[2,i] = nrow(all_unreported)
}
matrix_PPI_env_rep[1,] # 276 135 110 115 147 210 254 338 222
all_PPI_count = matrix_PPI_env_rep[1,] + matrix_PPI_env_rep[2,]
all_PPI_count # 6941 1300  640  534  554  662  570  619  513
ratio = matrix_PPI_env_rep[1,]/all_PPI_count
ratio # 0.03976372 0.10384615 0.17187500 0.21535581 0.26534296 0.31722054 0.44561404 0.54604200 0.43274854
ratio_reported = c("4.0%", "10.4%", "17.2%", "21.5%", "26.5%", "31.7%", "44.6%", "54.6%", "43.3%")
matrix_PPI_env_rep_reverse = matrix(0, nrow(matrix_PPI_env_rep), ncol(matrix_PPI_env_rep))
matrix_PPI_env_rep_reverse[1,] = matrix_PPI_env_rep[2,]
matrix_PPI_env_rep_reverse[2,] = matrix_PPI_env_rep[1,]
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2B_Number_environments_PPI_reproted.pdf", height = 5, width = 5)
barCenter = barplot(matrix_PPI_env_rep_reverse, horiz=F, beside=F, ylim=c(0,8000), ylab="Number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
                    col= apple_colors[c(2,1)], axisnames=F, border=NA)
legend("topright", legend=c("Previously reported", "Previously unreported"), 
       fill=apple_colors[c(1,2)], bty="n", border=FALSE)
text(x= barCenter, y = all_PPI_count + 300, labels = ratio_reported , 
     cex=0.6, xpd = TRUE, col= apple_colors[1]) # add cumulative number
text(x= barCenter, y = -500, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -1200, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

## Or make a pie plot for environment_number
environment_count = data.frame(table(environment_number))
## Environment number per positive PPI

library(ggplot2)
colfunc <- colorRampPalette(apple_colors[c(5,3,7)])
ggplot(environment_count, aes(x = factor(1), y = Freq, fill= environment_number)) +
  geom_bar(width =1, size = 0.1, color = "white", stat = "identity") + 
  geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 1.8)+
  coord_polar("y") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_fill_manual(values=colfunc(9),
                    name = "Number of environments",
                    breaks = 1:9,
                    labels = 1:9)

ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2/Figure2B_PPI_environment_distribution.pdf", width =5 , height = 5)

