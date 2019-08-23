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
#### Make a heatmap to show different groups of PPIs
setwd("~/Dropbox/PPiSeq_02/")

############## Raffinose environment remove the first 10 hour data very noisy

TECAN_data_AUC = function(data_file_name, map_name,time_window,time_window_low, matrix_name, bad_wells){
        ## robot_data_01 as MTX data
        robot_data_01= read.csv(data_file_name, skip = 19)
        robot_data_01 = robot_data_01[1:(nrow(robot_data_01)-2),]
        #Matix of ODs
        od_robot_01 = as.matrix(robot_data_01[, 3:ncol(robot_data_01)])
        map_index = 1:96
        od_robot_01 = od_robot_01[,which(!map_index %in% bad_wells)] # remove the data of bad wells
        #Vector of Times (s)
        y = as.vector(robot_data_01[,1])
        z = strsplit(y, "\t")
        time_robot_01 = 1:length(z)
        for(i in 1:length(time_robot_01)){
                time_robot_01[i] = as.numeric(z[[i]][2])/60
        }
        ## input the map 
        diploid_05_02 = read.csv(map_name)
        PPI_map = diploid_05_02[,1]
        PPI_map = PPI_map[which(!map_index %in% bad_wells)]
        PPI_unique = unique(PPI_map)
        PPI_unique = PPI_unique[which(PPI_unique != "Negative")]
        
        library(pracma)
        AUC_PPI = function(time, od_data, column_PPI){
                auc_PPI = rep(0, length(column_PPI))
                for(i in 1: length(column_PPI)){
                        auc_PPI[i] = trapz(time, od_data[,column_PPI[i]])
                }
                return(auc_PPI)
        }
        
        ## Calculate AUC 
        column_negative_PPI = which(PPI_map == "Negative")
        time_select_01 = time_robot_01[which(time_robot_01/60 <= time_window & time_robot_01/60 >= time_window_low)] ## only consider data before MTX_time_window_threshold
        od_robot_select_01 = od_robot_01[1:length(time_select_01),] # only consider data before MTX_time_window_threshold
        Negative_AUC_01 = AUC_PPI(time_select_01, od_robot_select_01, column_negative_PPI)
        
        
        ### Create a matrix that contain normalized AUCs (subtract negative) 
        ### for each PPI, and calculate mean and SD of these SDs
        auc_PPI_matrix = matrix(0, length(PPI_unique), 10)
        colnames(auc_PPI_matrix) = c("PPI","Mean_Neg_AUC", "SD_Neg_AUC", "Mean_raw_AUC", 
                                     "AUC_norm_01", "AUC_norm_02", "AUC_norm_03", 
                                     "Mean_norm_AUC", "SD_norm_AUC", "Raw_AUC_pvalue")
        auc_PPI_matrix[,1] = as.character(PPI_unique)
        auc_PPI_matrix[,2] = mean(Negative_AUC_01)
        auc_PPI_matrix[,3] = sd(Negative_AUC_01)
        Neg_AUC = mean(Negative_AUC_01)
        Mean_AUC = rep(0, nrow(auc_PPI_matrix))
        for(i in 1:length(PPI_unique)){
                column_each_PPI = which(PPI_map == PPI_unique[i])
                auc_each_PPI_01 = AUC_PPI(time_select_01, od_robot_select_01, column_each_PPI)
                Mean_AUC[i] = mean(auc_each_PPI_01)
        }
        auc_PPI_matrix[,4] = Mean_AUC
        Max_AUC = max(Mean_AUC)
        for(i in 1:length(PPI_unique)){
                column_each_PPI = which(PPI_map == PPI_unique[i])
                auc_each_PPI_01 = AUC_PPI(time_select_01, od_robot_select_01, column_each_PPI)
                auc_PPI_matrix[i,5:7] = (auc_each_PPI_01-Neg_AUC)/(Max_AUC - Neg_AUC)
                auc_PPI_matrix[i,8] = mean((auc_each_PPI_01-Neg_AUC)/(Max_AUC - Neg_AUC))
                auc_PPI_matrix[i,9] = sd((auc_each_PPI_01-Neg_AUC)/(Max_AUC - Neg_AUC))
                auc_PPI_matrix[i,10] = t.test(auc_each_PPI_01, Negative_AUC_01, alternative = "greater")$p.value
        }
        auc_PPI_matrix[,10] = p.adjust(auc_PPI_matrix[,10], "BH")
        write.csv(auc_PPI_matrix, matrix_name, row.names = F)
        #return(length(which(as.numeric(auc_PPI_matrix[,6]) <= 0.05 & as.numeric(auc_PPI_matrix[,5]) >= 0)))
        return(length(which(as.numeric(auc_PPI_matrix[,10]) <= 0.05)))
}


### Transcription environments
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/")
data_file_name = "transcription/data/2019-08-20_16-21_Diploid_transcription_01_SD_MTX_96_T7.txt"
map_name = "transcription/Final_map/Diploid_transcription_1.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "transcription/data/Diploid_transcription_SD_MTX_T7.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, time_window_low, matrix_name, bad_wells) #5 from 30


### Translation environments
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/")
data_file_name = "translation/data/2019-08-20_16-23_Diploid_translation_01_SD_MTX_96_T2.txt"
map_name = "translation/Final_map/Diploid_translation_1.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "translation/data/Diploid_translation_01_SD_MTX_T2.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, time_window_low, matrix_name, bad_wells) #5 from 30




