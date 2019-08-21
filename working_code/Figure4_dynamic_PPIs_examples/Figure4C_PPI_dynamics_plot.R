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

##### Input the Tecan data and calculate the AUC for each PPI

TECAN_data_AUC = function(data_file_name, map_name,time_window, matrix_name, bad_wells){
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
        time_select_01 = time_robot_01[which(time_robot_01/60 <= time_window)] ## only consider data before MTX_time_window_threshold
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

### SD environments
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-12_11-53_Diploid_carbo_01_DMSO_MTX_T7_96_T7.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_1.csv" 
time_window = 40
bad_wells = 0
matrix_name = "Diploid_carbo_01_SD_MTX_T7.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, matrix_name, bad_wells) #12 from 30


setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-12_11-53_Diploid_carbo_06_DMSO_MTX_96_T14.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_6.csv" 
time_window = 40
bad_wells = 0
matrix_name = "Diploid_carbo_06_SD_MTX_T14.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, matrix_name, bad_wells) #12 from 30

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-12_11-55_Diploid_carbo_07_DMSO_MTX_96_T2.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_7.csv" 
time_window = 40
bad_wells = 0
matrix_name = "Diploid_carbo_07_SD_MTX_T2.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, matrix_name, bad_wells) #4 from 30

### NaCl environments
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-16_14-01_Diploid_carbo_NaCl_MTX_01_96_T7.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_1.csv" 
time_window = 40
bad_wells = 0
matrix_name = "Diploid_carbo_01_NaCl_MTX_T7.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, matrix_name, bad_wells) #10 from 30


setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-16_15-02_Diploid_carbo_06_NaCl_MTX_96_T14.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_6.csv" 
time_window = 40
bad_wells = 0
matrix_name = "Diploid_carbo_06_NaCl_MTX_T14.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, matrix_name, bad_wells) #10 from 30

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-16_15-03_Diploid_carbo_07_NaCl_MTX_96_T2.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_7.csv" 
time_window = 40
bad_wells = 0
matrix_name = "Diploid_carbo_07_NaCl_MTX_T2.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, matrix_name, bad_wells) #3 from 30


##### Combine two environments
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
dip_01_SD = dataFrameReader_T("Diploid_carbo_01_SD_MTX_T7.csv")
dip_06_SD = dataFrameReader_T("Diploid_carbo_06_SD_MTX_T14.csv")
dip_07_SD = dataFrameReader_T("Diploid_carbo_07_SD_MTX_T2.csv")

dip_01_NaCl = dataFrameReader_T("Diploid_carbo_01_NaCl_MTX_T7.csv")
dip_06_NaCl = dataFrameReader_T("Diploid_carbo_06_NaCl_MTX_T14.csv")
dip_07_NaCl = dataFrameReader_T("Diploid_carbo_07_NaCl_MTX_T2.csv")

dip_SD = rbind(dip_01_SD, dip_06_SD, dip_07_SD)
dip_NaCl = rbind(dip_01_NaCl, dip_06_NaCl, dip_07_NaCl)

dip_SD_NaCl = matrix(0, nrow(dip_SD), 7)
colnames(dip_SD_NaCl) = c("PPI", "Dynamics_01", "Dynamics_02", "Dynamics_03",
                          "Mean_dynamics", "SD_dynamics", "p-value")
dip_SD_NaCl[,1] = as.character(dip_SD[,1])
for(i in 1:nrow(dip_SD_NaCl)){
        dynamic = as.numeric(dip_NaCl[i,5:7]) - as.numeric(dip_SD[i,5:7])
        dip_SD_NaCl[i, 2:4] = dynamic
        dip_SD_NaCl[i,5] = mean(dynamic)
        dip_SD_NaCl[i,6] = sd(dynamic)
        dip_SD_NaCl[i,7] = t.test(as.numeric(dip_NaCl[i,5:7]), as.numeric(dip_SD[i,5:7]),
                                  alternative = "two.sided")$p.value
}

csvWriter(dip_SD_NaCl, "NaCl_DMSO_dynamics_Tecan.csv")

##### Get fitness values from PPiSeq
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
PPI_fit = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
Tecan_SD_NaCl = csvReader_T("NaCl_DMSO_dynamics_Tecan.csv")
PPI_fit_tecan = PPI_fit[which(PPI_fit[,1] %in% Tecan_SD_NaCl[,1]),]
PPI_NaCl_SD = as.numeric(PPI_fit_tecan[,10]) - as.numeric(PPI_fit_tecan[,4])
PPI_NaCl_SD_fit = cbind(PPI_fit_tecan[,1], PPI_NaCl_SD)
PPI_NaCl_SD_Tecan = Tecan_SD_NaCl[match(PPI_NaCl_SD_fit[,1], Tecan_SD_NaCl[,1]),5]
PPI_NaCl_SD_fit_Tecan = cbind(PPI_NaCl_SD_fit, PPI_NaCl_SD_Tecan)

PPI_NaCl_SD_fit_Tecan = PPI_NaCl_SD_fit_Tecan[which(PPI_NaCl_SD_fit_Tecan[,1] != "YBL029C-A_YHR096C"),]
csvWriter(PPI_NaCl_SD_fit_Tecan, "NaCl_DMSO_Tecan_PPiSeq_comparison.csv")

Fit_comp = csvReader_T("NaCl_DMSO_Tecan_PPiSeq_comparison.csv")
cor(as.numeric(Fit_comp[,2]), as.numeric(Fit_comp[,3])) # 0.3382956

PPI_count= csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
PPI_SD_pos_NaCl_neg = PPI_count[which(PPI_count[,3] == "1" & PPI_count[,9] == "0"),]
PPI_SD_pos_NaCl_pos = PPI_count[which(PPI_count[,3] == "1" & PPI_count[,9] == "1"),]
PPI_SD_neg_NaCl_pos = PPI_count[which(PPI_count[,3] == "0" & PPI_count[,9] == "1"),]
PPI_SD_neg_NaCl_neg = PPI_count[which(PPI_count[,3] == "0" & PPI_count[,9] == "0"),]

Fit_comp_01 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_neg),] # 35
Fit_comp_02 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_pos),] # 17
Fit_comp_03 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_pos),] # 0
Fit_comp_04 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_neg),] # 29

Fitness_comparison = Fitness_comparison[which(Fitness_comparison[,1] %in% PPI_chosen[,1]),]
cor(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3])) # 0.2650078
cor(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3])) # 0.1057152
cor(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3])) # 0.5229693
cor(as.numeric(Fit_comp[,2]), as.numeric(Fit_comp[,3])) # 0.3382956

pdf("NaCl_DMSO_Tecan_PPiSeq.pdf", width = 5, height = 5)
plot(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3]), 
     xlab = " NaCl - SD by PPiSeq",
     ylab = " NaCl - SD by Tecan", 
     pch = 16, col = apple_colors[1],  xlim = c(-0.8, 1.0), ylim = c(-0.8, 1.2))
points(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3]), pch = 16, col = apple_colors[2])
#points(as.numeric(Fit_comp_03[,2]), as.numeric(Fit_comp_03[,3]), pch = 16, col = apple_colors[3])
points(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3]), pch = 16, col = apple_colors[4])
legend("topleft", legend= c("SD + NaCl -", "SD + NaCl +", "SD - NaCl -"), col = apple_colors[c(1,2,4)],
       pch = c(16, 16, 16), bty = "n")
dev.off()

### Only consider PPIs containing HXT5
Fit_comp = Fit_comp[grep("YHR096C", Fit_comp[,1]),] # 17
Fit_comp_01 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_neg),] # 35
Fit_comp_02 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_pos),] # 17
Fit_comp_03 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_pos),] # 0
Fit_comp_04 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_neg),] # 29

cor(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3])) # 0.2650078
cor(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3])) # 0.1057152
cor(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3])) # 0.5229693
cor(as.numeric(Fit_comp[,2]), as.numeric(Fit_comp[,3])) # 0.3382956

pdf("NaCl_DMSO_Tecan_PPiSeq_HXT5.pdf", width = 5, height = 5)
plot(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3]), 
     xlab = " NaCl - SD by PPiSeq",
     ylab = " NaCl - SD by Tecan", 
     pch = 16, col = apple_colors[1],  xlim = c(-0.8, 1.0), ylim = c(-0.8, 1.2))
points(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3]), pch = 16, col = apple_colors[2])
#points(as.numeric(Fit_comp_03[,2]), as.numeric(Fit_comp_03[,3]), pch = 16, col = apple_colors[3])
points(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3]), pch = 16, col = apple_colors[4])
legend("topleft", legend= c("SD + NaCl -", "SD + NaCl +", "SD - NaCl -"), col = apple_colors[c(1,2,4)],
       pch = c(16, 16, 16), bty = "n")
dev.off()

## Make a barplot to show HXT5 increasing

#Fit_comp_chosen = rbind(Fit_comp_01, Fit_comp_02)
PPI_HXT5 = t(Fit_comp)
PPI_HXT5= rbind(as.numeric(PPI_HXT5[2,]), as.numeric(PPI_HXT5[3,]))
PPI_HXT5 = PPI_HXT5[,order(PPI_HXT5[1,], decreasing = T)]
colnames(PPI_HXT5) = Fit_comp[,1]
pdf("barplot_HXT5_NaCl_SD.pdf", height =5, width = 5)
barplot(PPI_HXT5,horiz=F, beside=T, axisnames = F, ylab = "NaCl - SD", 
        ylim = c(-0.5, 0.8),
        main = "HXT5 containing PPIs", col = apple_colors[c(5,7)])
legend("topright", legend=c("PPiSeq", "TECAN"),
       fill=apple_colors[c(5,7)], cex=0.8, bty="n", border=FALSE)
dev.off()



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

setwd("~/Dropbox/PPiSeq_02/")
### Raffinose environments
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190820/")
data_file_name = "2019-08-18_15-44_Diploid_carbo_Raffinose_01_MTX_96_T7.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_1.csv" 
time_window = 48
time_window_low = 10
bad_wells = 0
matrix_name = "Diploid_carbo_01_Raffinose_MTX_T7.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, time_window_low, matrix_name, bad_wells) #1 from 30


setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190820/")
data_file_name = "2019-08-18_16-00_Diploid_carbo_06_Raffinose_MTX_96_T14.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_6.csv" 
time_window = 48
time_window_low = 0
bad_wells = 0
matrix_name = "Diploid_carbo_06_Raffinose_MTX_T14.csv"
TECAN_data_AUC(data_file_name, map_name, time_window, time_window_low, matrix_name, bad_wells) #0 from 30

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190820/")
data_file_name = "2019-08-18_16-01_Diploid_carbo_07_Raffinose_MTX_96_T2.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_7.csv" 
time_window = 48
time_window_low = 10
bad_wells = 0
matrix_name = "Diploid_carbo_07_Raffinose_MTX_T2.csv"
TECAN_data_AUC(data_file_name, map_name, time_window,time_window_low, matrix_name, bad_wells) #3 from 30

##### Combine two environments
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/")
dip_01_SD = dataFrameReader_T("20190819/Diploid_carbo_01_SD_MTX_T7.csv")
dip_06_SD = dataFrameReader_T("20190819/Diploid_carbo_06_SD_MTX_T14.csv")
dip_07_SD = dataFrameReader_T("20190819/Diploid_carbo_07_SD_MTX_T2.csv")

dip_01_Raff = dataFrameReader_T("20190820/Diploid_carbo_01_Raffinose_MTX_T7.csv")
dip_06_Raff = dataFrameReader_T("20190820/Diploid_carbo_06_Raffinose_MTX_T14.csv")
dip_07_Raff = dataFrameReader_T("20190820/Diploid_carbo_07_Raffinose_MTX_T2.csv")

dip_SD = rbind(dip_01_SD, dip_06_SD, dip_07_SD)
dip_Raff = rbind(dip_01_Raff, dip_06_Raff, dip_07_Raff)

dip_SD_Raff = matrix(0, nrow(dip_SD), 7)
colnames(dip_SD_Raff) = c("PPI", "Dynamics_01", "Dynamics_02", "Dynamics_03",
                          "Mean_dynamics", "SD_dynamics", "p-value")
dip_SD_Raff[,1] = as.character(dip_SD[,1])
for(i in 1:nrow(dip_SD_Raff)){
        dynamic = as.numeric(dip_Raff[i,5:7]) - as.numeric(dip_SD[i,5:7])
        dip_SD_Raff[i, 2:4] = dynamic
        dip_SD_Raff[i,5] = mean(dynamic)
        dip_SD_Raff[i,6] = sd(dynamic)
        dip_SD_Raff[i,7] = t.test(as.numeric(dip_Raff[i,5:7]), as.numeric(dip_SD[i,5:7]),
                                  alternative = "two.sided")$p.value
}

csvWriter(dip_SD_Raff, "20190820/Raffinose_DMSO_dynamics_Tecan.csv")

##### Get fitness values from PPiSeq
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190820/")
PPI_fit = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
Tecan_SD_NaCl = csvReader_T("Raffinose_DMSO_dynamics_Tecan.csv")
PPI_fit_tecan = PPI_fit[which(PPI_fit[,1] %in% Tecan_SD_NaCl[,1]),]
PPI_NaCl_SD = as.numeric(PPI_fit_tecan[,9]) - as.numeric(PPI_fit_tecan[,4])
PPI_NaCl_SD_fit = cbind(PPI_fit_tecan[,1], PPI_NaCl_SD)
PPI_NaCl_SD_Tecan = Tecan_SD_NaCl[match(PPI_NaCl_SD_fit[,1], Tecan_SD_NaCl[,1]),5]
PPI_NaCl_SD_fit_Tecan = cbind(PPI_NaCl_SD_fit, PPI_NaCl_SD_Tecan)

PPI_NaCl_SD_fit_Tecan = PPI_NaCl_SD_fit_Tecan[which(PPI_NaCl_SD_fit_Tecan[,1] != "YBL029C-A_YHR096C"),]
csvWriter(PPI_NaCl_SD_fit_Tecan, "Raffinose_DMSO_Tecan_PPiSeq_comparison.csv")

Fit_comp = csvReader_T("Raffinose_DMSO_Tecan_PPiSeq_comparison.csv")
cor(as.numeric(Fit_comp[,2]), as.numeric(Fit_comp[,3])) # 0.5345029

PPI_count= csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
PPI_SD_pos_NaCl_neg = PPI_count[which(PPI_count[,3] == "1" & PPI_count[,8] == "0"),]
PPI_SD_pos_NaCl_pos = PPI_count[which(PPI_count[,3] == "1" & PPI_count[,8] == "1"),]
PPI_SD_neg_NaCl_pos = PPI_count[which(PPI_count[,3] == "0" & PPI_count[,8] == "1"),]
PPI_SD_neg_NaCl_neg = PPI_count[which(PPI_count[,3] == "0" & PPI_count[,8] == "0"),]

Fit_comp_01 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_neg),] # 43
nrow(Fit_comp_01)
Fit_comp_02 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_pos),] # 9
nrow(Fit_comp_02)
Fit_comp_03 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_pos),] # 7
nrow(Fit_comp_03)
Fit_comp_04 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_neg),] # 23
nrow(Fit_comp_04)

#Fitness_comparison = Fitness_comparison[which(Fitness_comparison[,1] %in% PPI_chosen[,1]),]
cor(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3])) # 0.3173504
cor(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3])) # 0.9098217
cor(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3])) # 0.4000961
cor(as.numeric(Fit_comp[,2]), as.numeric(Fit_comp[,3])) # 0.4000961

pdf("Raffinose_DMSO_Tecan_PPiSeq.pdf", width = 5, height = 5)
plot(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3]), 
     xlab = " Raffinose - SD by PPiSeq",
     ylab = " Raffinose - SD by Tecan", 
     pch = 16, col = apple_colors[1],  xlim = c(-1,0.5), ylim = c(-2.5, 1.5))
points(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3]), pch = 16, col = apple_colors[2])
points(as.numeric(Fit_comp_03[,2]), as.numeric(Fit_comp_03[,3]), pch = 16, col = apple_colors[6])
points(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3]), pch = 16, col = apple_colors[4])
legend("topleft", legend= c("SD + Raffinose -", "SD + Raffinose +", "SD - Raffinose +","SD - Raffinose -"), 
       col = apple_colors[c(1,2,6,4)],
       pch = c(16, 16, 16, 16), bty = "n")
dev.off()

a = Fit_comp[which(as.numeric(Fit_comp[,2]) > 0 & as.numeric(Fit_comp[,3]) > 0),]
a[order(as.numeric(a[,3]), decreasing = T),]
intersect(a[,1], c(Fit_comp_02[,1], Fit_comp_03[,1]))
### Only consider PPIs containing HXT5
Fit_comp_HXT5 = Fit_comp[grep("YHR096C", Fit_comp[,1]),] # 17
## Make a barplot to show HXT5 increasing
#Fit_comp_chosen = rbind(Fit_comp_01, Fit_comp_02)
PPI_HXT5 = t(Fit_comp_HXT5)
PPI_HXT5= rbind(as.numeric(PPI_HXT5[2,]), as.numeric(PPI_HXT5[3,]))
PPI_HXT5 = PPI_HXT5[,order(PPI_HXT5[1,], decreasing = T)]
colnames(PPI_HXT5) = Fit_comp[,1]
pdf("barplot_HXT5_Raffinose_SD.pdf", height =5, width = 5)
barplot(PPI_HXT5,horiz=F, beside=T, axisnames = F, ylab = "Raffinose - SD", 
        #ylim = c(-0.5, 0.8),
        main = "HXT5 containing PPIs", col = apple_colors[c(5,7)])
legend("topright", legend=c("PPiSeq", "TECAN"),
       fill=apple_colors[c(5,7)], cex=0.8, bty="n", border=FALSE)
dev.off()

Fit_comp_HXT1 = Fit_comp[grep("YHR094C", Fit_comp[,1]),] # 17
PPI_HXT1 = t(Fit_comp_HXT1)
PPI_HXT1= rbind(as.numeric(PPI_HXT1[2,]), as.numeric(PPI_HXT1[3,]))
PPI_HXT1 = PPI_HXT1[,order(PPI_HXT1[1,], decreasing = T)]
colnames(PPI_HXT1) = Fit_comp_HXT1[,1]
pdf("barplot_HXT1_Raffinose_SD.pdf", height =5, width = 5)
barplot(PPI_HXT1,horiz=F, beside=T, axisnames = F, ylab = "Raffinose - SD", 
        #ylim = c(-0.5, 0.8),
        main = "HXT1 containing PPIs", col = apple_colors[c(5,7)])
legend("topright", legend=c("PPiSeq", "TECAN"),
       fill=apple_colors[c(5,7)], cex=0.8, bty="n", border=FALSE)
dev.off()

Fit_comp_HXT1 = Fit_comp[grep("YDR345C", Fit_comp[,1]),] # 17
PPI_HXT1 = t(Fit_comp_HXT1)
PPI_HXT1= rbind(as.numeric(PPI_HXT1[2,]), as.numeric(PPI_HXT1[3,]))
PPI_HXT1 = PPI_HXT1[,order(PPI_HXT1[1,], decreasing = T)]
colnames(PPI_HXT1) = Fit_comp_HXT1[,1]
pdf("barplot_HXT3_Raffinose_SD.pdf", height =5, width = 5)
barplot(PPI_HXT1,horiz=F, beside=T, axisnames = F, ylab = "Raffinose - SD", 
        #ylim = c(-0.5, 0.8),
        main = "HXT3 containing PPIs", col = apple_colors[c(5,7)])
legend("topright", legend=c("PPiSeq", "TECAN"),
       fill=apple_colors[c(5,7)], cex=0.8, bty="n", border=FALSE)
dev.off()

Fit_comp_HXT1 = Fit_comp[grep("YDR342C", Fit_comp[,1]),] # 17
PPI_HXT1 = t(Fit_comp_HXT1)
PPI_HXT1= rbind(as.numeric(PPI_HXT1[2,]), as.numeric(PPI_HXT1[3,]))
PPI_HXT1 = PPI_HXT1[,order(PPI_HXT1[1,], decreasing = T)]
colnames(PPI_HXT1) = Fit_comp_HXT1[,1]
pdf("barplot_HXT7_Raffinose_SD.pdf", height =5, width = 5)
barplot(PPI_HXT1,horiz=F, beside=T, axisnames = F, ylab = "Raffinose - SD", 
        #ylim = c(-0.5, 0.8),
        main = "HXT7 containing PPIs", col = apple_colors[c(5,7)])
legend("topright", legend=c("PPiSeq", "TECAN"),
       fill=apple_colors[c(5,7)], cex=0.8, bty="n", border=FALSE)
dev.off()



