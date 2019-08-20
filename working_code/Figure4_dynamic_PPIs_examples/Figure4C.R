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
tecan_1 = csvReader_T("Working_data/TECAN_validation/carbohydrate_transport/Diploid_carbohydrate_96well_1.csv")
tecan_6 = csvReader_T("Working_data/TECAN_validation/carbohydrate_transport/Diploid_carbohydrate_96well_6.csv")
tecan_7 = csvReader_T("Working_data/TECAN_validation/carbohydrate_transport/Diploid_carbohydrate_96well_7.csv")
#tecan_7= tecan_7[which(tecan_7[,1] != "YBL029C-A_YHR096C"),] # Lose this PPI during the growth
PPI_tecan = unique(c(tecan_1[,1], tecan_6[,1], tecan_7[,1])) #
PPI_fit = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
PPI_fit_tecan = PPI_fit[which(PPI_fit[,1] %in% PPI_tecan),]
PPI_NaCl_SD = as.numeric(PPI_fit_tecan[,10]) - as.numeric(PPI_fit_tecan[,4])
PPI_Raffinose_SD = as.numeric(PPI_fit_tecan[,9]) - as.numeric(PPI_fit_tecan[,4])
PPI_NaCl_SD_fit = cbind(PPI_fit_tecan[,1], PPI_NaCl_SD)
PPI_Raffinose_SD_fit = cbind(PPI_fit_tecan[,1], PPI_Raffinose_SD)

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
        auc_PPI_matrix = matrix(0, length(PPI_unique), 9)
        colnames(auc_PPI_matrix) = c("PPI","Mean_Neg_AUC", "SD_Neg_AUC","AUC_01", "AUC_02", "AUC_03", "Mean_AUC", "SD_AUC", "P-value")
        auc_PPI_matrix[,1] = as.character(PPI_unique)
        auc_PPI_matrix[,2] = mean(Negative_AUC_01)
        auc_PPI_matrix[,3] = sd(Negative_AUC_01)
        for(i in 1:length(PPI_unique)){
                column_each_PPI = which(PPI_map == PPI_unique[i])
                auc_each_PPI_01 = AUC_PPI(time_select_01, od_robot_select_01, column_each_PPI)
                auc_PPI_matrix[i,4:6] = auc_each_PPI_01- mean(Negative_AUC_01)
                auc_PPI_matrix[i,7] = mean(auc_each_PPI_01- mean(Negative_AUC_01))
                auc_PPI_matrix[i,8] = sd(auc_each_PPI_01 - mean(Negative_AUC_01))
                auc_PPI_matrix[i,9] = t.test(auc_each_PPI_01, Negative_AUC_01, alternative = "greater")$p.value
        }
        auc_PPI_matrix[,9] = p.adjust(auc_PPI_matrix[,9], "BH")
        write.csv(auc_PPI_matrix, matrix_name, row.names = F)
        #return(length(which(as.numeric(auc_PPI_matrix[,6]) <= 0.05 & as.numeric(auc_PPI_matrix[,5]) >= 0)))
        return(length(which(as.numeric(auc_PPI_matrix[,9]) <= 0.05)))
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
dip_01_SD = dataFrameReader_T("Diploid_carbo_01_SD_MTX_T7.csv")
dip_01_NaCl = dataFrameReader_T("Diploid_carbo_01_NaCl_MTX_T7.csv")
dip_06_SD = dataFrameReader_T("Diploid_carbo_06_SD_MTX_T14.csv")
dip_06_NaCl = dataFrameReader_T("Diploid_carbo_06_NaCl_MTX_T14.csv")
dip_07_SD = dataFrameReader_T("Diploid_carbo_07_SD_MTX_T2.csv")
dip_07_NaCl = dataFrameReader_T("Diploid_carbo_07_NaCl_MTX_T2.csv")

dip_01_NaCl_SD = dip_01_NaCl$Mean_AUC - dip_01_SD$Mean_AUC
dip_06_NaCl_SD = dip_06_NaCl$Mean_AUC - dip_06_SD$Mean_AUC
dip_07_NaCl_SD = dip_07_NaCl$Mean_AUC - dip_07_SD$Mean_AUC

NaCl_SD = rbind(cbind(as.character(dip_01_SD[,1]),dip_01_NaCl_SD),
                cbind(as.character(dip_06_SD[,1]), dip_06_NaCl_SD),
                cbind(as.character(dip_07_SD[,1]), dip_07_NaCl_SD))

#### Compare the two measurements (subtraction)
Fit_tecan = NaCl_SD[match(PPI_NaCl_SD_fit[,1], NaCl_SD[,1]),2]
Fitness_comparison = cbind(PPI_NaCl_SD_fit, Fit_tecan)
Fitness_comparison = Fitness_comparison[which(Fitness_comparison[,1] != "YBL029C-A_YHR096C"),]
colnames(Fitness_comparison) = c("PPI", "PPiSeq", "Tecan")
csvWriter(Fitness_comparison, "NaCl_DMSO_Tecan_PPiSeq_comparison.csv")

Fit_comp = csvReader_T("NaCl_DMSO_Tecan_PPiSeq_comparison.csv")
Fit_comp= Fit_comp[which(Fit_comp[,1] != "YBL029C-A_YHR096C"),]
PPI_count= csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
PPI_SD_pos_NaCl_neg = PPI_count[which(PPI_count[,3] == "1" & PPI_count[,9] == "0"),]
PPI_SD_pos_NaCl_pos = PPI_count[which(PPI_count[,3] == "1" & PPI_count[,9] == "1"),]
PPI_SD_neg_NaCl_pos = PPI_count[which(PPI_count[,3] == "0" & PPI_count[,9] == "1"),]
PPI_SD_neg_NaCl_neg = PPI_count[which(PPI_count[,3] == "0" & PPI_count[,9] == "0"),]

Fit_comp_01 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_neg),] # 44
Fit_comp_02 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_pos),] # 17
Fit_comp_03 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_pos),] # 0
Fit_comp_04 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_neg),] # 29

Fitness_comparison = Fitness_comparison[which(Fitness_comparison[,1] %in% PPI_chosen[,1]),]
cor(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3])) # 0.06328027
cor(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3])) # 0.5037215
cor(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3])) # 0.5085061
cor(as.numeric(Fit_comp[,2]), as.numeric(Fit_comp[,3])) # 0.1911889

pdf("NaC_DMSO_Tecan_PPiSeq.pdf", width = 5, height = 5)
plot(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3]), 
     xlab = " NaCl - SD by PPiSeq",
     ylab = " NaCl - SD by Tecan", 
     pch = 16, col = apple_colors[1], bty = "n", ylim = c(-250, 150), xlim = c(-0.4, 0.6))
points(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3]), pch = 16, col = apple_colors[2])
#points(as.numeric(Fit_comp_03[,2]), as.numeric(Fit_comp_03[,3]), pch = 16, col = apple_colors[3])
points(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3]), pch = 16, col = apple_colors[4])
legend("topleft", legend= c("SD + NaCl -", "SD + NaCl +", "SD - NaCl -"), col = apple_colors[c(1,2,4)],
       pch = c(16, 16, 16), bty = "n")
dev.off()


##### Combine two environments
dip_01_SD = dataFrameReader_T("Diploid_carbo_01_SD_MTX_T7.csv")
dip_01_NaCl = dataFrameReader_T("Diploid_carbo_01_NaCl_MTX_T7.csv")
dip_06_SD = dataFrameReader_T("Diploid_carbo_06_SD_MTX_T14.csv")
dip_06_NaCl = dataFrameReader_T("Diploid_carbo_06_NaCl_MTX_T14.csv")
dip_07_SD = dataFrameReader_T("Diploid_carbo_07_SD_MTX_T2.csv")
dip_07_NaCl = dataFrameReader_T("Diploid_carbo_07_NaCl_MTX_T2.csv")

dip_01_NaCl_SD = dip_01_NaCl$Mean_AUC / dip_01_SD$Mean_AUC
dip_06_NaCl_SD = dip_06_NaCl$Mean_AUC / dip_06_SD$Mean_AUC
dip_07_NaCl_SD = dip_07_NaCl$Mean_AUC / dip_07_SD$Mean_AUC

NaCl_SD = rbind(cbind(as.character(dip_01_SD[,1]),dip_01_NaCl_SD),
                cbind(as.character(dip_06_SD[,1]), dip_06_NaCl_SD),
                cbind(as.character(dip_07_SD[,1]), dip_07_NaCl_SD))

#### Compare the two measurements (subtraction)
PPI_fit = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
PPI_fit_tecan = PPI_fit[which(PPI_fit[,1] %in% PPI_tecan),]
PPI_NaCl_SD = as.numeric(PPI_fit_tecan[,10]) / as.numeric(PPI_fit_tecan[,4])
PPI_NaCl_SD_fit = cbind(PPI_fit_tecan[,1], PPI_NaCl_SD)

Fit_tecan = NaCl_SD[match(PPI_NaCl_SD_fit[,1], NaCl_SD[,1]),2]
Fitness_comparison = cbind(PPI_NaCl_SD_fit, Fit_tecan)
Fitness_comparison = Fitness_comparison[which(Fitness_comparison[,1] != "YBL029C-A_YHR096C"),]
colnames(Fitness_comparison) = c("PPI", "PPiSeq", "Tecan")
csvWriter(Fitness_comparison, "NaCl_DMSO_Tecan_PPiSeq_comparison_division.csv")

Fit_comp = csvReader_T("NaCl_DMSO_Tecan_PPiSeq_comparison_division.csv")
#Fit_comp= Fit_comp[which(Fit_comp[,1] != "YBL029C-A_YHR096C"),]
PPI_count= csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
PPI_SD_pos_NaCl_neg = PPI_count[which(PPI_count[,3] == "1" & PPI_count[,9] == "0"),]
PPI_SD_pos_NaCl_pos = PPI_count[which(PPI_count[,3] == "1" & PPI_count[,9] == "1"),]
PPI_SD_neg_NaCl_pos = PPI_count[which(PPI_count[,3] == "0" & PPI_count[,9] == "1"),]
PPI_SD_neg_NaCl_neg = PPI_count[which(PPI_count[,3] == "0" & PPI_count[,9] == "0"),]

Fit_comp_01 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_neg),] # 35
Fit_comp_02 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_pos_NaCl_pos),] # 17
Fit_comp_03 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_pos),] # 0
Fit_comp_04 = Fit_comp[which(Fit_comp[,1] %in% PPI_SD_neg_NaCl_neg),] # 29

cor(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3])) # 0.04489586
cor(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3])) # 0.5037215
cor(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3])) # 0.5169818


pdf("NaC_DMSO_Tecan_PPiSeq_division.pdf", width = 5, height = 5)
plot(as.numeric(Fit_comp_01[,2]), as.numeric(Fit_comp_01[,3]), 
     xlab = " NaCl - SD by PPiSeq",
     ylab = " NaCl - SD by Tecan", 
     pch = 16, col = apple_colors[1], bty = "n", ylim = c(-5,5), xlim = c(-4, 4))
points(as.numeric(Fit_comp_02[,2]), as.numeric(Fit_comp_02[,3]), pch = 16, col = apple_colors[2])
#points(as.numeric(Fit_comp_03[,2]), as.numeric(Fit_comp_03[,3]), pch = 16, col = apple_colors[3])
points(as.numeric(Fit_comp_04[,2]), as.numeric(Fit_comp_04[,3]), pch = 16, col = apple_colors[4])
legend("topleft", legend= c("SD + NaCl -", "SD + NaCl +", "SD - NaCl -"), col = apple_colors[c(1,2,4)],
       pch = c(16, 16, 16), bty = "n")
dev.off()