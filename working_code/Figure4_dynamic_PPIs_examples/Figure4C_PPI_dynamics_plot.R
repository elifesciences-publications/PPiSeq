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
        #Max_AUC = max(Mean_AUC)
        for(i in 1:length(PPI_unique)){
                column_each_PPI = which(PPI_map == PPI_unique[i])
                auc_each_PPI_01 = AUC_PPI(time_select_01, od_robot_select_01, column_each_PPI)
                auc_PPI_matrix[i,5:7] = (auc_each_PPI_01-Neg_AUC)/Neg_AUC
                auc_PPI_matrix[i,8] = mean((auc_each_PPI_01-Neg_AUC)/Neg_AUC)
                auc_PPI_matrix[i,9] = sd((auc_each_PPI_01-Neg_AUC)/Neg_AUC)
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

##### NaCl environment
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-16_14-01_Diploid_carbo_NaCl_MTX_01_96_T7.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_1.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "Diploid_carbo_01_NaCl_MTX_T7.csv"
TECAN_data_AUC(data_file_name, map_name, time_window,time_window_low, matrix_name, bad_wells) #11 from 30

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-16_15-02_Diploid_carbo_06_NaCl_MTX_96_T14.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_6.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "Diploid_carbo_06_NaCl_MTX_T14.csv"
TECAN_data_AUC(data_file_name, map_name, time_window,time_window_low, matrix_name, bad_wells) #11 from 30

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-16_15-03_Diploid_carbo_07_NaCl_MTX_96_T2.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_7.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "Diploid_carbo_07_NaCl_MTX_T2.csv"
TECAN_data_AUC(data_file_name, map_name, time_window,time_window_low, matrix_name, bad_wells) #11 from 30

#### SD environment
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-12_11-53_Diploid_carbo_01_DMSO_MTX_T7_96_T7.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_1.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "Diploid_carbo_01_SD_MTX_T7.csv"
TECAN_data_AUC(data_file_name, map_name, time_window,time_window_low, matrix_name, bad_wells) #12 from 30

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-12_11-53_Diploid_carbo_06_DMSO_MTX_96_T14.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_6.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "Diploid_carbo_06_SD_MTX_T14.csv"
TECAN_data_AUC(data_file_name, map_name, time_window,time_window_low, matrix_name, bad_wells) #13 from 30

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
data_file_name = "2019-08-12_11-55_Diploid_carbo_07_DMSO_MTX_96_T2.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/Final_map/Diploid_carbohydrate_7.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "Diploid_carbo_07_SD_MTX_T2.csv"
TECAN_data_AUC(data_file_name, map_name, time_window,time_window_low, matrix_name, bad_wells) #4 from 30

##### Combine two environments SD and Raffinose; SD and NaCl
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/")
dip_01_SD = dataFrameReader_T("20190819/Diploid_carbo_01_SD_MTX_T7.csv")
dip_06_SD = dataFrameReader_T("20190819/Diploid_carbo_06_SD_MTX_T14.csv")
dip_07_SD = dataFrameReader_T("20190819/Diploid_carbo_07_SD_MTX_T2.csv")

dip_01_Raff = dataFrameReader_T("20190820/Diploid_carbo_01_Raffinose_MTX_T7.csv")
dip_06_Raff = dataFrameReader_T("20190820/Diploid_carbo_06_Raffinose_MTX_T14.csv")
dip_07_Raff = dataFrameReader_T("20190820/Diploid_carbo_07_Raffinose_MTX_T2.csv")

dip_01_NaCl = dataFrameReader_T("20190819/Diploid_carbo_01_NaCl_MTX_T7.csv")
dip_06_NaCl = dataFrameReader_T("20190819/Diploid_carbo_06_NaCl_MTX_T14.csv")
dip_07_NaCl = dataFrameReader_T("20190819/Diploid_carbo_07_NaCl_MTX_T2.csv")

dip_SD = rbind(dip_01_SD, dip_06_SD, dip_07_SD)
dip_Raff = rbind(dip_01_Raff, dip_06_Raff, dip_07_Raff)
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
        dip_SD_NaCl[i,7] = t.test(as.numeric(dip_NaCl[i,5:7]), as.numeric(dip_NaCl[i,5:7]),
                                  alternative = "two.sided")$p.value
}

csvWriter(dip_SD_NaCl, "20190819/NaCl_DMSO_dynamics_Tecan.csv")


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

##### Get normalized fitness values from PPiSeq and calculate the mean dynamics, sd, and p-value
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190820/")
PPI_fit_SD = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/SD_PPI_barcodes_fitness_counts.csv")
PPI_fit_Raff = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/Raffinose_PPI_barcodes_fitness_counts.csv")
PPI_fit_NaCl = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/NaCl_PPI_barcodes_fitness_counts.csv")
Tecan = csvReader_T("Raffinose_DMSO_dynamics_Tecan.csv")

SD_chosen = PPI_fit_SD[which(PPI_fit_SD[,1] %in% c(Tecan[,1], "positive_DHFR","negative_non_DHFR")),]
Raff_chosen = PPI_fit_Raff[which(PPI_fit_Raff[,1] %in% c(Tecan[,1], "positive_DHFR","negative_non_DHFR")),]
NaCl_chosen = PPI_fit_NaCl[which(PPI_fit_NaCl[,1] %in% c(Tecan[,1], "positive_DHFR","negative_non_DHFR")),]

Fitness_normalization = function(DMSO_lineage){
        DMSO_DHFR_Pos= DMSO_lineage[which(DMSO_lineage[,1] == "positive_DHFR"),]
        DMSO_DHFR_Neg= DMSO_lineage[which(DMSO_lineage[,1] == "negative_non_DHFR"),]
        
        DMSO_DHFR_Pos_mean = mean(as.numeric(DMSO_DHFR_Pos[,4])) 
        DMSO_DHFR_Neg_mean = mean(as.numeric(DMSO_DHFR_Neg[,4])) 
        
        DMSO_lineage_chosen = DMSO_lineage[, c(1,2,3,4)]
        DMSO_lineage_chosen[,4] = (as.numeric(DMSO_lineage_chosen[,4]) - DMSO_DHFR_Neg_mean)/(DMSO_DHFR_Pos_mean - DMSO_DHFR_Neg_mean)
        return(DMSO_lineage_chosen)
}
SD_norm = Fitness_normalization(SD_chosen)
Raff_norm = Fitness_normalization(Raff_chosen)
NaCl_norm = Fitness_normalization(NaCl_chosen)
csvWriter(SD_norm, "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/PPiseq_norm_SD_Tecan.csv")
csvWriter(NaCl_norm, "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/PPiseq_norm_NaCl_Tecan.csv")
csvWriter(Raff_norm, "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190820/PPiseq_norm_Raffinose_Tecan.csv")

#### align the PPI score according to their barcodes
#SD_norm_overlap = SD_norm[which(SD_norm[,1] %in% Overlap),]
#Raff_norm_overlap = Raff_norm[which(Raff_norm[,1] %in% Overlap),]
#Raff_matched = Raff_norm[match(SD_norm_overlap[,3], Raff_norm_overlap[,3]),4]
#SD_Raff_norm = cbind(SD_norm_overlap, Raff_matched)
#colnames(SD_Raff_norm) = c(colnames(SD_norm_overlap), "Fit_Raff")
#csvWriter(SD_Raff_norm, "SD_Raff_normalized_Fitness.csv")
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/")
SD_norm = csvReader_T("20190819/PPiseq_norm_SD_Tecan.csv")
NaCl_norm = csvReader_T("20190819/PPiseq_norm_NaCl_Tecan.csv")
Raff_norm = csvReader_T("20190820/PPiseq_norm_Raffinose_Tecan.csv")

### SD an Raffinose
Overlap = intersect(SD_norm[,1], Raff_norm[,1]) # 82
Overlap = Overlap[which(!Overlap %in% c("positive_DHFR", "negative_non_DHFR"))]
SD_Raff_matrix = matrix(NA, length(Overlap), 6)
colnames(SD_Raff_matrix) = c("PPI", "BC_1", "BC_2", 
                             "Mean_diff", "SD_diff", "p-value")
SD_Raff_matrix[,1] = Overlap
for(i in 1:length(Overlap)){
        a = which(SD_norm[,1] == Overlap[i])
        b = which(Raff_norm[,1] == Overlap[i])
        SD_Raff_matrix[i,2] = length(a)
        SD_Raff_matrix[i,3] = length(b)
        fit_SD = as.numeric(SD_norm[a, 4])
        fit_Raff = as.numeric(Raff_norm[b,4])
        dynamics = rep(0, length(a) * length(b))
        z = 1
        for(m in 1:length(b)){
                for(n in 1:length(a)){
                        dynamics[z] = fit_Raff[m] - fit_SD[n]
                        z = z + 1
                }
        }
        SD_Raff_matrix[i,4] = mean(dynamics)
        SD_Raff_matrix[i,5] = sd(dynamics)
        if(length(a) > 1 & length(b) > 1){
                SD_Raff_matrix[i,6] = t.test(fit_SD, fit_Raff, alternative = "two.sided")$p.value
        }else{
                SD_Raff_matrix[i,6] = NA
        }
        
}
csvWriter(SD_Raff_matrix, "20190820/SD_Raff_normalized_PPiSeq_Fitness_summary.csv")

Tecan_matched = Tecan[match(SD_Raff_matrix[,1], Tecan[,1]),]
SD_Raff_final = cbind(Tecan_matched, SD_Raff_matrix) # Combine Tecan data and PPiSeq data
csvWriter(SD_Raff_final, "20190820/SD_Raff_Tecan_PPiSeq_all_comparison.csv")

### SD an NaCl
Overlap = intersect(SD_norm[,1], NaCl_norm[,1]) # 82
Overlap = Overlap[which(!Overlap %in% c("positive_DHFR", "negative_non_DHFR"))]
SD_NaCl_matrix = matrix(NA, length(Overlap), 6)
colnames(SD_NaCl_matrix) = c("PPI", "BC_1", "BC_2", 
                             "Mean_diff", "SD_diff", "p-value")
SD_NaCl_matrix[,1] = Overlap
for(i in 1:length(Overlap)){
        a = which(SD_norm[,1] == Overlap[i])
        b = which(NaCl_norm[,1] == Overlap[i])
        SD_NaCl_matrix[i,2] = length(a)
        SD_NaCl_matrix[i,3] = length(b)
        fit_SD = as.numeric(SD_norm[a, 4])
        fit_NaCl = as.numeric(NaCl_norm[b,4])
        dynamics = rep(0, length(a) * length(b))
        z = 1
        for(m in 1:length(b)){
                for(n in 1:length(a)){
                        dynamics[z] = fit_NaCl[m] - fit_SD[n]
                        z = z + 1
                }
        }
        SD_NaCl_matrix[i,4] = mean(dynamics)
        SD_NaCl_matrix[i,5] = sd(dynamics)
        if(length(a) > 1 & length(b) > 1){
                SD_NaCl_matrix[i,6] = t.test(fit_SD, fit_NaCl, alternative = "two.sided")$p.value
        }else{
                SD_NaCl_matrix[i,6] = NA
        }
        
}
csvWriter(SD_NaCl_matrix, "20190819/SD_NaCl_normalized_PPiSeq_Fitness_summary.csv")

#### Combine Tecan data and PPiSeq data
Tecan = csvReader_T("20190819/NaCl_DMSO_dynamics_Tecan.csv")
Tecan_matched = Tecan[match(SD_NaCl_matrix[,1], Tecan[,1]),]
SD_NaCl_final = cbind(Tecan_matched, SD_NaCl_matrix)
csvWriter(SD_NaCl_final, "20190819/SD_NaCl_Tecan_PPiSeq_all_comparison.csv")

############### Make a scatter plot to show the comparison between Tecan and PPiSeq
##### Separate each group by HXT1, HXT3, HXT5, and HXT7
##### Add error bars onto each spot (x-axis, y-axis)
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190820/")
SD_Raff_final = dataFrameReader_T("SD_Raff_Tecan_PPiSeq_all_comparison.csv")
name_exchange = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Systematic_standard_protein.csv")
PPI_split = split_string_vector(SD_Raff_final[,1])
protein_1 = name_exchange[match(PPI_split[,1], name_exchange[,1]),2]
protein_2 = name_exchange[match(PPI_split[,2], name_exchange[,1]),2]
label = rep("0", nrow(SD_Raff_final))
for(i in 1:length(protein_1)){
        if(protein_1[i] == "HXT1" | protein_2[i] == "HXT1"){
                label[i] = "HXT1"
        }
        else if(protein_1[i] == "HXT5" | protein_2[i] == "HXT5"){
                label[i] = "HXT5"
        }
        else if(protein_1[i] == "HXT3" | protein_2[i] == "HXT3"){
                label[i] = "HXT3"
        }
        else if(protein_1[i] == "HXT7" | protein_2[i] == "HXT7"){
                label[i] = "HXT7"
        }else if(protein_1[i] == "FPS1" | protein_2[i] == "FPS1"){
                label[i] = "FPS1"
        }else if(protein_1[i] == "HXT2" | protein_2[i] == "HXT2"){
                label[i] = "HXT2"
        }else{
                label[i] = "Other"
        }
}
SD_Raff_final = data.frame(protein_1, protein_2, SD_Raff_final, label)
P_value = rep("0", nrow(SD_Raff_final))
P_value[is.na(SD_Raff_final$p.value.1)] = "NA"
P_value[SD_Raff_final$p.value.1 <= 0.05] = "Significant"
P_value[SD_Raff_final$p.value.1 > 0.05] = "Non-significant"
SD_Raff_final$p.value.1 = P_value
#SD_Raff_final$p.value.1[which(SD_Raff_final$p.value.1 == 0)] = 1e-7
#SD_Raff_final$p.value.1 = -log10(SD_Raff_final$p.value.1)
PPiseq_min = SD_Raff_final$Mean_diff - SD_Raff_final$SD_diff
PPiseq_max = SD_Raff_final$Mean_diff + SD_Raff_final$SD_diff
Tecan_min = SD_Raff_final$Mean_dynamics - SD_Raff_final$SD_dynamics
Tecan_max = SD_Raff_final$Mean_dynamics + SD_Raff_final$SD_dynamics
SD_Raff_final = data.frame(SD_Raff_final, PPiseq_min, PPiseq_max, Tecan_min, Tecan_max)
SD_Raff_final$label = factor(SD_Raff_final$label, levels = c("HXT1", "HXT7", "HXT3", "HXT5", "HXT2", "FPS1", "Other"))
SD_Raff_final$p.value.1 = factor(SD_Raff_final$p.value.1, levels = c("Significant", "Non-significant", "NA"))
library(ggplot2)
#"HXT1" = "#7b3294", "HXT3" = "#c2a5cf", "HXT5" = "#d01c8b", "HXT7" = "#a6dba0", "Others" ="#008837"
cor(SD_Raff_final$Mean_diff, SD_Raff_final$Mean_dynamics, method = "spearman") #0.6020863
ggplot(data = SD_Raff_final, aes(x = Mean_diff, y = Mean_dynamics))+ 
        geom_point(aes(shape= p.value.1, col = label), size =3)+
        geom_errorbarh(aes(xmin = PPiseq_min, xmax = PPiseq_max), col = apple_colors[8], size = 0.2 )+
        geom_errorbar(aes(ymin = Tecan_min, ymax = Tecan_max), col = apple_colors[8], size = 0.2)+
        geom_vline(xintercept = 0, col = apple_colors[11], linetype = 2,size = 0.2)+
        geom_hline(yintercept = 0, col = apple_colors[11], linetype = 2, size = 0.2)+
        annotate("text", x = -0.3, y = 0.7, label = expression(paste("Spearman's ", italic(r), " = 0.6")),  
                 parse = TRUE, col = apple_colors[11]) +
        scale_shape_manual(name = "", values = c(16, 15, 17)) +
        scale_color_manual(name = "", values = c("#1b9e77","#e7298a", "#d95f02", "#7570b3", "#1f78b4", "#984ea3", "#CECED2"))+
        scale_y_continuous(name = "Fitness change in Raffinose by OD595",
                           limits=c(-0.5, 0.8),
                           breaks=seq(-0.5,0.8, by =0.1),
                           labels = seq(-0.5,0.8, by= 0.1)) +
        scale_x_continuous(name = "Fitness change in Rraffinose by PPiSeq", 
                           limits=c(-0.5, 0.3),
                           breaks=seq(-0.5,0.3, by =0.1),
                           labels = seq(-0.5,0.3, by= 0.1))+
        theme(legend.key=element_blank(), legend.text=element_text(size=10)) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure4/Figure4C_Raffinose_SD_Tecan_PPiSeq_comparison.pdf", width =6, height =5 )       

################ NaCl Environment
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/carbohydrate_transport/data/20190819/")
SD_Raff_final = dataFrameReader_T("SD_NaCl_Tecan_PPiSeq_all_comparison.csv")
name_exchange = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Systematic_standard_protein.csv")
PPI_split = split_string_vector(SD_Raff_final[,1])
protein_1 = name_exchange[match(PPI_split[,1], name_exchange[,1]),2]
protein_2 = name_exchange[match(PPI_split[,2], name_exchange[,1]),2]
label = rep("0", nrow(SD_Raff_final))
for(i in 1:length(protein_1)){
        if(protein_1[i] == "HXT1" | protein_2[i] == "HXT1"){
                label[i] = "HXT1"
        }
        else if(protein_1[i] == "HXT5" | protein_2[i] == "HXT5"){
                label[i] = "HXT5"
        }
        else if(protein_1[i] == "HXT3" | protein_2[i] == "HXT3"){
                label[i] = "HXT3"
        }
        else if(protein_1[i] == "HXT7" | protein_2[i] == "HXT7"){
                label[i] = "HXT7"
        }else if(protein_1[i] == "FPS1" | protein_2[i] == "FPS1"){
                label[i] = "FPS1"
        }else if(protein_1[i] == "HXT2" | protein_2[i] == "HXT2"){
                label[i] = "HXT2"
        }else{
                label[i] = "Other"
        }
}
SD_Raff_final = data.frame(protein_1, protein_2, SD_Raff_final, label)
P_value = rep("0", nrow(SD_Raff_final))
P_value[is.na(SD_Raff_final$p.value.1)] = "NA"
P_value[SD_Raff_final$p.value.1 <= 0.05] = "Significant"
P_value[SD_Raff_final$p.value.1 > 0.05] = "Non-significant"
SD_Raff_final$p.value.1 = P_value
#SD_Raff_final$p.value.1[which(SD_Raff_final$p.value.1 == 0)] = 1e-7
#SD_Raff_final$p.value.1 = -log10(SD_Raff_final$p.value.1)
PPiseq_min = SD_Raff_final$Mean_diff - SD_Raff_final$SD_diff
PPiseq_max = SD_Raff_final$Mean_diff + SD_Raff_final$SD_diff
Tecan_min = SD_Raff_final$Mean_dynamics - SD_Raff_final$SD_dynamics
Tecan_max = SD_Raff_final$Mean_dynamics + SD_Raff_final$SD_dynamics
SD_Raff_final = data.frame(SD_Raff_final, PPiseq_min, PPiseq_max, Tecan_min, Tecan_max)
SD_Raff_final$label = factor(SD_Raff_final$label, levels = c("HXT1", "HXT7", "HXT3", "HXT5", "HXT2", "FPS1", "Other"))
SD_Raff_final$p.value.1 = factor(SD_Raff_final$p.value.1, levels = c("Significant", "Non-significant", "NA"))
library(ggplot2)
#"HXT1" = "#7b3294", "HXT3" = "#c2a5cf", "HXT5" = "#d01c8b", "HXT7" = "#a6dba0", "Others" ="#008837"
cor(SD_Raff_final$Mean_diff, SD_Raff_final$Mean_dynamics, method = "spearman") #0.1559781
ggplot(data = SD_Raff_final, aes(x = Mean_diff, y = Mean_dynamics))+ 
        geom_point(aes(shape= p.value.1, col = label), size = 3)+
        geom_errorbarh(aes(xmin = PPiseq_min, xmax = PPiseq_max), col = apple_colors[8], size = 0.2 )+
        geom_errorbar(aes(ymin = Tecan_min, ymax = Tecan_max), col = apple_colors[8], size = 0.2)+
        geom_vline(xintercept = 0, col = apple_colors[11], linetype = 2,size = 0.2)+
        geom_hline(yintercept = 0, col = apple_colors[11], linetype = 2, size = 0.2)+
        annotate("text", x = -0.3, y = 0.2, label = expression(paste("Spearman's ", italic(r), " = 0.16")),  
                 parse = TRUE, col = apple_colors[11]) +
        scale_shape_manual(name = "", values = c(16, 15, 17)) +
        scale_color_manual(name = "", values = c("#1b9e77","#e7298a", "#d95f02", "#7570b3", "#1f78b4", "#984ea3", "#CECED2"))+
        scale_y_continuous(name = "Fitness change in Raffinose by OD595",
                           limits=c(-0.4, 0.3),
                           breaks=seq(-0.4,0.3, by =0.1),
                           labels = seq(-0.4,0.3, by= 0.1)) +
        scale_x_continuous(name = "Fitness change in Rraffinose by PPiSeq", 
                           limits=c(-0.5, 0.3),
                           breaks=seq(-0.5,0.3, by =0.1),
                           labels = seq(-0.5,0.3, by= 0.1))+
        theme(legend.key=element_blank(), legend.text=element_text(size=10)) +
        #guides(fill=guide_legend(title="Log10(Count)")) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure4/Figure4C_NaCl_SD_Tecan_PPiSeq_comparison.pdf", width =6, height =5 )       




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





