# Written by Zhimin
# modified by darach 190419
library(tidyverse)
library(pracma)

TECAN_data_AUC_correct <- function( data_file_name_MTX, data_file_name_non_MTX, 
                                    map_name, time_window_MTX, 
                                    time_window_no_MTX, 
                                    bad_wells) {
        # data_file_name_MTX: MTX TECAN data
        # data_file_name_non_MTX: non-MTX TECAN data
        # map_name: PPI name for each well
        # time_window_MTX: time point in MTX to calculate the AUC
        # time_window_no_MTX: time point in non-MTX to calculate the AUC
        # matrix_name: final output matrix name. Matrix contain all the information
        # bad wells: wells that have noisy data or empty data
#        
        ## robot_data_01 as MTX data
        robot_data_01= read.csv(data_file_name_MTX, skip = 19)
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
        ## robot_data_02 as non_MTX data
        robot_data_02= read.csv(data_file_name_non_MTX, skip = 19)
        robot_data_02 = robot_data_02[1:(nrow(robot_data_02)-2),]
        #Matix of ODs
        od_robot_02 = as.matrix(robot_data_02[, 3:ncol(robot_data_02)])
        od_robot_02 = od_robot_02[,which(!map_index %in% bad_wells)] # remove the data of bad wells
#        
        y = as.vector(robot_data_01[,1])
        z = strsplit(y, "\t")
        time_robot_02 = 1:length(z)
        for(i in 1:length(time_robot_02)){
                time_robot_02[i] = as.numeric(z[[i]][2])/60
        }
#        
        ## input the map 
        diploid_05_02 = read.csv(map_name)
        PPI_map = diploid_05_02[,1]
        PPI_map = PPI_map[which(!map_index %in% bad_wells)]
        PPI_unique = unique(PPI_map)
        PPI_unique = PPI_unique[which(PPI_unique != "Negative")]
#        
        AUC_PPI <- function(time, od_data, column_PPI){
                auc_PPI = rep(0, length(column_PPI))
                for(i in 1: length(column_PPI)){
                        auc_PPI[i] = trapz(time, od_data[,column_PPI[i]])
                }
                return(auc_PPI)
        }
#        
        ## Calculate AUC for negative wells in MTX and non_MTX 
        #MTX
        column_negative_PPI = which(PPI_map == "Negative")
        time_select_01 = time_robot_01[which(time_robot_01/60 <= time_window_MTX)] ## only consider data before MTX_time_window_threshold
        od_robot_select_01 = od_robot_01[1:length(time_select_01),] # only consider data before MTX_time_window_threshold
        Negative_AUC_01 = AUC_PPI(time_select_01, od_robot_select_01, column_negative_PPI)
        # non_MTX
        time_select_02 = time_robot_02[which(time_robot_02/60 <= time_window_no_MTX)] ## only consider data before non_MTX_time_window_threshold
        od_robot_select_02 = od_robot_02[1:length(time_select_02),] # only consider data before non_MTX_time_window_threshold
        Negative_AUC_02 = AUC_PPI(time_select_02, od_robot_select_02, column_negative_PPI)
        ## Subtract non_MTX_AUC from MTX_AUC
        Negative_AUC_final = Negative_AUC_01 - Negative_AUC_02
#        
        ### Create a matrix that contain mean of AUC_MTX - AUC_non_MTX for each PPI, and calculate one-tail Q value
        auc_PPI_matrix = matrix(0, length(PPI_unique), 12)
        colnames(auc_PPI_matrix) = c("PPI", "Mean_auc_MTX", "Mean_auc_non_MTX", 
            "Mean_auc_subtract", "SD_auc_subtract", "Mean_auc_neg_MTX", 
            "Mean_auc_neg_non_MTX","Mean_auc_neg_subtract", 
            "difference_subtract", "p-value", "one_tail_Q", "statistic")
        auc_PPI_matrix[,1] = as.character(PPI_unique)
        auc_PPI_matrix[,6] = mean(Negative_AUC_01)
        auc_PPI_matrix[,7] = mean(Negative_AUC_02)
        auc_PPI_matrix[,8] = mean(Negative_AUC_final)
        for(i in 1:length(PPI_unique)){
                column_each_PPI = which(PPI_map == PPI_unique[i])
                auc_each_PPI_01 = AUC_PPI(time_select_01, od_robot_select_01, column_each_PPI)
                auc_each_PPI_02 = AUC_PPI(time_select_02, od_robot_select_02, column_each_PPI)
                auc_each_PPI = auc_each_PPI_01 - auc_each_PPI_02
                auc_PPI_matrix[i,2] = mean(auc_each_PPI_01)
                auc_PPI_matrix[i,3] = mean(auc_each_PPI_02)
                auc_PPI_matrix[i,4] = mean(auc_each_PPI)
                auc_PPI_matrix[i,5] = sd(auc_each_PPI)
                auc_PPI_matrix[i,9] = as.numeric(auc_PPI_matrix[i,4]) - as.numeric(auc_PPI_matrix[i,8])
                auc_PPI_matrix[i,10] = t.test(auc_each_PPI, Negative_AUC_final, alternative = "greater")$p.value
                auc_PPI_matrix[i,11] = t.test(auc_each_PPI, Negative_AUC_final, alternative = "greater")$p.value
                auc_PPI_matrix[i,12] = t.test(auc_each_PPI, Negative_AUC_final, alternative = "greater")$statistic
        }
        auc_PPI_matrix[,11] = p.adjust(auc_PPI_matrix[,11], "BH")
        return(as.data.frame(auc_PPI_matrix))
        #write.csv(auc_PPI_matrix, matrix_name, row.names = F)
        #return(length(which(as.numeric(auc_PPI_matrix[,6]) <= 0.05 & as.numeric(auc_PPI_matrix[,5]) >= 0)))
        #return(length(which(as.numeric(auc_PPI_matrix[,11]) <= 0.05)))
}

list_to_proc <-       
    list(
data.frame(
    data_file_name_01 = "20190303_TECAN/2019-03-01_17-08_D_01_1_t5_96_T5.txt",
    data_file_name_02 =  "20190311_TECAN/2019-03-11_16-30_D_01_1_no_MTx_96_T5.txt",
    map_name = "Diploid_01_PPI_200_1.csv",
    time_window_01 = 40,
    time_window_02 = 25,
    bad_wells = 0
                        ),
data.frame(
#### Diploid_01_1 TECAN-T5 
## MTX data: 20190303_TECAN; non-MTX data: 20190311_TECAN
## Bad wells: MTX: no, non-MTX: no
data_file_name_01 = "20190303_TECAN/2019-03-01_17-08_D_01_1_t5_96_T5.txt",
data_file_name_02 =  "20190311_TECAN/2019-03-11_16-30_D_01_1_no_MTx_96_T5.txt",
map_name = "Diploid_01_PPI_200_1.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #15 from 30
data.frame(
#### Diploid_01_2 TECAN-T14 
## MTX data: 20190220_TECAN; non-MTX data: 20190310_TECAN
## Bad wells: MTX: no, non-MTX: no
data_file_name_01 = "20190220_TECAN/2019-02-18_14-33_Diploid-01-02-T14_96_T14.txt",
data_file_name_02 =  "20190310_TECAN/2019-03-10_09-13_D_01_2_no_MTX_96_T14.txt",
map_name = "Diploid_01_PPI_200_2.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #12 from 30
data.frame(
#### Diploid_01_3 TECAN-T15 
## MTX data: 20190220_TECAN; non-MTX data: 20190310_TECAN
## Bad wells: MTX: 52-54 one PPI, non-MTX: 52:54 one PPI
data_file_name_01 = "20190220_TECAN/2019-02-18_14-34_Diploid-01-03-T15_96_T15.txt",
data_file_name_02 =  "20190310_TECAN/2019-03-10_09-15_D_01_3_No_MTX_96_T15.txt",
map_name = "Diploid_01_PPI_200_3.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(c(1,13,25,37,49,52:54, 61,73,85))
                       ), #12 from 29
data.frame(
#### Diploid_01_4 TECAN-T2 
## MTX data: 20190220_TECAN; non-MTX data: 20190310_TECAN
## Bad wells: MTX: 94:96 negative control, non-MTX: 94:96 negative control
data_file_name_01 = "20190220_TECAN/2019-02-18_14-36_Diploid-01-04-T2_96_T2.txt",
data_file_name_02 =  "20190310_TECAN/2019-03-10_09-19_D_01_4_No_MTX_96_T2.txt",
map_name = "Diploid_01_PPI_200_4.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(94:96)
                       ), #12 from 30
data.frame(
#### Diploid_01_5 TECAN-T5 
## MTX data: 20190220_TECAN; non-MTX data: 20190310_TECAN
## Bad wells: MTX: 16:18 one PPI, non-MTX: 16:18 one PPI
data_file_name_01 = "20190220_TECAN/2019-02-18_14-38_Diploid-01-05-T5_96_T5.txt",
data_file_name_02 =  "20190310_TECAN/2019-03-10_09-20_D_01_5_No_MTX_96_T5.txt",
map_name = "Diploid_01_PPI_200_5.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(16:18)
                       ), #16 from 29
data.frame(
#### Diploid_01_6 TECAN-T7 
## MTX data: 20190222_TECAN; non-MTX data: 20190306_TECAN
## Bad wells: no, non-MTX: no
data_file_name_01 = "20190222_TECAN/2019-02-20_16-18_D_01_06_T7_96_T7.txt",
data_file_name_02 =  "20190306_TECAN/2019-03-04_11-50_D_01_6_no_MTX_t7_96_T7.txt",
map_name = "Diploid_01_PPI_200_6.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #10 from 30
data.frame(
#### Diploid_01_7 TECAN-T14
## MTX data: 20190222_TECAN; non-MTX data: 201903011_TECAN
## Bad wells:55:57 (one PPI), 64:93 (empty) non-MTX: 55:57 (one PPI), 64:93 (empty)
data_file_name_01 = "20190222_TECAN/2019-02-20_16-19_D_01_07_t14_96_T14.txt",
data_file_name_02 =  "20190311_TECAN/2019-03-11_16-26_D_01_7_No_MTX_96_T14.txt",
map_name = "Diploid_01_PPI_200_7.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(c(55:57, 64:93))
                       ), #12 from 19
data.frame(
#### Diploid_02_03_04_01 TECAN-T15
## MTX data: 20190222_TECAN; non-MTX data: 201903011_TECAN
## Bad wells:T15 non-MTX: T15
data_file_name_01 = "20190222_TECAN/2019-02-20_16-21_D_02_03_04_1_t15_96_T15.txt",
data_file_name_02 =  "20190311_TECAN/2019-03-11_16-28_D_02_03_04_1_no_MTX_96_T15.txt",
map_name = "Diploid_02_03_04_PPI_1.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(c(1,13,25,37,49, 61,73,85))
                       ), #24 from 30
data.frame(
#### Diploid_02_03_04_02 MTX TECAN -07 and non-MTXTECAN-T2
## MTX data: 20190222_TECAN; non-MTX data: 201903011_TECAN
## Bad wells:73 non-MTX: 73
data_file_name_01 = "20190316_TECAN/2019-03-16_09-48_D_02_03_04_2_MTX_96_T7.txt",
data_file_name_02 =  "20190311_TECAN/2019-03-11_16-29_D_02_03_04_2_No_MTX_96_T2.txt",
map_name = "Diploid_02_03_04_PPI_2.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(73)
                       ), #19 from 30
data.frame(
#### Diploid_02_03_04_03 TECAN-T5
## MTX data: 20190316_TECAN; non-MTX data: 201903013_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190222_TECAN/2019-02-20_16-24_D_02_03_04_3_t5_96_T5.txt",
data_file_name_02 =  "20190313_TECAN/2019-03-13_09-14_D_02_03_04_3_no_MTx_96_T5.txt",
map_name = "Diploid_02_03_04_PPI_3.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #26 from 30
data.frame(
#### Diploid_02_03_04_04 TECAN-T7
## MTX data: 20190227_TECAN; non-MTX data: 201903014_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190227_TECAN/2019-02-25_17-18_D_02_03_04_4_T7_96_T7.txt",
data_file_name_02 =  "20190314_TECAN/2019-03-14_16-06_D_02_03_04_4_no_MTX_96_T7.txt",
map_name = "Diploid_02_03_04_PPI_4.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #26 from 30
data.frame(
#### Diploid_02_03_04_05 TECAN-T14
## MTX data: 20190227_TECAN; non-MTX data: 20190306_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190227_TECAN/2019-02-25_17-19_D_02_03_04_5_t14_96_T14.txt",
data_file_name_02 =  "20190306_TECAN/2019-03-04_11-51_D_02_03_04_5_no_MTX_t14_96_T14.txt",
map_name = "Diploid_02_03_04_PPI_5.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #29 from 30
data.frame(
#### Diploid_02_03_04_06 TECAN-T15
## MTX data: 20190227_TECAN; non-MTX data: 20190313_TECAN
## Bad wells:T15 non-MTX: T15
data_file_name_01 = "20190227_TECAN/2019-02-25_17-20_D_02_03_04_6_t15_96_T15.txt",
data_file_name_02 =  "20190313_TECAN/2019-03-13_09-12_D_02_03_04_6_no_MTX_96_T15.txt",
map_name = "Diploid_02_03_04_PPI_6.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(c(1,13,25,37,49, 61,73,85))
                       ), #27 from 30
data.frame(
#### Diploid_02_03_04_07 TECAN-T2
## MTX data: 20190227_TECAN; non-MTX data: 20190313_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190227_TECAN/2019-02-25_17-21_D_02_03_04_7_t2_96_T2.txt",
data_file_name_02 =  "20190313_TECAN/2019-03-13_09-13_D_02_03_04_7_No_MTX_96_T2.txt",
map_name = "Diploid_02_03_04_PPI_7.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #25 from 30
data.frame(
#### Diploid_02_03_04_08 TECAN-T5
## MTX data: 20190227_TECAN; non-MTX data: 20190314_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190227_TECAN/2019-02-25_17-22_D_02_03_04_8_t5_96_T5.txt",
data_file_name_02 =  "20190314_TECAN/2019-03-14_16-12_D_02_03_04_8_no_MTX_96_T5.txt",
map_name = "Diploid_02_03_04_PPI_8.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #23 from 30
data.frame(
#### Diploid_02_03_04_09 TECAN-T7
## MTX data: 20190301_TECAN; non-MTX data: 20190311_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190301_TECAN/2019-02-27_17-07_D_02_03_04_9_T7_96_T7.txt",
data_file_name_02 =  "20190311_TECAN/2019-03-11_16-25_D_02_03_04_9_No_MTX_96_T7.txt",
map_name = "Diploid_02_03_04_PPI_9.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #28 from 30
data.frame(
#### Diploid_02_03_04_10 TECAN-T14
## MTX data: 20190301_TECAN; non-MTX data: 20190313_TECAN
## Bad wells:10:12 (one PPI) non-MTX: 10:12 (one PPI)
data_file_name_01 = "20190301_TECAN/2019-02-27_17-08_D_02_03_04_10_t14_96_T14.txt",
data_file_name_02 =  "20190313_TECAN/2019-03-13_09-10_D_02_03_04_10_No_MTX_96_T14.txt",
map_name = "Diploid_02_03_04_PPI_10.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(10:12)
                       ), #25 from 29
data.frame(
#### Diploid_05_1 TECAN-T15
## MTX data: 20190301_TECAN; non-MTX data: 20190314_TECAN
## Bad wells:T15 non-MTX: T15
data_file_name_01 = "20190301_TECAN/2019-02-27_17-09_D_05_01_t15_96_T15.txt",
data_file_name_02 =  "20190314_TECAN/2019-03-14_16-08_D_05_1_no_MTX_96_T15.txt",
map_name = "Diploid_05_all_PPI_1.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(c(1,13,25,37,49, 61,73,85))
                       ), #28 from 30
data.frame(
#### Diploid_05_2 TECAN-T2
## MTX data: 20190301_TECAN; non-MTX data: 20190314_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190301_TECAN/2019-02-27_17-10_D_05_02_t2_96_T2.txt",
data_file_name_02 =  "20190314_TECAN/2019-03-14_16-11_D_05_2_no_MTX_96_T2.txt",
map_name = "Diploid_05_all_PPI_2.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #26 from 30
data.frame(
#### Diploid_05_3 TECAN-T5
## MTX data: 20190301_TECAN; non-MTX data: 20190316_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190301_TECAN/2019-02-27_17-14_D_05_03_t5_96_T5.txt",
data_file_name_02 =  "20190316_TECAN/2019-03-16_09-49_D_05_3_no_MTX_96_T5.txt",
map_name = "Diploid_05_all_PPI_3.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #29 from 30
data.frame(
#### Diploid_05_4 TECAN-T7
## MTX data: 20190303_TECAN; non-MTX data: 20190313_TECAN
## Bad wells:0 non-MTX: 0
data_file_name_01 = "20190303_TECAN/2019-03-01_17-04_D_05_4_T7_96_T7.txt",
data_file_name_02 =  "20190313_TECAN/2019-03-13_09-06_D_5_4_No_MTX_96_T7.txt",
map_name = "Diploid_05_all_PPI_4.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(0)
                       ), #29 from 30
data.frame(
#### Diploid_05_5 TECAN-T14
## MTX data: 20190303_TECAN; non-MTX data: 20190314_TECAN
## Bad wells:91:93 (one PPI) non-MTX: 91:93(one PPI)
data_file_name_01 = "20190303_TECAN/2019-03-01_17-05_D_05_5_t14_96_T14.txt",
data_file_name_02 =  "20190314_TECAN/2019-03-14_16-07_D_05_5_no_MTX_96_T14.txt",
map_name = "Diploid_05_all_PPI_5.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(91:93)
                       ), #26 from 29
data.frame(
#### Diploid_05_6 TECAN-T15
## MTX data: 20190303_TECAN; non-MTX data: 20190316_TECAN
## Bad wells:T15 non-MTX: T15
data_file_name_01 = "20190303_TECAN/2019-03-01_17-06_D_05_6_t15_96_T15.txt",
data_file_name_02 =  "20190316_TECAN/2019-03-16_09-46_D_05_6_no_MTX_96_T15.txt",
map_name = "Diploid_05_all_PPI_6.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(c(1,13,25,37,49, 61,73,85))
                       ), #30 from 30
data.frame(
#### Diploid_05_7 TECAN-T2
## MTX data: 20190303_TECAN; non-MTX data: 20190316_TECAN
## Bad wells:64:93 (empty) non-MTX: 64:93(empty)
data_file_name_01 = "20190303_TECAN/2019-03-01_17-07_D_05_7_t2_96_T2.txt",
data_file_name_02 =  "20190316_TECAN/2019-03-16_09-45_D_05_7_no_MTX_96_T2.txt",
map_name = "Diploid_05_all_PPI_7.csv", 
time_window_01 = 40,
time_window_02 = 25,
bad_wells = c(64:93)
                       ) #14 from 20
)


z <- bind_rows(list_to_proc) %>% nest(bad_wells,.key="bad_wells") %>%
    mutate(bad_wells=map(bad_wells,function(x){list(x$bad_wells)})) %>% 
    mutate(datar=pmap(
            list(
                str_c("data/tecan_validation_assays/tecan_data/",
                    data_file_name_01), 
                str_c("data/tecan_validation_assays/tecan_data/",
                    data_file_name_02), 
                str_c("data/tecan_validation_assays/plate_strain_map/",
                    map_name), 
                time_window_01, time_window_02, 
                bad_wells),
            TECAN_data_AUC_correct)
        ) %>% unnest(datar) %>%
    select(PPI,`p-value`,one_tail_Q,statistic,ends_with("MTX"),ends_with("substract"))

write_csv(z,"tmp/tecan_validation_statistics.csv")


