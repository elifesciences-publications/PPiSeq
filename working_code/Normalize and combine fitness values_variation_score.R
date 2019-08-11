########### This script contain 5 steps of analysis: 
##########  (1) Normalize fitness of each lineage by DHFR(+) and DHFR(-) control strains (mean)
##########  (2) Put normized fitness values of the same PPI onto the same row
##########  (3) Combine normalized mean fitness values of different environments onto the same row
#########   (4) Put the fitness of PPI of two orientations across different environments onto the same row
#########   (5) Take the mean fitness (remove NA) for the same PPI with two orientations and calculate the variation score
source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

##### (1) Normalize fitenss values for each barcode by range of positive and negative fitness values
Fitness_normalization = function(DMSO_lineage){
        DMSO_DHFR_Pos= DMSO_lineage[which(DMSO_lineage[,1] == "positive_DHFR"),]
        DMSO_DHFR_Neg= DMSO_lineage[which(DMSO_lineage[,1] == "negative_non_DHFR"),]
        
        DMSO_DHFR_Pos_mean = mean(as.numeric(DMSO_DHFR_Pos[,4])) 
        DMSO_DHFR_Neg_mean = mean(as.numeric(DMSO_DHFR_Neg[,4])) 
        
        DMSO_lineage_chosen = DMSO_lineage[, c(1,2,4)]
        DMSO_lineage_chosen[,3] = (as.numeric(DMSO_lineage_chosen[,3]) - DMSO_DHFR_Neg_mean)/(DMSO_DHFR_Pos_mean - DMSO_DHFR_Neg_mean)
        return(DMSO_lineage_chosen)
}

##### (2) Put fitness values of the same PPI onto the same line
cluster_same_PPI = function(DMSO_lineage_norm_file, DMSO_multiple, output_file){
        PPI_multiple_RRS = DMSO_multiple[grep("Neg_PPI", DMSO_multiple[,1]),] #97
        PPI_multiple_PRS = DMSO_multiple[grep("Pos_PPI", DMSO_multiple[,1]),] #108
        PPI_multiple_pos = DMSO_multiple[grep("positive_DHFR", DMSO_multiple[,1]),] # 1
        PPI_multiple_neg = DMSO_multiple[grep("negative_non_DHFR", DMSO_multiple[,1]),] # 1
        PPI_multiple_control = rbind(PPI_multiple_RRS, PPI_multiple_PRS, PPI_multiple_pos, PPI_multiple_neg)
        PPI_multiple_select = DMSO_multiple[which(!DMSO_multiple[,1] %in% PPI_multiple_control[,1]),] 
        DMSO_lineage_norm= DMSO_lineage_norm_file[which(DMSO_lineage_norm_file[,1] %in% PPI_multiple_select[,1]), ] 
        
        PPI_unique= unique(PPI_multiple_select[,1])
        PPI_indiv_matrix= matrix(0, length(PPI_unique), 7)
        PPI_indiv_matrix[,1]= PPI_unique
        PPI_indiv_matrix[,2] = as.numeric(PPI_multiple_select[match(as.character(PPI_indiv_matrix[,1]), 
                                                                    as.character(PPI_multiple_select[,1])),2])
        index = 0
        for (i in 1:length(PPI_unique)){
                if(as.numeric(PPI_indiv_matrix[i,2]) == 4){
                        PPI_indiv_matrix[i, 4]= DMSO_lineage_norm[index + 4 -3, 3]
                        PPI_indiv_matrix[i, 5]= DMSO_lineage_norm[index + 4 -2, 3]
                        PPI_indiv_matrix[i, 6]= DMSO_lineage_norm[index + 4 -1, 3]
                        PPI_indiv_matrix[i, 7]= DMSO_lineage_norm[index + 4, 3] 
                        PPI_indiv_matrix[i,3] = mean(as.numeric(PPI_indiv_matrix[i,4:7]))
                        index = index + 4
                }else if(as.numeric(PPI_indiv_matrix[i,2] == 3)){
                        PPI_indiv_matrix[i, 4]= DMSO_lineage_norm[index + 3 -2, 3]
                        PPI_indiv_matrix[i, 5]= DMSO_lineage_norm[index + 3 -1, 3]
                        PPI_indiv_matrix[i, 6]= DMSO_lineage_norm[index + 3, 3]
                        PPI_indiv_matrix[i,3] = mean(as.numeric(PPI_indiv_matrix[i,4:6]))
                        index = index + 3
                }else if(as.numeric(PPI_indiv_matrix[i,2] == 2)){
                        PPI_indiv_matrix[i, 4]= DMSO_lineage_norm[index + 2 -1, 3]
                        PPI_indiv_matrix[i, 5]= DMSO_lineage_norm[index + 2, 3]
                        PPI_indiv_matrix[i,3] = mean(as.numeric(PPI_indiv_matrix[i,4:5]))
                        index = index + 2
                }
        }
        colnames(PPI_indiv_matrix)= c("PPI", "Barcodes", "Mean_fitness", "fit01", "fit02", "fit03", "fit04")
        csvWriter(PPI_indiv_matrix, output_file)
        #return(PPI_indiv_matrix)
}

### (1) and (2) normalize each PPI and put the normalized fitness onto the same row
setwd("/Volumes/zmliu_02/PPiseq/Combine_environments/")
SD_lineage = "lineage_fitness_files/SD_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/SD_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/SD_normalized_multiple.csv"
Normalize_environment = function(SD_lineage, SD_multiple, SD_multiple_normal){
        SD_fit = csvReader_T(SD_lineage)
        SD_fit_norm = Fitness_normalization(SD_fit)
        SD_multiple = csvReader_T(SD_multiple)
        cluster_same_PPI(SD_fit_norm, SD_multiple, SD_multiple_normal)
}

# SD
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

# SD2
SD_lineage = "lineage_fitness_files/SD2_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/SD2_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/SD2_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

# Raffinose
SD_lineage = "lineage_fitness_files/Raffinose_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/Raffinose_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/Raffinose_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

# NaCl
SD_lineage = "lineage_fitness_files/NaCl_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/NaCl_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/NaCl_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

#HU
SD_lineage = "lineage_fitness_files/Hydroxyurea_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/HU_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/Hydroxyurea_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

#H2O2
SD_lineage = "lineage_fitness_files/H2O2_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/H2O2_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/H2O2_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

#Forskolin
SD_lineage = "lineage_fitness_files/Forskolin_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/Forskolin_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/Forskolin_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

#FK506
SD_lineage = "lineage_fitness_files/FK506_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/FK506_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/FK506_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

# Doxorubicin
SD_lineage = "lineage_fitness_files/Doxorubicin_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/Dox_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/Doxorubicin_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

# Cold_16C
SD_lineage = "lineage_fitness_files/Cold_16C_PPI_barcodes_fitness_counts.csv"
SD_multiple = "PPI_multiple_files/16C_PPI_multiple_p.values.csv"
SD_multiple_normal = "Normalized_multiple_files/Cold_16C_normalized_multiple.csv"
Normalize_environment(SD_lineage, SD_multiple, SD_multiple_normal)

##### (3) Combine different environments and put the same PPI across different environments on the same row
setwd("/Volumes/zmliu_02/PPiseq/Combine_environments/Normalized_multiple_files/")

# order SD,SD2, Forskolin, FK506, Raffinose, NaCl, H2O2, Dox, 16C, HU
DMSO_norm = csvReader_T("SD_normalized_multiple.csv")
DMSO2_norm = csvReader_T('SD2_normalized_multiple.csv')
Forskolin_norm = csvReader_T("Forskolin_normalized_multiple.csv")
FK506_norm = csvReader_T("FK506_normalized_multiple.csv")
Raffinose_norm = csvReader_T("Raffinose_normalized_multiple.csv")
NaCl_norm = csvReader_T("NaCl_normalized_multiple.csv")
H2O2_norm = csvReader_T("H2O2_normalized_multiple.csv")
Dox_norm = csvReader_T("Doxorubicin_normalized_multiple.csv")
cold_norm = csvReader_T("Cold_16C_normalized_multiple.csv")
HU_norm = csvReader_T("Hydroxyurea_normalized_multiple.csv")

#### Combine all the PPIs that has been barcoded twice in one of environments
all = unique(c(DMSO_norm[,1], DMSO2_norm[,1], Forskolin_norm[,1], FK506_norm[,1], Raffinose_norm[,1], 
                   NaCl_norm[,1], H2O2_norm[,1], Dox_norm[,1], cold_norm[,1], HU_norm[,1])) #16431

matrix_all = matrix(0, length(all), 11)
matrix_all[,1] = all
matrix_all[,2] = as.numeric(DMSO_norm[match(all, DMSO_norm[,1]),3])
matrix_all[,3] = as.numeric(DMSO2_norm[match(all, DMSO2_norm[,1]),3])
matrix_all[,4] = as.numeric(H2O2_norm[match(all, H2O2_norm[,1]),3])
matrix_all[,5] = as.numeric(HU_norm[match(all, HU_norm[,1]),3])
matrix_all[,6] = as.numeric(Dox_norm[match(all, Dox_norm[,1]),3])
matrix_all[,7] = as.numeric(Forskolin_norm[match(all, Forskolin_norm[,1]),3])
matrix_all[,8] = as.numeric(Raffinose_norm[match(all, Raffinose_norm[,1]),3])
matrix_all[,9] = as.numeric(NaCl_norm[match(all, NaCl_norm[,1]),3])
matrix_all[,10] = as.numeric(cold_norm[match(all, cold_norm[,1]),3])
matrix_all[,11] = as.numeric(FK506_norm[match(all, FK506_norm[,1]),3])

colnames(matrix_all) = c("PPI", "SD", "SD2","H2O2", "HU", "Dox", "Forskolin", 
                         "Raffinose", "NaCl", "16C", "FK506")
csvWriter(matrix_all, "All_PPI_environments_normalized_fit.csv")

######## Only consider positive PPIs in one of environments
DMSO_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/SD_Pos_PPI_real.csv")
DMSO2_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/SD2_Pos_PPI_real.csv")
Forskolin_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/Forskolin_Pos_PPI_real.csv")
FK506_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/FK506_Pos_PPI_real.csv")
Raffinose_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/Raffinose_Pos_PPI_real.csv")
NaCl_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/NaCl_Pos_PPI_real.csv")
H2O2_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/H2O2_Pos_PPI_real.csv")
Dox_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/Doxorubicin_Pos_PPI_real.csv")
cold_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/Cold_16C_Pos_PPI_real.csv")
HU_real = csvReader_T("/Volumes/zmliu_02/PPiseq/Combine_environments/Positive_PPI_remove_promiscuous/Hydroxyurea_Pos_PPI_real.csv")


DMSO_norm_pos = DMSO_norm[which(DMSO_norm[,1] %in% DMSO_real[,1]),c(1,3)]
DMSO2_norm_pos = DMSO2_norm[which(DMSO2_norm[,1] %in% DMSO2_real[,1]),c(1,3)]
Forskolin_norm_pos = Forskolin_norm[which(Forskolin_norm[,1] %in% Forskolin_real[,1]),c(1,3)]
FK506_norm_pos = FK506_norm[which(FK506_norm[,1] %in% FK506_real[,1]),c(1,3)]
Raffinose_norm_pos = Raffinose_norm[which(Raffinose_norm[,1] %in% Raffinose_real[,1]),c(1,3)]
NaCl_norm_pos = NaCl_norm[which(NaCl_norm[,1] %in% NaCl_real[,1]),c(1,3)]
H2O2_norm_pos = H2O2_norm[which(H2O2_norm[,1] %in% H2O2_real[,1]),c(1,3)]
Dox_norm_pos = Dox_norm[which(Dox_norm[,1] %in% Dox_real[,1]),c(1,3)]
cold_norm_pos = cold_norm[which(cold_norm[,1] %in% cold_real[,1]),c(1,3)]
HU_norm_pos = HU_norm[which(HU_norm[,1] %in% HU_real[,1]),c(1,3)]

all_pos = unique(c(DMSO_norm_pos[,1], DMSO2_norm_pos[,1], Forskolin_norm_pos[,1], FK506_norm_pos[,1], Raffinose_norm_pos[,1], 
                   NaCl_norm_pos[,1], H2O2_norm_pos[,1], Dox_norm_pos[,1], cold_norm_pos[,1], HU_norm_pos[,1])) #16431

matrix_pos = matrix(0, length(all_pos), 11)
matrix_pos[,1] = all_pos
matrix_pos[,2] = as.numeric(DMSO_norm_pos[match(all_pos, DMSO_norm_pos[,1]),2])
matrix_pos[,3] = as.numeric(DMSO2_norm_pos[match(all_pos, DMSO2_norm_pos[,1]),2])
matrix_pos[,4] = as.numeric(H2O2_norm_pos[match(all_pos, H2O2_norm_pos[,1]),2])
matrix_pos[,5] = as.numeric(HU_norm_pos[match(all_pos, HU_norm_pos[,1]),2])
matrix_pos[,6] = as.numeric(Dox_norm_pos[match(all_pos, Dox_norm_pos[,1]),2])
matrix_pos[,7] = as.numeric(Forskolin_norm_pos[match(all_pos, Forskolin_norm_pos[,1]),2])
matrix_pos[,8] = as.numeric(Raffinose_norm_pos[match(all_pos, Raffinose_norm_pos[,1]),2])
matrix_pos[,9] = as.numeric(NaCl_norm_pos[match(all_pos, NaCl_norm_pos[,1]),2])
matrix_pos[,10] = as.numeric(cold_norm_pos[match(all_pos, cold_norm_pos[,1]),2])
matrix_pos[,11] = as.numeric(FK506_norm_pos[match(all_pos, FK506_norm_pos[,1]),2])

colnames(matrix_pos) = c("PPI", "SD", "SD2","H2O2", "HU", "Dox", "Forskolin", 
                         "Raffinose", "NaCl", "16C", "FK506")
csvWriter(matrix_pos, "Pos_PPI_normalized_fit.csv")

##############################################################
# Keep the primary normalized fitness values for PPIs even not detected in an evironment

matrix_pos = matrix(0, length(all_pos), 11)
matrix_pos[,1] = all_pos
matrix_pos[,2] = as.numeric(DMSO_norm[match(all_pos, DMSO_norm[,1]),3])
matrix_pos[,3] = as.numeric(DMSO2_norm[match(all_pos, DMSO2_norm[,1]),3])
matrix_pos[,4] = as.numeric(H2O2_norm[match(all_pos, H2O2_norm[,1]),3])
matrix_pos[,5] = as.numeric(HU_norm[match(all_pos, HU_norm[,1]),3])
matrix_pos[,6] = as.numeric(Dox_norm[match(all_pos, Dox_norm[,1]),3])
matrix_pos[,7] = as.numeric(Forskolin_norm[match(all_pos, Forskolin_norm[,1]),3])
matrix_pos[,8] = as.numeric(Raffinose_norm[match(all_pos, Raffinose_norm[,1]),3])
matrix_pos[,9] = as.numeric(NaCl_norm[match(all_pos, NaCl_norm[,1]),3])
matrix_pos[,10] = as.numeric(cold_norm[match(all_pos, cold_norm[,1]),3])
matrix_pos[,11] = as.numeric(FK506_norm[match(all_pos, FK506_norm[,1]),3])

colnames(matrix_pos) = c("PPI", "SD", "SD2","H2O2", "HU", "Dox", "Forskolin", 
                         "Raffinose", "NaCl", "16C", "FK506")
csvWriter(matrix_pos, "Pos_PPI_normalized_fit_primary.csv")

### (4) Put the fitness of PPI of two orientations across different environments onto the same row

setwd("/Volumes/zmliu_02/PPiseq/Combine_environments/Normalized_multiple_files/")
PPI_fit_norm = csvReader_T("Pos_PPI_normalized_fit_primary.csv") # 16431
## Make the PPI name for the count and normalized fitness consistent
PPI_count = csvReader_T("PPI_environment_count_summary.csv") 
PPI_dup = mark_duplicates_fast(PPI_fit_norm[,1]) # 15656
PPI_dup_01 = c("0", "0")
for(i in 1:nrow(PPI_count)){
        if(PPI_count[i,1] %in% PPI_dup[,1]){
                index = which(PPI_dup[,1] == PPI_count[i,1])
                PPI_dup_01 = rbind(PPI_dup_01, PPI_dup[index,])
        }else if (PPI_count[i,1] %in% PPI_dup[,2]){
                index = which(PPI_dup[,2] == PPI_count[i,1])
                PPI_dup_01 = rbind(PPI_dup_01, PPI_dup[index, c(2,1)])
        }
}

PPI_dup = PPI_dup_01[2:nrow(PPI_dup_01),]
matrix = matrix(NA, nrow(PPI_dup), 2* ncol(PPI_fit_norm))
for(i in 1:nrow(PPI_dup)){
        if (PPI_dup[i,2] != "0"){
                matrix[i,1:11] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
                matrix[i,12:22] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,2]),]
        }else{
                matrix[i,1:11] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
        }   
        
}
colnames(matrix) = c("PPI", "SD", "SD2", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", 
                     "NaCl", "16C", "FK506", "PPI_opposite", "SD", "SD2" ,"H2O2", "HU",
                     "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
mean_fitness = rep(0, nrow(matrix))
for(i in 1:length(mean_fitness)){
        mean_fitness[i] = mean(na.omit(as.numeric(matrix[i, c(2:11, 13:22)])))
}

matrix_final = cbind(matrix[,1], mean_fitness, matrix[,2:ncol(matrix)])
csvWriter(matrix_final, "Normalized_fitness_PPI_all_primary.csv")

## only consider positive PPIs, all negative PPIs have fitness of 0
PPI_fit_norm = csvReader_T("Pos_PPI_normalized_fit.csv") # 16431
PPI_count = csvReader_T("PPI_environment_count_summary.csv") 
PPI_dup = mark_duplicates_fast(PPI_fit_norm[,1]) # 15656
PPI_dup_01 = c("0", "0")
for(i in 1:nrow(PPI_count)){
        if(PPI_count[i,1] %in% PPI_dup[,1]){
                index = which(PPI_dup[,1] == PPI_count[i,1])
                PPI_dup_01 = rbind(PPI_dup_01, PPI_dup[index,])
        }else if (PPI_count[i,1] %in% PPI_dup[,2]){
                index = which(PPI_dup[,2] == PPI_count[i,1])
                PPI_dup_01 = rbind(PPI_dup_01, PPI_dup[index, c(2,1)])
        }
}

PPI_dup = PPI_dup_01[2:nrow(PPI_dup_01),]
matrix = matrix(NA, nrow(PPI_dup), 2* ncol(PPI_fit_norm))
for(i in 1:nrow(PPI_dup)){
        if (PPI_dup[i,2] != "0"){
                matrix[i,1:11] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
                matrix[i,12:22] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,2]),]
        }else{
                matrix[i,1:11] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
        }   
        
}
colnames(matrix) = c("PPI", "SD", "SD2", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", 
                     "NaCl", "16C", "FK506", "PPI_opposite", "SD", "SD2" ,"H2O2", "HU",
                     "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
mean_fitness = rep(0, nrow(matrix))
for(i in 1:length(mean_fitness)){
        mean_fitness[i] = mean(na.omit(as.numeric(matrix[i, c(2:11, 13:22)])))
}

matrix_final = cbind(matrix[,1], mean_fitness, matrix[,2:ncol(matrix)])
csvWriter(matrix_final, "Normalized_fitness_PPI_all_pos.csv")

### Consider normalzied negative fitness values to be 0
PPI_fit_norm = csvReader_T("Pos_PPI_normalized_fit_primary.csv") # 16431
for (i in 1:nrow(PPI_fit_norm)){
        for (j in 2:ncol(PPI_fit_norm)){
                if ((!is.na(PPI_fit_norm[i,j])) & as.numeric(PPI_fit_norm[i,j]) < 0){
                        PPI_fit_norm[i,j] = 0
                }
        }
}
PPI_count = csvReader_T("PPI_environment_count_summary.csv") 
PPI_dup = mark_duplicates_fast(PPI_fit_norm[,1]) # 15656
PPI_dup_01 = c("0", "0")
for(i in 1:nrow(PPI_count)){
        if(PPI_count[i,1] %in% PPI_dup[,1]){
                index = which(PPI_dup[,1] == PPI_count[i,1])
                PPI_dup_01 = rbind(PPI_dup_01, PPI_dup[index,])
        }else if (PPI_count[i,1] %in% PPI_dup[,2]){
                index = which(PPI_dup[,2] == PPI_count[i,1])
                PPI_dup_01 = rbind(PPI_dup_01, PPI_dup[index, c(2,1)])
        }
}

PPI_dup = PPI_dup_01[2:nrow(PPI_dup_01),]
matrix = matrix(NA, nrow(PPI_dup), 2* ncol(PPI_fit_norm))
for(i in 1:nrow(PPI_dup)){
        if (PPI_dup[i,2] != "0"){
                matrix[i,1:11] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
                matrix[i,12:22] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,2]),]
        }else{
                matrix[i,1:11] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
        }   
        
}
colnames(matrix) = c("PPI", "SD", "SD2", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", 
                     "NaCl", "16C", "FK506", "PPI_opposite", "SD", "SD2" ,"H2O2", "HU",
                     "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
mean_fitness = rep(0, nrow(matrix))
for(i in 1:length(mean_fitness)){
        mean_fitness[i] = mean(na.omit(as.numeric(matrix[i, c(2:11, 13:22)])))
}

matrix_final = cbind(matrix[,1], mean_fitness, matrix[,2:ncol(matrix)])
csvWriter(matrix_final, "Normalized_fitness_PPI_all_neg_zero.csv")


########### (5)Caculate an variation score for each PPI based on their fitness values across different environments
#### Only consider positive PPIs all negative PPIs have fitness of 0
PPI_norm = csvReader_T("Normalized_fitness_PPI_all_pos.csv")
PPI_norm_matrix = matrix(0, nrow(PPI_norm), 11)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
        for(j in 3:12){
                mean_fit = mean(as.numeric(na.omit(PPI_norm[i,c(j, j + 11)])))
                if(is.na(mean_fit)){
                        PPI_norm_matrix[i, j-1] = 0
                }
                else{
                        PPI_norm_matrix[i,j-1] = mean_fit
                }
                
        }  
}
variation_score = rep(0, nrow(PPI_norm_matrix))
environment_number = rep(0, nrow(PPI_norm_matrix))
for(i in 1:length(variation_score)){
        fitness = as.numeric(PPI_norm_matrix[i,2:11])
        variation_score[i] = sd(fitness)/mean(fitness)
        environment_number[i] = length(which(fitness != 0))
}

# Include the variation score into the matrix, and output the matrix
PPI_norm_matrix_final = cbind(PPI_norm_matrix[,1],environment_number, variation_score,
                              PPI_norm_matrix[,2:11])
colnames(PPI_norm_matrix_final) = c("PPI", "Environment_number", "Variation_score",
                                    "SD","SD2", "H2O2", "HU", "Dox", "Forskolin", "Raffinose",
                                    "NaCl", "16C", "FK506")
csvWriter(PPI_norm_matrix_final, "Variation_score_PPI_environment_pos.csv")
cor(environment_number, variation_score, method="spearman") # -0.9057337


###Use the normalized value and Keep the primary normalized fitess even a PPI is called negative
PPI_norm = csvReader_T("Normalized_fitness_PPI_all_primary.csv")
PPI_count = csvReader_T("PPI_environment_count_summary.csv")
PPI_norm_matrix = matrix(0, nrow(PPI_norm), 11)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
        for(j in 3:12){
                mean_fit = mean(as.numeric(na.omit(PPI_norm[i,c(j, j + 11)])))
                if(is.na(mean_fit)){
                        PPI_norm_matrix[i, j-1] = 0
                }
                else{
                        PPI_norm_matrix[i,j-1] = mean_fit
                }
                
        }  
}
variation_score = rep(0, nrow(PPI_norm_matrix))
for(i in 1:length(variation_score)){
        fitness = as.numeric(PPI_norm_matrix[i,2:11])
        variation_score[i] = sd(fitness)/mean(fitness)
}
environment_number = as.numeric(PPI_count[match(PPI_norm_matrix[,1], PPI_count[,1]),2])
# Include the variation score into the matrix, and output the matrix
PPI_norm_matrix_final = cbind(PPI_norm_matrix[,1],environment_number, variation_score,
                              PPI_norm_matrix[,2:11])
colnames(PPI_norm_matrix_final) = c("PPI", "Environment_number", "Variation_score",
                                    "SD","SD2", "H2O2", "HU", "Dox", "Forskolin", "Raffinose",
                                    "NaCl", "16C", "FK506")
csvWriter(PPI_norm_matrix_final, "Variation_score_PPI_environment_primary.csv")
cor(environment_number, variation_score, method="spearman") # -0.6790267


########### Consider negative fitness values to be zero
PPI_norm = csvReader_T("Normalized_fitness_PPI_all_neg_zero.csv")
PPI_count = csvReader_T("PPI_environment_count_summary.csv")
PPI_norm_matrix = matrix(0, nrow(PPI_norm), 11)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
        for(j in 3:12){
                mean_fit = mean(as.numeric(na.omit(PPI_norm[i,c(j, j + 11)])))
                if(is.na(mean_fit)){
                        PPI_norm_matrix[i, j-1] = 0
                }
                else{
                        PPI_norm_matrix[i,j-1] = mean_fit
                }
                
        }  
}
variation_score = rep(0, nrow(PPI_norm_matrix))
for(i in 1:length(variation_score)){
        fitness = as.numeric(PPI_norm_matrix[i,2:11])
        variation_score[i] = sd(fitness)/mean(fitness)
}
environment_number = as.numeric(PPI_count[match(PPI_norm_matrix[,1], PPI_count[,1]),2])
# Include the variation score into the matrix, and output the matrix
PPI_norm_matrix_final = cbind(PPI_norm_matrix[,1],environment_number, variation_score,
                              PPI_norm_matrix[,2:11])
colnames(PPI_norm_matrix_final) = c("PPI", "Environment_number", "Variation_score",
                                    "SD","SD2", "H2O2", "HU", "Dox", "Forskolin", "Raffinose",
                                    "NaCl", "16C", "FK506")
csvWriter(PPI_norm_matrix_final, "Variation_score_PPI_environment_neg_zero.csv")
cor(environment_number, variation_score, method="spearman") # -0.7038279

