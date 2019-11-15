### This script contains several steps of analysis

source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments/lineage_fitness_files/")
DMSO = dataFrameReader_F("SD_known_PPI_fitness_barcodes_sorted.csv")
colnames(DMSO) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error", 
                   "Counts_G0", "Counts_G3", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(DMSO, "SD_PPI_barcodes_fitness_counts.csv")

DMSO_2 = dataFrameReader_F("SD2_known_PPI_fitness_barcodes_sorted.csv")
colnames(DMSO_2) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                   "Counts_G0", "Counts_G3", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(DMSO_2, "SD2_PPI_barcodes_fitness_counts.csv")

H2O2 = dataFrameReader_F("H2O2_known_PPI_fitness_barcodes_sorted.csv")
colnames(H2O2) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                   "Counts_G0", "Counts_G3", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(H2O2, "H2O2_PPI_barcodes_fitness_counts.csv")

HU = dataFrameReader_F("Hydroxyurea_known_PPI_fitness_barcodes_sorted.csv")
colnames(HU) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                   "Counts_G0", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(HU, "Hydroxyurea_PPI_barcodes_fitness_counts.csv")

Dox = dataFrameReader_F("Doxorubicin_known_PPI_fitness_barcodes_sorted.csv")
colnames(Dox) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                   "Counts_G0", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(Dox, "Doxorubicin_PPI_barcodes_fitness_counts.csv")

Forskolin = dataFrameReader_F("Forskolin_known_PPI_fitness_barcodes_sorted.csv")
colnames(Forskolin) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                  "Counts_G0", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(Forskolin, "Forskolin_PPI_barcodes_fitness_counts.csv")

Raffinose = dataFrameReader_F("Raffinose_known_PPI_fitness_barcodes_sorted.csv")
colnames(Raffinose) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                        "Counts_G0", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(Raffinose, "Raffinose_PPI_barcodes_fitness_counts.csv")

FK506 = dataFrameReader_F("FK506_known_PPI_fitness_barcodes_sorted.csv")
colnames(FK506) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                        "Counts_G0", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(FK506, "FK506_PPI_barcodes_fitness_counts.csv")

NaCl = dataFrameReader_F("NaCl_known_PPI_fitness_barcodes_sorted.csv")
colnames(NaCl) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                    "Counts_G0", "Counts_G6", "Counts_G12", "Counts_G18")
csvWriter(NaCl, "NaCl_PPI_barcodes_fitness_counts.csv")

cold = dataFrameReader_F("Cold_16C_known_PPI_fitness_barcodes_sorted.csv")
colnames(cold) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error(d)", 
                   "Counts_G0", "Counts_G6", "Counts_G12", "Counts_G18")
csvWriter(cold, "cold_16C_PPI_barcodes_fitness_counts.csv")

##### After finishing adding headers into these files, delete the primary ones

####(2) Summarize the promiscuous proteins and remove all PPIs that contain promiscuous proteins reported in two environments
#### Check the number of environments in which a promiscuous protein is detected
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Promiscuous_PPIs/")
DMSO = csvReader_T("SD_Promiscuous_proteins.csv")
DMSO_2 = csvReader_T("SD2_Promiscuous_proteins.csv")
H2O2 = csvReader_T("H2O2_Promiscuous_proteins.csv")
HU = csvReader_T("Hydroxyurea_Promiscuous_proteins.csv")
Dox = csvReader_T("Doxorubicin_Promiscuous_proteins.csv")
Forskolin = csvReader_T("Forskolin_Promiscuous_proteins.csv")
Raffinose = csvReader_T("Raffinose_Promiscuous_proteins.csv")
NaCl = csvReader_T("NaCl_Promiscuous_proteins.csv")
FK506 = csvReader_T("FK506_Promiscuous_proteins.csv")
cold = csvReader_T("Cold_16C_Promiscuous_proteins.csv")
extract_promiscuous_protein = function(DMSO){
        temp = unique(c(split_string_vector(DMSO[,1])[,1:2]))
        temp_HO = temp[grep("HO", temp)]
        temp_unique = temp[which(!temp %in% temp_HO)]
        return(temp_unique)
}
DMSO_bad = extract_promiscuous_protein(DMSO)
DMSO_2_bad = extract_promiscuous_protein(DMSO_2)
H2O2_bad = extract_promiscuous_protein(H2O2)
HU_bad = extract_promiscuous_protein(HU)
Dox_bad = extract_promiscuous_protein(Dox)
Forskolin_bad = extract_promiscuous_protein(Forskolin)
Raffinose_bad = extract_promiscuous_protein(Raffinose)
NaCl_bad = extract_promiscuous_protein(NaCl)
FK506_bad = extract_promiscuous_protein(FK506)
cold_bad = extract_promiscuous_protein(cold)

final_bad = unique(c(DMSO_bad, DMSO_2_bad, H2O2_bad, HU_bad, Dox_bad, Forskolin_bad,
                     Raffinose_bad, NaCl_bad, FK506_bad, cold_bad)) # 59

PPI_list = list(DMSO_bad,DMSO_2_bad, H2O2_bad, HU_bad, Dox_bad, Forskolin_bad,
                Raffinose_bad, NaCl_bad, FK506_bad, cold_bad)
matrix_bad = matrix(0, length(final_bad), 12)
colnames(matrix_bad) = c("PPI","Positive_sum","SD","SD2","H2O2", "HU", "Dox", "Forskolin",
                         "Raffinose", "NaCl", "FK506", "Cold")
matrix_bad[,1] = final_bad
for(i in 1:nrow(matrix_bad)){
        for(j in 1:length(PPI_list)){
                if(matrix_bad[i,1] %in% PPI_list[[j]]){
                        matrix_bad[i,j +2] = 1
                }
        }
        matrix_bad[i,2] = sum(as.numeric(matrix_bad[i,3: ncol(matrix_bad)]))
}

matrix_bad_order = matrix_bad[order(matrix_bad[,2], decreasing = T),]
csvWriter(matrix_bad_order, "Promiscuous_protein_summary.csv")

### Remove PPIs that contain a protein that is identified as promiscuous protein in >= 2 environments
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/")
DMSO_pos = csvReader_T("SD_Pos_PPI_real.csv") # 6071
DMSO_2_pos = csvReader_T("SD2_Pos_PPI_real.csv") # 5322
H2O2_pos = csvReader_T("H2O2_Pos_PPI_real.csv") # 6176
HU_pos = csvReader_T("Hydroxyurea_Pos_PPI_real.csv") # 7804
Dox_pos = csvReader_T("Doxorubicin_Pos_PPI_real.csv") # 2482
Forskolin_pos = csvReader_T("Forskolin_Pos_PPI_real.csv") # 4612
NaCl_pos = csvReader_T("NaCl_Pos_PPI_real.csv") # 1337
Raffinose_pos = csvReader_T("Raffinose_Pos_PPI_real.csv") # 4920
FK506_pos = csvReader_T("FK506_Pos_PPI_real.csv") #3766
Cold_pos = csvReader_T("Cold_16C_Pos_PPI_real.csv") # 2604
promiscuous_protein = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Promiscuous_PPIs/Promiscuous_protein_summary.csv")
promiscuous_protein_chosen = promiscuous_protein[which(as.numeric(promiscuous_protein[,2])>=2),1]
remove_PPI_protein = function(DMSO_pos, promiscuous_protein_chosen){
        check = rep(0, nrow(DMSO_pos))
        PPI = split_string_vector(DMSO_pos[,1])
        for(i in 1:nrow(DMSO_pos)){
                if (PPI[i,1] %in% promiscuous_protein_chosen | PPI[i,2] %in% promiscuous_protein_chosen){
                        check[i] = 1
                }
        }
        DMSO_pos_filter = DMSO_pos[which(check != 1),]
        return(DMSO_pos_filter)
}
DMSO_pos_filter = remove_PPI_protein(DMSO_pos, promiscuous_protein_chosen) # 6065 from 6071
DMSO_2_pos_filter = remove_PPI_protein(DMSO_2_pos, promiscuous_protein_chosen) # 5317 from 5322
H2O2_pos_filter = remove_PPI_protein(H2O2_pos, promiscuous_protein_chosen) # 6164 from 6176
HU_pos_filter = remove_PPI_protein(HU_pos, promiscuous_protein_chosen) # 6484 from 7804
Dox_pos_filter = remove_PPI_protein(Dox_pos, promiscuous_protein_chosen) # 2307 from 2482
Forskolin_pos_filter = remove_PPI_protein(Forskolin_pos, promiscuous_protein_chosen) # 4424 from 4612
NaCl_pos_filter = remove_PPI_protein(NaCl_pos, promiscuous_protein_chosen) # 1274 from 1337
Raffinose_pos_filter = remove_PPI_protein(Raffinose_pos, promiscuous_protein_chosen) # 4915 from 4920
FK506_pos_filter = remove_PPI_protein(FK506_pos, promiscuous_protein_chosen) # 3766 from 3766
cold_pos_filter = remove_PPI_protein(Cold_pos, promiscuous_protein_chosen) # 2577 from 2604

csvWriter(DMSO_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/SD_Pos_PPI_real.csv")
csvWriter(DMSO_2_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/SD2_Pos_PPI_real.csv")
csvWriter(H2O2_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/H2O2_Pos_PPI_real.csv")
csvWriter(HU_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Hydroxyurea_Pos_PPI_real.csv")
csvWriter(Dox_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Doxorubicin_Pos_PPI_real.csv")
csvWriter(Forskolin_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Forskolin_Pos_PPI_real.csv")
csvWriter(NaCl_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/NaCl_Pos_PPI_real.csv")
csvWriter(Raffinose_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Raffinose_Pos_PPI_real.csv")
csvWriter(FK506_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/FK506_Pos_PPI_real.csv")
csvWriter(cold_pos_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Cold_16C_Pos_PPI_real.csv")


#(3) Make a file that contain mean fitness, sd, p.value, Q.value, positive (1: positive), control(1:control)
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_multiple_files/")
Combine_matrix = function(DMSO, DMSO_pos_control, DMSO_pos_real, output){
        DMSO_pos_PPI = c(DMSO_pos_control[,1], DMSO_pos_real[,1])
        pos_index = rep(0, nrow(DMSO))
        pos_index[which(as.character(DMSO[,1]) %in% DMSO_pos_PPI)] = 1
        length(which(pos_index == 1))
        DMSO_label = cbind(DMSO, pos_index)
        colnames(DMSO_label) = c("PPI", "Number_of_barcodes", "Mean_fitness","SD","t", "P_value", "FDR_adjusted_P_value", "Positive")
        csvWriter(DMSO_label, output)
}

#### SD
DMSO = dataFrameReader_T("SD_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/SD_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/SD_Pos_PPI_real.csv")
output = "SD_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

#### SD_2
DMSO = dataFrameReader_T("SD2_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/SD2_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/SD2_Pos_PPI_real.csv")
output = "SD2_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

### H2O2
DMSO = dataFrameReader_T("H2O2_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/H2O2_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/H2O2_Pos_PPI_real.csv")
output = "H2O2_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

### HU
DMSO = dataFrameReader_T("Hydroxyurea_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/Hydroxyurea_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Hydroxyurea_Pos_PPI_real.csv")
output = "Hydroxyurea_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

### Dox
DMSO = dataFrameReader_T("Doxorubicin_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/Doxorubicin_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Doxorubicin_Pos_PPI_real.csv")
output = "Doxorubicin_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

### Forskolin
DMSO = dataFrameReader_T("Forskolin_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/Forskolin_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Forskolin_Pos_PPI_real.csv")
output = "Forskolin_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

### Raffinose
DMSO = dataFrameReader_T("Raffinose_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/Raffinose_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Raffinose_Pos_PPI_real.csv")
output = "Raffinose_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

### NaCl
DMSO = dataFrameReader_T("NaCl_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/NaCl_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/NaCl_Pos_PPI_real.csv")
output = "NaCl_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

### 16C
DMSO = dataFrameReader_T("Cold_16C_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/Cold_16C_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/Cold_16C_Pos_PPI_real.csv")
output = "Cold_16C_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

### FK506
DMSO = dataFrameReader_T("FK506_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/PPI_each_environment/FK506_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Positive_PPI_remove_promiscuous/FK506_Pos_PPI_real.csv")
output = "FK506_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)
