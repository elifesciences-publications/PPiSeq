apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
## This script describes all the dynamic threshold in each environment
setwd("/Volumes//zmliu_02/PPiseq_03/")
SD_merge = read.csv("DMSO_merge/reference_set/p_value_2/Threshold_metrics_final.csv") 
FK506 = read.csv("FK506/reference_set/p_value/Threshold_metrics_final.csv")
H2O2 = read.csv("H2O2/reference_set/p_value/Threshold_metrics_final.csv")
HU = read.csv("HU/reference_set/p_value/Threshold_metrics_final.csv")
NaCl = read.csv("NaCl_0.4M/reference_set/p_value/Threshold_metrics_final.csv") 
Forskolin = read.csv("Forskolin/reference_set/p_value/Threshold_metrics_final.csv")
Raffinose = read.csv("Raffinose/reference_set/p_value/Threshold_metrics_final.csv")
Dox = read.csv("Dox/reference_set/p_value/Threshold_metrics_final.csv")
cold = read.csv("16C/reference_set/p_value/Threshold_metrics_final.csv")

### Find out the dynamic thresholds which render the maximum MCC
SD_merge_MCC = SD_merge[which(SD_merge$MCC == max(SD_merge$MCC)),]
FK506_MCC = FK506[which(FK506$MCC == max(FK506$MCC)),]
H2O2_MCC = H2O2[which(H2O2$MCC == max(H2O2$MCC)),]
HU_MCC = HU[which(HU$MCC == max(HU$MCC)),]
NaCl_MCC = NaCl[which(NaCl$MCC == max(NaCl$MCC)),]
Forskolin_MCC = Forskolin[which(Forskolin$MCC == max(Forskolin$MCC)),]
Raffinose_MCC = Raffinose[which(Raffinose$MCC == max(Raffinose$MCC)),]
Dox_MCC = Dox[which(Dox$MCC == max(Dox$MCC)),]
cold_MCC = cold[which(cold$MCC == max(cold$MCC)),]

MCC_matrix = rbind(SD_merge_MCC, FK506_MCC, H2O2_MCC, HU_MCC, NaCl_MCC, 
                   Forskolin_MCC, Raffinose_MCC, Dox_MCC, cold_MCC)

Condition = c("SD", "FK506", "H2O2", "Hydroxyurea", "NaCl", "Forskolin", "Raffinose", "Doxorubicin", "Cold_16C")
MCC_matrix = cbind(Condition, MCC_matrix)

# A function that can find out the maximum FPR that is smaller than 0.001
max_FPR = function(vector, FPR_chosen){
  vector_chosen = vector[which(vector <= FPR_chosen)]
  return(max(vector_chosen))
}

SD_merge_chosen = SD_merge[which(SD_merge$FPR == max_FPR(SD_merge$FPR, 0.001)),] # 0.79
FK506_chosen = FK506[which(FK506$FPR == max_FPR(FK506$FPR, 0.001)),] # 0.78
H2O2_chosen = H2O2[which(H2O2$FPR == max_FPR(H2O2$FPR, 0.001)),] # 0.80
HU_chosen = HU[which(HU$FPR == max_FPR(HU$FPR, 0.001)),] # 0.80
NaCl_chosen = NaCl[which(NaCl$FPR == max_FPR(NaCl$FPR, 0.001)),] # 0.76
Forskolin_chosen = Forskolin[which(Forskolin$FPR == max_FPR(Forskolin$FPR, 0.001)),] # 0.77
Raffinose_chosen = Raffinose[which(Raffinose$Combined_PPV == 0.56),] # arbitrarily choose 0.56
Dox_chosen = Dox[which(Dox$FPR == max_FPR(Dox$FPR, 0.001)),] # 0.82
cold_chosen = cold[which(cold$FPR == max_FPR(cold$FPR, 0.001)),] # 0.51

FPR_matrix = rbind(SD_merge_chosen, FK506_chosen, H2O2_chosen, HU_chosen, NaCl_chosen, 
                   Forskolin_chosen, Raffinose_chosen, Dox_chosen, cold_chosen)
FPR_matrix = cbind(Condition, FPR_matrix)
write.csv(FPR_matrix, "/Volumes/zmliu_02/PPiseq_03/strict_threshold/FPR_threshold_environments.csv", row.names = F, quote = F)

# Add one more column for each threshold and then plot the difference
FPR_matrix$Threshold = "low_FPR"
MCC_matrix$Threshold = "balance"
# Merge two matrix
two_threshold = rbind(FPR_matrix, MCC_matrix)
write.csv(two_threshold, "/Volumes/zmliu_02/PPiseq_03/strict_threshold/MCC_FPR_threshold_environments.csv", row.names = F, quote = F)

## plot the comparison between two thresholds
library(ggplot2)
ggplot(data = two_threshold, aes(x= Condition, y = FPR, fill = Threshold, color = Threshold)) +
  geom_bar(stat="identity", position="dodge")+
  scale_color_manual(values = apple_colors[c(5,7)]) +
  scale_fill_manual(values = apple_colors[c(5,7)])+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("strict_threshold/FPR_two_threshold.pdf", width = 5, height =4)

### PPV
ggplot(data = two_threshold, aes(x= Condition, y = PPV, fill = Threshold, color = Threshold)) +
  geom_bar(stat="identity", position="dodge")+
  scale_color_manual(values = apple_colors[c(5,7)]) +
  scale_fill_manual(values = apple_colors[c(5,7)])+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("strict_threshold/PPV_two_threshold.pdf", width = 5, height =4)

### Number of PPIs
ggplot(data = two_threshold, aes(x= Condition, y = all, fill = Threshold, color = Threshold)) +
  geom_bar(stat="identity", position="dodge")+
  scale_color_manual(values = apple_colors[c(5,7)]) +
  scale_fill_manual(values = apple_colors[c(5,7)])+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("strict_threshold/Number_PPIs_two_threshold.pdf", width = 5, height =4)


#######################################################################################
source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
## Calling PPIs with strict threshold in each environment
PPI_filtering = function(coeff_matrix, best_metrics_threshold, p_threshold, p_loc){
  ### Use the coefficients that give the best_metrics_threshold
  #best_metrics_threshold = metrics_matrix[which(metrics_matrix$MCC == max(metrics_matrix$MCC)),1]

  coeff_select = as.numeric(coeff_matrix[which(coeff_matrix[,1] == best_metrics_threshold), 2: ncol(coeff_matrix)])
  
  ### Only consider PPIs with positive fitness values
  PPI_multiple_pos_fit = PPI_multiple[which(PPI_multiple[,3] > 0),] 
  fitness_multiple = PPI_multiple_pos_fit[,3]
  Q_value_multiple = PPI_multiple_pos_fit[,p_loc]
  log10_Q_value = log10(Q_value_multiple)
  log10_Q_value[log10_Q_value == -Inf] = -350
  
  ### Split PPIs into two groups based on p-values
  index_1 = which(log10_Q_value >= p_threshold)
  index_2 = which(log10_Q_value < p_threshold)
  fitness_multiple_1 = fitness_multiple[index_1]
  fitness_multiple_2 = fitness_multiple[index_2] 
  
  ### Binomial relationships between p-value and fitness
  pos_PPI = rep(0, nrow(PPI_multiple_pos_fit))
  a = log10_Q_value[index_1]
  Q_matrix_cal = coeff_select[1]/(1 + exp((coeff_select[2] - a)/coeff_select[3]))
  pos_PPI[index_1] = (fitness_multiple_1 >= Q_matrix_cal)
  pos_PPI_min = min(PPI_multiple_pos_fit[which(pos_PPI == 1), 3]) # minimum fitness 
  pos_PPI[index_2] = (fitness_multiple_2 >= pos_PPI_min)
  
  ### Get the PPIs that pass the filter
  PPI_called = PPI_multiple_pos_fit[which(pos_PPI == 1),] 
  
  ### Separate the control strains 
  PRS_PPI = PPI_called[grep("Pos", PPI_called[,1]),] 
  RRS_PPI = PPI_called[grep("Neg", PPI_called[,1]),]
  pos_PPI_DHFR = PPI_called[grep("pos", PPI_called[,1]),]
  pos_PPI_control = rbind(PRS_PPI, RRS_PPI, pos_PPI_DHFR)
  csvWriter(pos_PPI_control, "Pos_PPI_control.csv")
  
  ### Separate the PPIs containing promiscuous proteins
  fragment_PPI = PPI_called[grep("HO:", PPI_called[,1]),]
  csvWriter(fragment_PPI, "Promiscuous_proteins.csv")
  
  pos_PPI_real = PPI_called[which(!PPI_called[,1] %in% pos_PPI_control[,1]),] 
  
  ### Remove PPIs that contain promiscuous proteins
  if (nrow(fragment_PPI) >= 1){
      fragment_protein = unique(as.vector(split_string_vector(fragment_PPI[,1])))
      control_select = rep(0, nrow(pos_PPI_real))
      for(i in 1: nrow(pos_PPI_real)){
          PPI = split_string(pos_PPI_real[i,1])
          if(PPI[1] %in% fragment_protein | PPI[2] %in% fragment_protein){
              control_select[i] = 1
          }
      }
      pos_PPI_real_filter = pos_PPI_real[which(control_select == 0),] 
  }else{
    pos_PPI_real_filter = pos_PPI_real
  }
  csvWriter(pos_PPI_real_filter, "Pos_PPI_real.csv")
  return (best_metrics_threshold) ## return the best threshold 
}

FPR_threshold = read.csv("/Volumes/zmliu_02/PPiseq_03/strict_threshold/FPR_threshold_environments.csv")
## SD_merge
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/DMSO_merge/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/reference_set/p_value_2/coefficient_matrix.csv")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "SD"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.79

## FK506
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/FK506/counts/")
PPI_multiple = read.csv("PPI_multiple_p.values.csv", header = T, sep = ",")
coeff_matrix = read.csv("/Volumes/zmliu_02/PPiseq_03/FK506/reference_set/p_value/coefficient_matrix.csv", header = T, sep = ",")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "FK506"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.78

## H2O2
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/H2O2/counts/")
PPI_multiple = read.csv("PPI_multiple_p.values.csv", header = T, sep = ",")
coeff_matrix = read.csv("/Volumes/zmliu_02/PPiseq_03/H2O2/reference_set/p_value/coefficient_matrix.csv", header = T, sep = ",")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "H2O2"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.8

## Hydroxyurea
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/HU/counts/")
PPI_multiple = read.csv("PPI_multiple_p.values.csv", header = T, sep = ",")
coeff_matrix = read.csv("/Volumes/zmliu_02/PPiseq_03/HU/reference_set/p_value/coefficient_matrix.csv", header = T, sep = ",")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "Hydroxyurea"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.8

## NaCl
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/NaCl/counts/")
PPI_multiple = read.csv("PPI_multiple_p.values.csv", header = T, sep = ",")
coeff_matrix = read.csv("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/reference_set/p_value/coefficient_matrix.csv", header = T, sep = ",")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "NaCl"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.76

## Forskolin
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Forskolin/counts/")
PPI_multiple = read.csv("PPI_multiple_p.values.csv", header = T, sep = ",")
coeff_matrix = read.csv("/Volumes/zmliu_02/PPiseq_03/Forskolin/reference_set/p_value/coefficient_matrix.csv", header = T, sep = ",")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "Forskolin"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.77

## Raffinose
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Raffinose/counts/")
PPI_multiple = read.csv("PPI_multiple_p.values.csv", header = T, sep = ",")
coeff_matrix = read.csv("/Volumes/zmliu_02/PPiseq_03/Raffinose/reference_set/p_value/coefficient_matrix.csv", header = T, sep = ",")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "Raffinose"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.56

## Dox
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Dox/counts/")
PPI_multiple = read.csv("PPI_multiple_p.values.csv", header = T, sep = ",")
coeff_matrix = read.csv("/Volumes/zmliu_02/PPiseq_03/Dox/reference_set/p_value/coefficient_matrix.csv", header = T, sep = ",")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "Doxorubicin"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.82

## 16C
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/16C/counts/")
PPI_multiple = read.csv("PPI_multiple_p.values.csv", header = T, sep = ",")
coeff_matrix = read.csv("/Volumes/zmliu_02/PPiseq_03/16C/reference_set/p_value/coefficient_matrix.csv", header = T, sep = ",")
best_threshold = FPR_threshold[which(FPR_threshold$Condition == "Cold_16C"),2]
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, best_threshold, p_threshold, p_loc) # 0.51

###########################################################
## Summarize the promiscuous proteins from different environments and remove PPIs containing these proteins
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Combine_environments_SD_merge/Promiscuous_PPIs/")
DMSO = csvReader_T("SD_merge_Promiscuous_proteins.csv")
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
H2O2_bad = extract_promiscuous_protein(H2O2)
HU_bad = extract_promiscuous_protein(HU)
Dox_bad = extract_promiscuous_protein(Dox)
#Forskolin_bad = extract_promiscuous_protein(Forskolin)
Raffinose_bad = extract_promiscuous_protein(Raffinose)
NaCl_bad = extract_promiscuous_protein(NaCl)
FK506_bad = extract_promiscuous_protein(FK506)
cold_bad = extract_promiscuous_protein(cold)

final_bad = unique(c(DMSO_bad, H2O2_bad, HU_bad, Dox_bad, 
                     Raffinose_bad, NaCl_bad, FK506_bad, cold_bad)) # 22

PPI_list = list(DMSO_bad, H2O2_bad, HU_bad, Dox_bad, 
                Raffinose_bad, NaCl_bad, FK506_bad, cold_bad)
matrix_bad = matrix(0, length(final_bad), 10)
colnames(matrix_bad) = c("PPI","Positive_sum","SD_merge","H2O2", "HU", "Dox", 
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
csvWriter(matrix_bad_order, "Promiscuous_protein_summary_SD_merge.csv")

setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Combine_environments_SD_merge/")
DMSO_pos = csvReader_T("PPI_each_environment/SD_merge_Pos_PPI_real.csv") # 2382
H2O2_pos = csvReader_T("PPI_each_environment/H2O2_Pos_PPI_real.csv") # 2075
HU_pos = csvReader_T("PPI_each_environment/Hydroxyurea_Pos_PPI_real.csv") # 2211
Dox_pos = csvReader_T("PPI_each_environment/Doxorubicin_Pos_PPI_real.csv") # 2060
Forskolin_pos = csvReader_T("PPI_each_environment/Forskolin_Pos_PPI_real.csv") # 1914
NaCl_pos = csvReader_T("PPI_each_environment/NaCl_Pos_PPI_real.csv") # 1952
Raffinose_pos = csvReader_T("PPI_each_environment/Raffinose_Pos_PPI_real.csv") # 2586
FK506_pos = csvReader_T("PPI_each_environment/FK506_Pos_PPI_real.csv") # 1984
Cold_pos = csvReader_T("PPI_each_environment/Cold_16C_Pos_PPI_real.csv") # 2146
promiscuous_protein = csvReader_T("Promiscuous_PPIs/Promiscuous_protein_summary_SD_merge.csv")
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
DMSO_pos_filter = remove_PPI_protein(DMSO_pos, promiscuous_protein_chosen) # 2274 from 2382
H2O2_pos_filter = remove_PPI_protein(H2O2_pos, promiscuous_protein_chosen) # 2069 from 2075
HU_pos_filter = remove_PPI_protein(HU_pos, promiscuous_protein_chosen) # 2178 from 2211
Dox_pos_filter = remove_PPI_protein(Dox_pos, promiscuous_protein_chosen) # 1966 from 2060
Forskolin_pos_filter = remove_PPI_protein(Forskolin_pos, promiscuous_protein_chosen) # 1696 from 1914
NaCl_pos_filter = remove_PPI_protein(NaCl_pos, promiscuous_protein_chosen) # 1865 from 1952
Raffinose_pos_filter = remove_PPI_protein(Raffinose_pos, promiscuous_protein_chosen) # 2527 from 2586
FK506_pos_filter = remove_PPI_protein(FK506_pos, promiscuous_protein_chosen) # 1772 from 1984
cold_pos_filter = remove_PPI_protein(Cold_pos, promiscuous_protein_chosen) # 1710 from 2146

csvWriter(DMSO_pos_filter, "Positive_PPI_remove_promiscuous/SD_merge_Pos_PPI_real.csv")
csvWriter(H2O2_pos_filter, "Positive_PPI_remove_promiscuous/H2O2_Pos_PPI_real.csv")
csvWriter(HU_pos_filter, "Positive_PPI_remove_promiscuous/Hydroxyurea_Pos_PPI_real.csv")
csvWriter(Dox_pos_filter, "Positive_PPI_remove_promiscuous/Doxorubicin_Pos_PPI_real.csv")
csvWriter(Forskolin_pos_filter, "Positive_PPI_remove_promiscuous/Forskolin_Pos_PPI_real.csv")
csvWriter(NaCl_pos_filter, "Positive_PPI_remove_promiscuous/NaCl_Pos_PPI_real.csv")
csvWriter(Raffinose_pos_filter, "Positive_PPI_remove_promiscuous/Raffinose_Pos_PPI_real.csv")
csvWriter(FK506_pos_filter, "Positive_PPI_remove_promiscuous/FK506_Pos_PPI_real.csv")
csvWriter(cold_pos_filter, "Positive_PPI_remove_promiscuous/Cold_16C_Pos_PPI_real.csv")

######################################################################################################
### Generate a count summary for each PPI (single_orientation) across different environments
#### If a PPI has two oritentations and either of them (or both) is positive in a environment, that PPI is positive in that environment
# The basic idea of the code is to attach positive environments to each PPI

setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Combine_environments_SD_merge/Positive_PPI_remove_promiscuous")

DMSO_pos = csvReader_T("SD_merge_Pos_PPI_real.csv")
H2O2_pos = csvReader_T("H2O2_Pos_PPI_real.csv")
HU_pos = csvReader_T("Hydroxyurea_Pos_PPI_real.csv")
Dox_pos = csvReader_T("Doxorubicin_Pos_PPI_real.csv")
Forskolin_pos = csvReader_T("Forskolin_Pos_PPI_real.csv")
Raffinose_pos = csvReader_T("Raffinose_Pos_PPI_real.csv")
NaCl_pos = csvReader_T("NaCl_Pos_PPI_real.csv")
cold_pos = csvReader_T("Cold_16C_Pos_PPI_real.csv")
FK506_pos = csvReader_T("FK506_Pos_PPI_real.csv")

PPI_list = list(DMSO_pos[,1],  H2O2_pos[,1], HU_pos[,1], Dox_pos[,1], Forskolin_pos[,1],
                Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1]) # store PPIs of one environment as an element in a list
all_PPI = unique(c(DMSO_pos[,1],  H2O2_pos[,1], HU_pos[,1], Dox_pos[,1], Forskolin_pos[,1],
                   Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1])) #6390
all_PPI_unique = mark_duplicates_fast(all_PPI) # 6071; 319 duplicates
all_PPI_matrix = matrix(0, nrow(all_PPI_unique), 10)
all_PPI_matrix[,1] = all_PPI_unique[,1]
colnames(all_PPI_matrix)= c("PPI","SD", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
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
csvWriter(all_PPI_matrix_final, "PPI_environment_count_summary_SD_merge.csv")
### Remove single
PPI_count = csvReader_T("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Combine_environments_SD_merge/Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge.csv")
PPI_SD = PPI_count[which(PPI_count[,3] == "1"),]
SD1_single = csvReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO/counts/PPI_1_p.values.csv") # 189787
SD2_single = csvReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_2/counts/PPI_1_p.values.csv") # 168705
common_single = intersect(SD1_single[,1], SD2_single[,1]) # 100412
SD_pos_single = intersect(common_single, PPI_SD[,1]) # 290
PPI_count_filter = PPI_count[which(!PPI_count[,1] %in% SD_pos_single),] # 5928
csvWriter(PPI_count_filter, "/Volumes/zmliu_02/PPiseq_03/strict_threshold/Combine_environments_SD_merge/Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge_filter_strict_threshold.csv")

### Move the PPI_count_filter file into dropbox PPiSeq_02/Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter_strict_threshold.csv

################################################################ Calculate the variation score for each PPI
##### (3) Combine different environments and put the same PPI across different environments on the same row
## The normalized fitness values are the same with previous calculation
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Combine_environments_SD_merge/")
DMSO_matrix = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/SD_merge_normalized_multiple.csv")
Forskolin_norm = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/Forskolin_normalized_multiple.csv")
FK506_norm = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/FK506_normalized_multiple.csv")
Raffinose_norm = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/Raffinose_normalized_multiple.csv")
NaCl_norm = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/NaCl_normalized_multiple.csv")
H2O2_norm = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/H2O2_normalized_multiple.csv")
Dox_norm = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/Doxorubicin_normalized_multiple.csv")
cold_norm = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/Cold_16C_normalized_multiple.csv")
HU_norm = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Normalized_multiple_files/Hydroxyurea_normalized_multiple.csv")

### Only consider the positive PPIs detected in any environment
DMSO_real = csvReader_T("Positive_PPI_remove_promiscuous/SD_merge_Pos_PPI_real.csv")  
Forskolin_real = csvReader_T("Positive_PPI_remove_promiscuous/Forskolin_Pos_PPI_real.csv")
FK506_real = csvReader_T("Positive_PPI_remove_promiscuous/FK506_Pos_PPI_real.csv")
Raffinose_real = csvReader_T("Positive_PPI_remove_promiscuous/Raffinose_Pos_PPI_real.csv")
NaCl_real = csvReader_T("Positive_PPI_remove_promiscuous/NaCl_Pos_PPI_real.csv")
H2O2_real = csvReader_T("Positive_PPI_remove_promiscuous/H2O2_Pos_PPI_real.csv")
Dox_real = csvReader_T("Positive_PPI_remove_promiscuous/Doxorubicin_Pos_PPI_real.csv")
cold_real = csvReader_T("Positive_PPI_remove_promiscuous/Cold_16C_Pos_PPI_real.csv")
HU_real = csvReader_T("Positive_PPI_remove_promiscuous/Hydroxyurea_Pos_PPI_real.csv")


all_pos = unique(c(DMSO_real[,1], Forskolin_real[,1], FK506_real[,1], Raffinose_real[,1], 
                   NaCl_real[,1], H2O2_real[,1], Dox_real[,1], cold_real[,1], HU_real[,1])) #6390

# Keep the primary normalized fitness values for PPIs even not detected in an evironment
matrix_pos = matrix(0, length(all_pos), 10)
matrix_pos[,1] = all_pos
matrix_pos[,2] = as.numeric(DMSO_matrix[match(all_pos, DMSO_matrix[,1]),3])
matrix_pos[,3] = as.numeric(H2O2_norm[match(all_pos, H2O2_norm[,1]),3])
matrix_pos[,4] = as.numeric(HU_norm[match(all_pos, HU_norm[,1]),3])
matrix_pos[,5] = as.numeric(Dox_norm[match(all_pos, Dox_norm[,1]),3])
matrix_pos[,6] = as.numeric(Forskolin_norm[match(all_pos, Forskolin_norm[,1]),3])
matrix_pos[,7] = as.numeric(Raffinose_norm[match(all_pos, Raffinose_norm[,1]),3])
matrix_pos[,8] = as.numeric(NaCl_norm[match(all_pos, NaCl_norm[,1]),3])
matrix_pos[,9] = as.numeric(cold_norm[match(all_pos, cold_norm[,1]),3])
matrix_pos[,10] = as.numeric(FK506_norm[match(all_pos, FK506_norm[,1]),3])

colnames(matrix_pos) = c("PPI", "SD", "H2O2", "HU", "Dox", "Forskolin", 
                         "Raffinose", "NaCl", "16C", "FK506")
csvWriter(matrix_pos, "Normalized_multiple_files/Pos_PPI_normalized_fit_primary_SD_merge_strict_threshold.csv")

### (4) Put the fitness of PPI of two orientations across different environments onto the same row
setwd("/Volumes/zmliu_02/PPiseq_03/strict_threshold/Combine_environments_SD_merge/")
PPI_fit_norm = csvReader_T("Normalized_multiple_files/Pos_PPI_normalized_fit_primary_SD_merge_strict_threshold.csv") # 6390
## Make the PPI name for the count and normalized fitness consistent
PPI_count = csvReader_T("Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge.csv") # 6071
PPI_dup = mark_duplicates_fast(PPI_fit_norm[,1]) # 6218
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
    matrix[i,1:10] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
    matrix[i,11:20] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,2]),]
  }else{
    matrix[i,1:10] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
  }   
  
}
colnames(matrix) = c("PPI", "SD", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", 
                     "NaCl", "16C", "FK506", "PPI_opposite", "SD", "H2O2", "HU",
                     "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
mean_fitness = rep(0, nrow(matrix))
for(i in 1:length(mean_fitness)){
  mean_fitness[i] = mean(na.omit(as.numeric(matrix[i, c(2:10, 12:20)])))
}

matrix_final = cbind(matrix[,1], mean_fitness, matrix[,2:ncol(matrix)])
csvWriter(matrix_final, "Normalized_multiple_files/Normalized_fitness_PPI_all_primary_SD_merge.csv")

########### Consider any negative fitness values of any orientation to be zero

PPI_norm = csvReader_T("Normalized_multiple_files/Normalized_fitness_PPI_all_primary_SD_merge.csv")
PPI_count = csvReader_T("Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge.csv")
PPI_norm_matrix = matrix(0, nrow(PPI_norm), 10)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
  for(j in 3:11){
    a = as.numeric(PPI_norm[i, j])
    b = as.numeric(PPI_norm[i, j + 10])
    if ((!is.na(a)) & a < 0){ ### Consider negative values to be zero
      a = 0
    }
    if ((!is.na(b)) & b < 0){
      b = 0
    }
    mean_fit = mean(as.numeric(na.omit(c(a,b))))
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
  fitness = as.numeric(PPI_norm_matrix[i,2:10])
  variation_score[i] = sd(fitness)/mean(fitness)
}
environment_number = as.numeric(PPI_count[match(PPI_norm_matrix[,1], PPI_count[,1]),2])
# Include the variation score into the matrix, and output the matrix
PPI_norm_matrix_final = cbind(PPI_norm_matrix[,1],environment_number, variation_score,
                              PPI_norm_matrix[,2:10])
colnames(PPI_norm_matrix_final) = c("PPI", "Environment_number", "Variation_score",
                                    "SD", "H2O2", "HU", "Dox", "Forskolin", "Raffinose",
                                    "NaCl", "16C", "FK506")

csvWriter(PPI_norm_matrix_final, "Normalized_multiple_files/Variation_score_PPI_environment_neg_zero_SD_merge.csv")
cor(environment_number, variation_score, method="spearman") # -0.7404101

### Remove single barcoded PPI in SD environment
PPI_count_filter = csvReader_T("Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge_filter_strict_threshold.csv") 
PPI_vScore = csvReader_T("Normalized_multiple_files/Variation_score_PPI_environment_neg_zero_SD_merge.csv")
PPI_vScore_filter = PPI_vScore[which(PPI_vScore[,1] %in% PPI_count_filter[,1]),] # 5781
csvWriter(PPI_vScore_filter, "Normalized_multiple_files/Variation_score_PPI_environment_neg_zero_SD_merge_filter_strict_threshold.csv")

### Move the PPI_vScore_filter file into dropbox PPiSeq_02/Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter_strict_threshold.csv


