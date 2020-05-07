source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

Fitness_normalization = function(DMSO_lineage){
  DMSO_DHFR_Pos= DMSO_lineage[which(DMSO_lineage[,1] == "positive_DHFR"),]
  DMSO_DHFR_Neg= DMSO_lineage[which(DMSO_lineage[,1] == "negative_non_DHFR"),]
  
  DMSO_DHFR_Pos_mean = mean(as.numeric(DMSO_DHFR_Pos[,4])) 
  DMSO_DHFR_Neg_mean = mean(as.numeric(DMSO_DHFR_Neg[,4])) 
  
  DMSO_lineage[,4] = (as.numeric(DMSO_lineage[,4]) - DMSO_DHFR_Neg_mean)/(DMSO_DHFR_Pos_mean - DMSO_DHFR_Neg_mean)
  return(DMSO_lineage)
}

########################## (1) Normalize the fitness values since fitness range of two environments are not the same
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments/")
SD_lineage = "lineage_fitness_files/SD_PPI_barcodes_fitness_counts.csv"
SD_fit = csvReader_T(SD_lineage)
SD_fit_norm = Fitness_normalization(SD_fit)

SD2_lineage = "lineage_fitness_files/SD2_PPI_barcodes_fitness_counts.csv"
SD2_fit = csvReader_T(SD2_lineage)
SD2_fit_norm = Fitness_normalization(SD2_fit)

SD_all = rbind(SD_fit_norm, SD2_fit_norm)
SD_all[,2] = SD_all[,1]
SD_all_final = cbind(SD_all[,1:3], rep("DMSO", nrow(SD_all)), as.numeric(SD_all[,4]), 
                     as.numeric(SD_all[,5]),  as.numeric(SD_all[,6]), as.numeric(SD_all[,7]), 
                     as.numeric(SD_all[,8]),  as.numeric(SD_all[,9]), as.numeric(SD_all[,10]))

colnames(SD_all_final) = c("Systematic_PPI", "Standard_PPI", "barcode", "Environment", 
                           "Fitness", "Distance", "G0", "G3", "G6", "G9", "G12")
csvWriter(SD_all_final, "/Volumes/zmliu_02/PPiseq_03/DMSO_merge/counts/SD_merged_PPI_barcodes.csv")

######################## (2) Python code (Sort_PPI) to sort PPIs (put the same PPI right next to each other)

######################## (3) Calculate the p-value and mean fitness for each PPI
p_value_calculate = function(PPI_fitness, P_value_output_01, P_value_output_02){
  PPI_indiv= csvReader_F(PPI_fitness)
  MATa_DHFR12 = csvReader_T("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/PPI_barcodes/MATa_genome_combine.csv")
  MATalpha_DHFR3 = csvReader_T("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/PPI_barcodes/MATalpha_genome_combine.csv")
  PPI_DHFR12 = PPI_indiv[which(PPI_indiv[,3] %in% MATa_DHFR12[,3]),]
  PPI_DHFR3 = PPI_indiv[which(PPI_indiv[,3] %in% MATalpha_DHFR3[,3]),]
  #PPI_negative_DHFR= PPI_indiv[which(PPI_indiv[,1] == "negative_non_DHFR"),] 
  #PPI_negative = rbind(PPI_DHFR12, PPI_DHFR3, PPI_negative_DHFR)
  PPI_negative = rbind(PPI_DHFR12, PPI_DHFR3) 
  PPI_negative_fitness= as.numeric(PPI_negative[,4])
  # plot histrogram of fitness values of negative control group
  pdf("Distribution of negative fitness.pdf", width = 5, height = 5)
  hist(PPI_negative_fitness, breaks = seq(-1.5, 1.5, by = 0.05), xlim = c(-0.6, 0.6),
       xlab = "Fitness of DHFR fragment strains", ylab = "Frequency", main = NA)
  dev.off()
  
  ### Create a matrix containing PPI, barcode count, and other metrics
  PPI_unique= unique(PPI_indiv[,1])
  matrix_summary= matrix(0, length(PPI_unique), 7)
  matrix_summary[,1]= PPI_unique
  colnames(matrix_summary)= c("PPI", "barcode_counts", "mean", "sd", "t", "p.value", "FDR_adjusted_value")
  b=1
  for (i in 1: nrow(matrix_summary)){
    a= as.numeric(PPI_indiv[b,2])
    if(a > 1){
      matrix_summary[i,2]= a
      matrix_summary[i,3]= mean(as.numeric(PPI_indiv[b:(b+a-1), 4]))
      matrix_summary[i,4]= sd(as.numeric(PPI_indiv[b:(b+a-1),4]))
      t_test = t.test(as.numeric(PPI_indiv[(b:(b+a-1)),4]), PPI_negative_fitness, alternative = "greater")
      matrix_summary[i,5]= t_test$statistic
      matrix_summary[i,6]= t_test$p.value
      b= b+a
    }
    else{
      matrix_summary[i, 2]= a
      matrix_summary[i, 3]= as.numeric(PPI_indiv[b,4])
      matrix_summary[i, 4]= "NA"
      matrix_summary[i, 5]= "NA"
      matrix_summary[i, 6] = "NA"
      b= b+a
    }
    
  }
  matrix_PPI_01= matrix_summary[which(matrix_summary[,5] == "NA"),]
  matrix_filtered= matrix_summary[which(!matrix_summary[,5] == "NA"),]
  matrix_filtered[,7]= p.adjust(as.numeric(matrix_filtered[,6]), method= "fdr")
  
  csvWriter(matrix_PPI_01, P_value_output_01)
  csvWriter(matrix_filtered, P_value_output_02)
  
}

setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

################################ (4) Based on the data generate reference sets
yeast_PPI = read.delim("/Volumes/zmliu_02/PPiseq/DMSO_2/reference_set/BIOGRID-ORGANISM-3.5.165.tab2/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.165.tab2.txt", header = T)
pos_PPI= yeast_PPI[which(yeast_PPI$Experimental.System.Type == "physical"), c(6:7, 13)] # 164992
# generate positive PPIs with two directions (bait-prey or prey-bait)
pos_PPI_01= paste(pos_PPI[,1], pos_PPI[,2], sep="_")
pos_PPI_02= paste(pos_PPI[,2], pos_PPI[,1], sep="_")
pos_PPI_both= unique(c(pos_PPI_01, pos_PPI_02)) # 215884

### Gene included in PCA collection
yeast_gene= as.matrix(read.csv('/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/reference_set/Yeast_genes.csv', header = F)) # PCA library doe not contain all the genes
yeast_PPI= expand.grid(yeast_gene[,1], yeast_gene[,1])# 32,216,976 two directions
yeast_PPI_filtered = yeast_PPI[which(yeast_PPI$Var1 != yeast_PPI$Var2),] # 32211300, remove sefl-interacting PPIs
All_PPI= paste(yeast_PPI_filtered[,1], yeast_PPI_filtered[,2], sep="_")
all_PPI_filtered = All_PPI[which(!All_PPI %in% pos_PPI_both)] #  32009066 filter the positive PPIs (193364)

### PPI summary generated by BioGRID
PPI_summary = dataFrameReader_T("/Volumes/zmliu_02/PPiseq/H2O2/reference_set/one_binary/multiple_validated_PPI.csv")

generating_reference_sets = function(all_PPI_filtered, PPI_summary, PPI_multiple, 
                                     number_PPI, sample_number, data_set_number){
  Random_reference <- function(PPI_multiple, all_PPI_filtered, size, sampling_number){
    PPI_PPiseq = PPI_multiple[,1]
    RRS_size= 0
    RRS= character(length=0)
    while (RRS_size < size) {
      yeast_PPI_random= all_PPI_filtered[sample(1: length(all_PPI_filtered), sampling_number, replace=F)]
      RRS= unique(c(RRS, yeast_PPI_random))
      RRS_overlap= intersect(RRS, PPI_PPiseq)
      RRS_size= length(RRS_overlap)
    } 
    RRS_duplicate_marked = mark_duplicates_fast(RRS_overlap)
    return (RRS_duplicate_marked[,1])
  }
  
  Random_reference_generation = function(PPI_multiple, all_PPI_filtered, number_PPI, sample_number, data_set_number){
    neg_set = vector("list", data_set_number)
    for(i in 1:data_set_number){
      PPI_neg = Random_reference(PPI_multiple, all_PPI_filtered, number_PPI, sample_number)
      PPI_neg_matrix = PPI_multiple[which(PPI_multiple[,1] %in% PPI_neg),]
      neg_set[[i]] = PPI_neg_matrix
    }
    return (neg_set)
  }
  
  PPI_3_assay= PPI_summary[which(PPI_summary$methods >= 3 & PPI_summary$publications >= 3 & PPI_summary$Binary_assy == 1),1] # 1968
  PPI_pos_3_assay = match_both_direction(PPI_multiple, PPI_3_assay) # 633
  print(nrow(PPI_pos_3_assay)) #633
  pos_set = list(PPI_pos_3_assay)
  csvWriter(PPI_pos_3_assay[,1], "Positive_reference_set.csv")
  
  neg_set = Random_reference_generation (PPI_multiple, all_PPI_filtered, number_PPI, sample_number, data_set_number)
  
  pos_name = "PPI_pos_3_assay"
  neg_name = vector("character", data_set_number)
  for (i in 1:length(neg_name)){
    neg_name[i] = paste("neg", as.character(number_PPI), as.character(i), sep = "_")
  }
  
  for (i in 1:length(pos_set)){
    for (j in 1:length(neg_set)){
      a = as.matrix(pos_set[[i]])
      b = as.matrix(neg_set[[j]])
      csvWriter(b[,1], paste(neg_name[j],"reference.csv", sep = "_"))
      matrix = rbind(a,b)
      matrix[1:nrow(a), 1] = "Positive"
      matrix[(nrow(a) +1): (nrow(a) + nrow(b)), 1] = "Negative"
      colnames(matrix) = colnames(PPI_multiple)
      csvWriter(matrix, paste(pos_name[i], neg_name[j],"reference.csv", sep= "_"))
    }
  }
}
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/reference_set/p_value")
number_PPI = 6.4e4
sample_number = 6.4e4
data_set_number = 50
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/counts/PPI_multiple_p.values.csv") 
generating_reference_sets(all_PPI_filtered, PPI_summary, PPI_multiple, number_PPI, sample_number, data_set_number)

############################# (5) Defining thresholds to call PPIs
## import the PPI_calling_sigmoid from 6_threshold_script
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/reference_set/p_value/")
Neg_number_PPI = 6.4e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.73, by = 0.01))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)
##### test a different PPV range
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/reference_set/p_value_2/")
Neg_number_PPI = 6.4e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.73, by = 0.01), c(0.79, 0.80, 0.84, 0.86))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

####### Output the coefficient matrix
PPV_file_name = rep("A", length(specific_PPV))
for(i in 1:length(specific_PPV)){
  PPV_file_name[i] = paste("PPV", as.character(specific_PPV[i]), "threshold.csv", sep= "_")
}
coeff_calling = function(PPV_file_name, PPV){
  coeff_matrix = matrix(0, length(PPV_file_name), 4)
  coeff_matrix[,1] = PPV
  for (i in 1:length(PPV_file_name)){
    file_name= PPV_file_name[i]
    PPV_data = dataFrameReader_T(file_name)
    fitness_threshold = PPV_data[,3]
    Q_value_threshold = PPV_data[,2]
    fitmodel_2 = nls(fitness_threshold ~ SSlogis(Q_value_threshold, Asym, xmid, scal))
    coeff_matrix[i, 2: ncol(coeff_matrix)] = coef(fitmodel_2)
  }
  colnames(coeff_matrix) = c("threshold", "asym", "xmid", "scal")
  csvWriter(coeff_matrix, "coefficient_matrix.csv")
}
## write coefficients into a matrix
coeff_calling(PPV_file_name, specific_PPV)
coeff_matrix = csvReader_T("coefficient_matrix.csv")

############################# (6) Calling PPIs
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.7
















############################# (7) Transform data and summarize after filtering promiscuous
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/lineage_fitness_files/")
DMSO = dataFrameReader_F("SD_merge_known_PPI_fitness_barcodes_sorted.csv")
colnames(DMSO) = c("PPI", "Number_of_Barcodes", "Barcode_sequences", "Fitness", "Fitness_estimaion_error", 
                   "Counts_G0", "Counts_G3", "Counts_G6", "Counts_G9", "Counts_G12")
csvWriter(DMSO, "SD_merge__PPI_barcodes_fitness_counts.csv")

#### Summarize the promiscuous proteins and remove all PPIs that contain promiscuous proteins reported in two environments
#### Check the number of environments in which a promiscuous protein is detected
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Promiscuous_PPIs/")
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
Forskolin_bad = extract_promiscuous_protein(Forskolin)
Raffinose_bad = extract_promiscuous_protein(Raffinose)
NaCl_bad = extract_promiscuous_protein(NaCl)
FK506_bad = extract_promiscuous_protein(FK506)
cold_bad = extract_promiscuous_protein(cold)

final_bad = unique(c(DMSO_bad, H2O2_bad, HU_bad, Dox_bad, Forskolin_bad,
                     Raffinose_bad, NaCl_bad, FK506_bad, cold_bad)) # 51

PPI_list = list(DMSO_bad, H2O2_bad, HU_bad, Dox_bad, Forskolin_bad,
                Raffinose_bad, NaCl_bad, FK506_bad, cold_bad)
matrix_bad = matrix(0, length(final_bad), 11)
colnames(matrix_bad) = c("PPI","Positive_sum","SD_merge","H2O2", "HU", "Dox", "Forskolin",
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
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/")
DMSO_pos = csvReader_T("PPI_each_environment/SD_merge_Pos_PPI_real.csv") # 5059
H2O2_pos = csvReader_T("PPI_each_environment/H2O2_Pos_PPI_real.csv") # 4341
HU_pos = csvReader_T("PPI_each_environment/Hydroxyurea_Pos_PPI_real.csv") # 4967
Dox_pos = csvReader_T("PPI_each_environment/Doxorubicin_Pos_PPI_real.csv") # 3387
Forskolin_pos = csvReader_T("PPI_each_environment/Forskolin_Pos_PPI_real.csv") # 6103
NaCl_pos = csvReader_T("PPI_each_environment/NaCl_Pos_PPI_real.csv") # 2482
Raffinose_pos = csvReader_T("PPI_each_environment/Raffinose_Pos_PPI_real.csv") # 4452
FK506_pos = csvReader_T("PPI_each_environment/FK506_Pos_PPI_real.csv") # 3130
Cold_pos = csvReader_T("PPI_each_environment/Cold_16C_Pos_PPI_real.csv") # 5558
promiscuous_protein = csvReader_T("Promiscuous_PPIs/Promiscuous_protein_summary.csv")
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
DMSO_pos_filter = remove_PPI_protein(DMSO_pos, promiscuous_protein_chosen) # 5041 from 5059
H2O2_pos_filter = remove_PPI_protein(H2O2_pos, promiscuous_protein_chosen) # 4331 from 4341
HU_pos_filter = remove_PPI_protein(HU_pos, promiscuous_protein_chosen) # 4255 from 4967
Dox_pos_filter = remove_PPI_protein(Dox_pos, promiscuous_protein_chosen) # 3187 from 3387
Forskolin_pos_filter = remove_PPI_protein(Forskolin_pos, promiscuous_protein_chosen) # 5356 from 6103
NaCl_pos_filter = remove_PPI_protein(NaCl_pos, promiscuous_protein_chosen) # 2374 from 2482
Raffinose_pos_filter = remove_PPI_protein(Raffinose_pos, promiscuous_protein_chosen) # 4442 from 4452
FK506_pos_filter = remove_PPI_protein(FK506_pos, promiscuous_protein_chosen) # 3004 from 3130
cold_pos_filter = remove_PPI_protein(Cold_pos, promiscuous_protein_chosen) # 4586 from 5558

csvWriter(DMSO_pos_filter, "Positive_PPI_remove_promiscuous/SD_merge_Pos_PPI_real.csv")
csvWriter(H2O2_pos_filter, "Positive_PPI_remove_promiscuous/H2O2_Pos_PPI_real.csv")
csvWriter(HU_pos_filter, "Positive_PPI_remove_promiscuous/Hydroxyurea_Pos_PPI_real.csv")
csvWriter(Dox_pos_filter, "Positive_PPI_remove_promiscuous/Doxorubicin_Pos_PPI_real.csv")
csvWriter(Forskolin_pos_filter, "Positive_PPI_remove_promiscuous/Forskolin_Pos_PPI_real.csv")
csvWriter(NaCl_pos_filter, "Positive_PPI_remove_promiscuous/NaCl_Pos_PPI_real.csv")
csvWriter(Raffinose_pos_filter, "Positive_PPI_remove_promiscuous/Raffinose_Pos_PPI_real.csv")
csvWriter(FK506_pos_filter, "Positive_PPI_remove_promiscuous/FK506_Pos_PPI_real.csv")
csvWriter(cold_pos_filter, "Positive_PPI_remove_promiscuous/Cold_16C_Pos_PPI_real.csv")

# Make a file that contain mean fitness, sd, t, p.value, Q.value, positive (1: positive PPIs including controls)
# The promiscuous proteins did not change therefore I do not need relabel PPIs in each environment
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/")
Combine_matrix = function(DMSO, DMSO_pos_control, DMSO_pos_real, output){
  DMSO_pos_PPI = c(DMSO_pos_control[,1], DMSO_pos_real[,1])
  pos_index = rep(0, nrow(DMSO))
  pos_index[which(as.character(DMSO[,1]) %in% DMSO_pos_PPI)] = 1
  length(which(pos_index == 1))
  DMSO_label = cbind(DMSO, pos_index)
  colnames(DMSO_label) = c("PPI", "Number_of_barcodes", "Mean_fitness","sd", "t", "P_value", "FDR_adjusted_P_value", "Positive")
  csvWriter(DMSO_label, output)
}

DMSO = dataFrameReader_T("PPI_multiple_files/SD_merge_PPI_multiple_p.values.csv")
DMSO_pos_control = csvReader_T("PPI_each_environment/SD_merge_Pos_PPI_control.csv")
DMSO_pos_real = csvReader_T("Positive_PPI_remove_promiscuous/SD_merge_Pos_PPI_real.csv")
output = "PPI_multiple_files/SD_merge_mean_fitness_positive.csv"
Combine_matrix(DMSO, DMSO_pos_control, DMSO_pos_real, output)

###################### (8) Generate a count summary for each PPI (single_orientation) across different environments
#### If a PPI has two oritentations and either of them (or both) is positive in a environment, that PPI is positive in that environment
# The basic idea of the code is to attach positive environments to each PPI 
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Positive_PPI_remove_promiscuous")
DMSO_pos = csvReader_T("SD_merge_Pos_PPI_real.csv") # 5041
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
                   Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1])) #14463
all_PPI_unique = mark_duplicates_fast(all_PPI) # 13764
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
#csvWriter(all_PPI_matrix_final, "Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")

#############################(9) contain 5 steps of analysis: 
##########  (1) Normalize fitness of each lineage by DHFR(+) and DHFR(-) control strains (mean) (already done)
##########  (2) Put normized fitness values of the same PPI onto the same row
##########  (3) Combine normalized mean fitness values of different environments onto the same row
#########   (4) Put the fitness of PPI of two orientations across different environments onto the same row
#########   (5) Take the mean fitness (remove NA) for the same PPI with two orientations and calculate the variation score

### (1) and (2) normalize each PPI and put the normalized fitness onto the same row
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/")
DMSO_multiple = csvReader_T("PPI_multiple_files/SD_merge_PPI_multiple_p.values.csv")
## Because it has already been normalized, all I need to do here is to remove control strains
PPI_multiple_RRS = DMSO_multiple[grep("Neg_PPI", DMSO_multiple[,1]),] #97
PPI_multiple_PRS = DMSO_multiple[grep("Pos_PPI", DMSO_multiple[,1]),] #108
PPI_multiple_pos = DMSO_multiple[grep("positive_DHFR", DMSO_multiple[,1]),] # 1
PPI_multiple_neg = DMSO_multiple[grep("negative_non_DHFR", DMSO_multiple[,1]),] # 1
PPI_multiple_control = rbind(PPI_multiple_RRS, PPI_multiple_PRS, PPI_multiple_pos, PPI_multiple_neg)
PPI_multiple_select = DMSO_multiple[which(!DMSO_multiple[,1] %in% PPI_multiple_control[,1]),]
csvWriter(PPI_multiple_select, "Normalized_multiple_files/SD_merge_normalized_multiple.csv")

DMSO_matrix = csvReader_T("Normalized_multiple_files/SD_merge_normalized_multiple.csv")
Forskolin_norm = csvReader_T("Normalized_multiple_files/Forskolin_normalized_multiple.csv")
FK506_norm = csvReader_T("Normalized_multiple_files/FK506_normalized_multiple.csv")
Raffinose_norm = csvReader_T("Normalized_multiple_files/Raffinose_normalized_multiple.csv")
NaCl_norm = csvReader_T("Normalized_multiple_files/NaCl_normalized_multiple.csv")
H2O2_norm = csvReader_T("Normalized_multiple_files/H2O2_normalized_multiple.csv")
Dox_norm = csvReader_T("Normalized_multiple_files/Doxorubicin_normalized_multiple.csv")
cold_norm = csvReader_T("Normalized_multiple_files/Cold_16C_normalized_multiple.csv")
HU_norm = csvReader_T("Normalized_multiple_files/Hydroxyurea_normalized_multiple.csv")

#### Combine all the PPIs that has been barcoded twice in one of environments
all = unique(c(DMSO_matrix[,1], Forskolin_norm[,1], FK506_norm[,1], Raffinose_norm[,1], 
               NaCl_norm[,1], H2O2_norm[,1], Dox_norm[,1], cold_norm[,1], HU_norm[,1])) #1660648

matrix_all = matrix(0, length(all), 10)
matrix_all[,1] = all
matrix_all[,2] = as.numeric(DMSO_matrix[match(all, DMSO_matrix[,1]),3])
#matrix_all[,3] = as.numeric(DMSO2_norm[match(all, DMSO2_norm[,1]),3])
matrix_all[,3] = as.numeric(H2O2_norm[match(all, H2O2_norm[,1]),3])
matrix_all[,4] = as.numeric(HU_norm[match(all, HU_norm[,1]),3])
matrix_all[,5] = as.numeric(Dox_norm[match(all, Dox_norm[,1]),3])
matrix_all[,6] = as.numeric(Forskolin_norm[match(all, Forskolin_norm[,1]),3])
matrix_all[,7] = as.numeric(Raffinose_norm[match(all, Raffinose_norm[,1]),3])
matrix_all[,8] = as.numeric(NaCl_norm[match(all, NaCl_norm[,1]),3])
matrix_all[,9] = as.numeric(cold_norm[match(all, cold_norm[,1]),3])
matrix_all[,10] = as.numeric(FK506_norm[match(all, FK506_norm[,1]),3])

colnames(matrix_all) = c("PPI", "SD", "H2O2", "HU", "Dox", "Forskolin", 
                         "Raffinose", "NaCl", "16C", "FK506")
csvWriter(matrix_all, "Normalized_multiple_files/All_PPI_environments_normalized_fit.csv")

######## Only consider positive PPIs in one of environments
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/")
DMSO_real = csvReader_T("Positive_PPI_remove_promiscuous/SD_merge_Pos_PPI_real.csv") # 5410
Forskolin_real = csvReader_T("Positive_PPI_remove_promiscuous/Forskolin_Pos_PPI_real.csv")
FK506_real = csvReader_T("Positive_PPI_remove_promiscuous/FK506_Pos_PPI_real.csv")
Raffinose_real = csvReader_T("Positive_PPI_remove_promiscuous/Raffinose_Pos_PPI_real.csv")
NaCl_real = csvReader_T("Positive_PPI_remove_promiscuous/NaCl_Pos_PPI_real.csv")
H2O2_real = csvReader_T("Positive_PPI_remove_promiscuous/H2O2_Pos_PPI_real.csv")
Dox_real = csvReader_T("Positive_PPI_remove_promiscuous/Doxorubicin_Pos_PPI_real.csv")
cold_real = csvReader_T("Positive_PPI_remove_promiscuous/Cold_16C_Pos_PPI_real.csv")
HU_real = csvReader_T("Positive_PPI_remove_promiscuous/Hydroxyurea_Pos_PPI_real.csv")

DMSO_norm_pos = DMSO_matrix[which(DMSO_matrix[,1] %in% DMSO_real[,1]),c(1,3)]
Forskolin_norm_pos = Forskolin_norm[which(Forskolin_norm[,1] %in% Forskolin_real[,1]),c(1,3)]
FK506_norm_pos = FK506_norm[which(FK506_norm[,1] %in% FK506_real[,1]),c(1,3)]
Raffinose_norm_pos = Raffinose_norm[which(Raffinose_norm[,1] %in% Raffinose_real[,1]),c(1,3)]
NaCl_norm_pos = NaCl_norm[which(NaCl_norm[,1] %in% NaCl_real[,1]),c(1,3)]
H2O2_norm_pos = H2O2_norm[which(H2O2_norm[,1] %in% H2O2_real[,1]),c(1,3)]
Dox_norm_pos = Dox_norm[which(Dox_norm[,1] %in% Dox_real[,1]),c(1,3)]
cold_norm_pos = cold_norm[which(cold_norm[,1] %in% cold_real[,1]),c(1,3)]
HU_norm_pos = HU_norm[which(HU_norm[,1] %in% HU_real[,1]),c(1,3)]

all_pos = unique(c(DMSO_norm_pos[,1], Forskolin_norm_pos[,1], FK506_norm_pos[,1], Raffinose_norm_pos[,1], 
                   NaCl_norm_pos[,1], H2O2_norm_pos[,1], Dox_norm_pos[,1], cold_norm_pos[,1], HU_norm_pos[,1])) #14463

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
csvWriter(matrix_pos, "Normalized_multiple_files/Pos_PPI_normalized_fit_primary.csv")

### Put the fitness of PPI of two orientations across different environments onto the same row
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/")
PPI_fit_norm = csvReader_T("Normalized_multiple_files/Pos_PPI_normalized_fit_primary.csv") # 14532
## Make the PPI name for the count and normalized fitness consistent
PPI_count = csvReader_T("Positive_PPI_remove_promiscuous/PPI_environment_count_summary.csv") 
PPI_dup = mark_duplicates_fast(PPI_fit_norm[,1]) # 13829
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
csvWriter(matrix_final, "Normalized_multiple_files/Normalized_fitness_PPI_all_primary.csv")

########### Consider any negative fitness values of any orientation to be zero

PPI_norm = csvReader_T("Normalized_multiple_files/Normalized_fitness_PPI_all_primary.csv")
PPI_count = csvReader_T("Positive_PPI_remove_promiscuous/PPI_environment_count_summary.csv")
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
cor(environment_number, variation_score, method="spearman") # -0.7216613

################################### Filter out PPIs that were called by merge two SD environments
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/")
PPI_heatmap = dataFrameReader_T("Normalized_multiple_files/Variation_score_PPI_environment_neg_zero_SD_merge.csv")
PPI_count = csvReader_T("Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge.csv")
PPI_SD = PPI_count[which(PPI_count[,3] == "1"),]
### Remove single
SD1_single = csvReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO/counts/PPI_1_p.values.csv") # 189787
SD2_single = csvReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_2/counts/PPI_1_p.values.csv") # 168705
common_single = intersect(SD1_single[,1], SD2_single[,1]) # 100412
SD_pos_single = intersect(common_single, PPI_SD[,1])
PPI_heatmap_filter = PPI_heatmap[which(PPI_heatmap[,1] %in% SD_pos_single & as.numeric(PPI_heatmap[,2]) == 1),]
PPI_heatmap_remaining = PPI_heatmap[which(!PPI_heatmap[,1] %in% PPI_heatmap_filter[,1]),]
csvWriter(PPI_heatmap_remaining, "Normalized_multiple_files/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")


#####################################################################################
### Functions with a large number of codes

### Function of calling PPIs
PPI_filtering = function(coeff_matrix, metrics_matrix, p_threshold, p_loc){
  ### Use the coefficients that give the maximum MCC
  best_metrics_threshold = metrics_matrix[which(metrics_matrix$MCC == max(metrics_matrix$MCC)),1]
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
  fragment_protein = unique(as.vector(split_string_vector(fragment_PPI[,1])))
  control_select = rep(0, nrow(pos_PPI_real))
  for(i in 1: nrow(pos_PPI_real)){
    PPI = split_string(pos_PPI_real[i,1])
    if(PPI[1] %in% fragment_protein | PPI[2] %in% fragment_protein){
      control_select[i] = 1
    }
  }
  pos_PPI_real_filter = pos_PPI_real[which(control_select == 0),] 
  csvWriter(pos_PPI_real_filter, "Pos_PPI_real.csv")
  return (best_metrics_threshold) ## return the best threshold 
}


######### Function to make the dynamic threshold
PPI_calling_sigmoid = function(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold){
  reference_name_generation = function(Neg_number_PPI, Neg_ref_number){
    pos_name = "PPI_pos_3_assay"
    neg_name = vector("character", Neg_ref_number)
    for (i in 1:length(neg_name)){
      neg_name[i] = paste("neg", as.character(Neg_number_PPI), as.character(i), sep = "_")
    }
    reference_name = vector("character", length(pos_name)*length(neg_name))
    start_count = 1
    for (i in 1:length(pos_name)){
      for (j in 1:length(neg_name)){
        reference_name[start_count] = paste(pos_name[i], neg_name[j], "reference.csv", sep= "_")
        start_count = start_count + 1
      }
    }
    return(reference_name)
  }
  ### Create reference name
  reference_name = reference_name_generation(Neg_number_PPI, Neg_ref_number)
  
  ### Function to generate a matrix of PPV
  ROC_line_matrix = function(matrix_ref, Fitness, Q_values, output){
    number_Pos = length(which(matrix_ref[,1] == "Positive")) 
    number_Neg = length(which(matrix_ref[,1] == "Negative")) 
    TPR_matrix = matrix(0, length(Fitness), length(Q_values))
    FPR_matrix = matrix(0, length(Fitness), length(Q_values))
    PPV_matrix = matrix(0, length(Fitness), length(Q_values))
    rownames(TPR_matrix) = as.character(Fitness)
    colnames(TPR_matrix) = as.character(Q_values)
    rownames(FPR_matrix) = as.character(Fitness)
    colnames(FPR_matrix) = as.character(Q_values)
    rownames(PPV_matrix) = as.character(Fitness)
    colnames(PPV_matrix) = as.character(Q_values)
    for (i in 1:length(Fitness)){
      for (j in 1:length(Q_values)){
        positive_threshold = matrix_ref[which(matrix_ref[,3] >= Fitness[i] &
                                                log10(matrix_ref[,p_loc]) <= Q_values[j]),1]
        true_positive = length(which(positive_threshold == "Positive"))
        false_positive = length(which(positive_threshold == "Negative"))
        TPR_matrix[i, j] = true_positive/number_Pos
        specificity = (number_Neg - false_positive)/number_Neg
        FPR_matrix[i, j] = 1- specificity
        PPV_matrix[i, j] = true_positive/(true_positive + false_positive)
      }
    }
    if(output == "FPR"){
      return(FPR_matrix)
    }
    else if (output == "TPR"){
      return(TPR_matrix)
    }
    else if (output == "PPV"){
      return(PPV_matrix)
    }
    else{
      return("ERROR")
    }
    
  }
  
  #Find the combinations of Q-value and fitness that give the same FPR, or TPR, or PPV
  Q_fit_PPV_threshold = function(Q_fit_matrix, PPV_threshold){
    Q_fit_matrix_threshold = matrix(NA, ncol(Q_fit_matrix), 3)
    for (i in 1: ncol(Q_fit_matrix)){
      Q_fit_matrix_threshold[i,2] = colnames(Q_fit_matrix)[i]
      for (j in 1: nrow(Q_fit_matrix)){
        if (is.na(Q_fit_matrix[j,i])){
          next
        }
        else {
          if (Q_fit_matrix[j,i] >= PPV_threshold){
            Q_fit_matrix_threshold[i,3] = rownames(Q_fit_matrix)[j]
            break
          }  
        }
      }
    }
    
    Q_fit_matrix_threshold[,1]= paste("PPV",as.character(PPV_threshold), sep = "_")
    colnames(Q_fit_matrix_threshold) = c("PPV", "Q_value", "Fitness")
    return(Q_fit_matrix_threshold)
  }
  ##### Create matrix of PPV with two metrics and find out the combination for each PPV
  for (k in 1: length(reference_name)){
    
    matrix_ref = dataFrameReader_T(reference_name[k])
    PPV_matrix = ROC_line_matrix(matrix_ref, Fitness, p_values, "PPV")
    csvWriter_rownames(PPV_matrix, paste("Fitness_Q_value_PPV", as.character(k), "data.csv", sep = "_"))
    
    PPV_bin = specific_PPV
    matrix_different_PPV = rep(0, 3)
    for (m in 1:length(PPV_bin)){
      matrix_PPV_threshold = Q_fit_PPV_threshold(PPV_matrix, PPV_bin[m])
      matrix_different_PPV = rbind(matrix_different_PPV, matrix_PPV_threshold)
    }
    matrix_different_PPV = matrix_different_PPV[2:nrow(matrix_different_PPV),]
    csvWriter(matrix_different_PPV, paste("Fitness_Q_value_PPV_threshold", as.character(k), "data.csv", sep = "_"))
  }
  
  ######## Derive an equation of fitness and Q-values that can best describe 50 random reference sets
  PPV_threshold_files = rep("Name", 50)
  for (i in 1:50){
    PPV_threshold_files[i] = paste("Fitness_Q_value_PPV_threshold", as.character(i), "data.csv", sep = "_")
  }
  
  library(ggplot2)
  plot_threshold = function(specific_PPV, output){
    PPV_threshold_all = rep(0,3)
    for (k in 1:length(PPV_threshold_files)){
      file = csvReader_T(PPV_threshold_files[k])
      a = file[which(split_string_vector(file[,1])[,2] == as.character(specific_PPV)),]
      for (m in 1:nrow(a)){
        a[m,1] = paste(a[m,1], as.character(k), sep = "_")  
      }
      PPV_threshold_all= rbind(PPV_threshold_all, a)
    }
    PPV_threshold_all = PPV_threshold_all[2:nrow(PPV_threshold_all),]
    csvWriter(PPV_threshold_all, paste("PPV", as.character(specific_PPV), "threshold.csv", sep= "_"))
    
    PPV_threshold_all = csvReader_T(paste("PPV", as.character(specific_PPV), "threshold.csv", sep= "_"))
    fitness_threshold = as.numeric(PPV_threshold_all[,3])
    Q_value_threshold = as.numeric(PPV_threshold_all[,2])
    fitmodel = nls(fitness_threshold  ~ SSlogis( Q_value_threshold, Asym, xmid, scal))
    
    coeff = coef(fitmodel)
    simulated_Q = unique(Q_value_threshold)
    simulated_fit = rep(0, length(simulated_Q))
    for(i in 1:length(simulated_Q)){
      m = simulated_Q[i]
      simulated_fit[i] = coeff[1]/(1+exp((coeff[2] - m)/coeff[3]))
    }
    
    simulated_data = cbind(rep("Fitted line", length(simulated_Q)), rep("Simulation", length(simulated_Q)), simulated_Q, simulated_fit)
    colnames(simulated_data) = c("Data_type", "Lineage", "Q_value", "Fitness")
    csvWriter(simulated_data, paste("PPV", as.character(specific_PPV), "simulated_threshold.csv", sep= "_"))
    
    real_data = cbind(rep("Threshold", nrow(PPV_threshold_all)), PPV_threshold_all)
    colnames(real_data) = c("Data_type", "Lineage", "Q_value", "Fitness")
    csvWriter(real_data, paste("PPV", as.character(specific_PPV), "threshold_plot.csv", sep= "_"))
    
    real_data = dataFrameReader_T(paste("PPV", as.character(specific_PPV), "threshold_plot.csv", sep= "_"))
    simulated_data = dataFrameReader_T(paste("PPV", as.character(specific_PPV), "simulated_threshold.csv", sep= "_"))
    
    ggplot() +
      geom_point(aes(Q_value, Fitness, color = Data_type), real_data, alpha = 0.5)+
      geom_line(aes(Q_value, Fitness, color = Data_type), simulated_data) +
      scale_color_manual('', breaks = c("Threshold", "Fitted line"),
                         values = apple_colors[c(7,8)], 
                         guide = guide_legend(override.aes = list(linetype = c(NA,1),
                                                                  shape = c(16, NA)))) +
      scale_x_continuous(name= "Log10(P-value)",
                         limits = c(p_threshold, 0),
                         breaks = seq(p_threshold, 0, by = 0.2),
                         labels = seq(p_threshold, 0, by = 0.2)) +
      #scale_y_continuous(name = "Fitness",
      #limits = c(0, 0.6),
      #breaks = seq(0, 0.6, by = 0.1),
      #labels = seq(0, 0.6, by = 0.1)) +
      
      theme(legend.key=element_blank()) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text.x = element_text(size = 10, color = "black"), 
            axis.text.y.left = element_text(size = 10, color = "black")) + 
      theme(text = element_text(size=12))
    
    ggsave(output, width = 7, height =5, device = "pdf")
    
    
  }
  ###### Output the pdf plot for the fitting line
  output = rep("A", length(specific_PPV))
  for(i in 1:length(specific_PPV)){
    output[i] = paste("PPV", as.character(specific_PPV[i]), "threshold.pdf", sep= "_")
  }
  
  for(i in 1:length(specific_PPV)){
    plot_threshold(specific_PPV[i], output[i])
  }
  
  ####### Output the coefficient matrix
  PPV_file_name = rep("A", length(specific_PPV))
  for(i in 1:length(specific_PPV)){
    PPV_file_name[i] = paste("PPV", as.character(specific_PPV[i]), "threshold.csv", sep= "_")
  }
  coeff_calling = function(PPV_file_name, PPV){
    coeff_matrix = matrix(0, length(PPV_file_name), 4)
    coeff_matrix[,1] = PPV
    for (i in 1:length(PPV_file_name)){
      file_name= PPV_file_name[i]
      PPV_data = dataFrameReader_T(file_name)
      fitness_threshold = PPV_data[,3]
      Q_value_threshold = PPV_data[,2]
      fitmodel_2 = nls(fitness_threshold ~ SSlogis(Q_value_threshold, Asym, xmid, scal))
      coeff_matrix[i, 2: ncol(coeff_matrix)] = coef(fitmodel_2)
    }
    colnames(coeff_matrix) = c("threshold", "asym", "xmid", "scal")
    csvWriter(coeff_matrix, "coefficient_matrix.csv")
  }
  ## write coefficients into a matrix
  coeff_calling(PPV_file_name, specific_PPV)
  coeff_matrix = csvReader_T("coefficient_matrix.csv")
  
  #### Use the coefficients to describe the relationship between fitness and Q_value
  #### Output the number of PPIs, overlap with PCA, BioGRID, and PRS, RRS
  PPI_multiple_pos_fit = PPI_multiple[which(PPI_multiple[,3] > 0),] 
  fitness_multiple = PPI_multiple_pos_fit[,3]
  Q_value_multiple = PPI_multiple_pos_fit[,p_loc]
  log10_Q_value = log10(Q_value_multiple)
  log10_Q_value[log10_Q_value == -Inf] = -350
  index_1 = which(log10_Q_value >= p_threshold)
  index_2 = which(log10_Q_value < p_threshold)
  fitness_multiple_1 = fitness_multiple[index_1]#612685
  fitness_multiple_2 = fitness_multiple[index_2] #5584
  PPI_paper= read.delim(file = "~/Desktop/Big_PPiSeq_paper/PCA_rearray/PPI_set_science.txt", sep=" ")
  PPI_paper_Michnick= paste(PPI_paper[,1], PPI_paper[,3], sep="_") # 2770
  
  PPI_summary = csvReader_T("/Volumes/zmliu_02/PPiseq/H2O2/reference_set/one_binary/multiple_validated_PPI.csv")
  PPI_BIOGRID = PPI_summary[,1]
  
  count_summary = matrix(0, nrow(coeff_matrix), 10)
  count_summary[,1] = coeff_matrix[,1]
  for(i in 1:nrow(coeff_matrix)){
    coeff_select = as.numeric(coeff_matrix[i, 2:ncol(coeff_matrix)])
    pos_PPI = rep(0, nrow(PPI_multiple_pos_fit))
    a = log10_Q_value[index_1]
    Q_matrix_cal = coeff_select[1]/(1 + exp((coeff_select[2] - a)/coeff_select[3]))
    pos_PPI[index_1] = (fitness_multiple_1 >= Q_matrix_cal)
    count_summary[i,2] = length(which(pos_PPI == 1))
    pos_PPI_non_significant = PPI_multiple_pos_fit[which(pos_PPI == 1), 1:2]
    count_summary[i,3] = nrow(match_both_direction(pos_PPI_non_significant, PPI_paper_Michnick))
    count_summary[i,4] = nrow(match_both_direction(pos_PPI_non_significant, PPI_BIOGRID))
    pos_PPI_min = min(PPI_multiple_pos_fit[which(pos_PPI == 1), 3])
    count_summary[i,5] = pos_PPI_min
    pos_PPI[index_2] = (fitness_multiple_2 >= pos_PPI_min)
    count_summary[i,6] = length(which(pos_PPI == 1))
    PPI_called = PPI_multiple_pos_fit[which(pos_PPI == 1),1:2]
    count_summary[i,7] = nrow(match_both_direction(PPI_called, PPI_paper_Michnick))
    count_summary[i,8] = nrow(match_both_direction(PPI_called, PPI_BIOGRID))
    PRS_PPI = PPI_called[grep("Pos", PPI_called[,1]),1]
    RRS_PPI = PPI_called[grep("Neg", PPI_called[,1]),1]
    if(length(PRS_PPI) > 0){
      count_summary[i,9] = length(unique(split_string_vector_line(PRS_PPI)[,1]))  
    }else{count_summary[i,9] = 0}
    
    if(length(RRS_PPI) > 0){
      count_summary[i,10] = length(unique(split_string_vector_line(RRS_PPI)[,1]))
    }else{
      count_summary[i,10] = 0
    }
    
  }
  colnames(count_summary) = c("threshold", "larger_1.3","large_1.3_overlap_PCA", "large_1.3_overlap_BioGRID",
                              "fitness_threshold","all", "Overlap_Michnick",
                              "Overlap_BioGRID", "PRS_count", "RRS_count")
  
  csvWriter(count_summary, "Positive_PPI_number_threshold_1.csv")
  
  ##### Calculate the FPR, TPR, real_PPV under different PPV thresholds
  
  rate_calculate_PPV = function(matrix_ref, coeff, fpr){ # fpr =1: calculate FPR for coefficents = TPR lines
    fitness = matrix_ref[,3]
    Q_value = log10(matrix_ref[,p_loc])
    Q_value[which(Q_value == -Inf)] = -350
    index_1 = which((fitness >= coeff[1]/(1+ exp((coeff[2] - Q_value)/coeff[3]))) & Q_value >= p_threshold)
    minimum_fitness = min(fitness[index_1])
    index_2 = which(fitness >= minimum_fitness & Q_value < p_threshold)
    index = c(index_1, index_2)
    Pos_PPI = matrix_ref[index,1]
    number_Pos = length(which(matrix_ref[,1] == "Positive")) 
    number_Neg = length(which(matrix_ref[,1] == "Negative"))
    true_pos = length(which(Pos_PPI == "Positive"))
    false_pos = length(which(Pos_PPI == "Negative"))
    if (fpr == 1){
      FPR = false_pos/number_Neg
      return(FPR)
    }else{
      FPR = false_pos/number_Neg
      TPR = true_pos/number_Pos
      PPV = true_pos/(true_pos + false_pos)
      return (c(FPR, TPR, PPV, true_pos, false_pos))
    }
  } 
  # Calculate Mathew correlaltion coefficient which combine precision and recall
  MCC_calculate = function(matrix_ref, coeff){ 
    library(mccr)
    fitness = matrix_ref[,3]
    Q_value = log10(matrix_ref[,p_loc])
    Q_value[which(Q_value == -Inf)] = -350
    index_1 = which((fitness >= coeff[1]/(1+ exp((coeff[2] - Q_value)/coeff[3]))) & Q_value >= p_threshold)
    if (length(index_1) >0){
      minimum_fitness = min(fitness[index_1])
      index_2 = which(fitness >= minimum_fitness & Q_value < p_threshold)
      index = c(index_1, index_2)
      predict_value = rep(0, nrow(matrix_ref))
      predict_value[index] = 1
    }
    else{
      predict_value = rep(0, nrow(matrix_ref))
    }
    
    number_Pos = length(which(matrix_ref[,1] == "Positive")) 
    number_Neg = length(which(matrix_ref[,1] == "Negative"))
    act_value = c(rep(1, number_Pos), rep(0, number_Neg))
    mcc_value = mccr(act_value, predict_value)
    return(mcc_value)
  } 
  
  # output a matrix that contain various metrics measured with reference sets
  matrix_final = rep(0,11)
  for (k in 1:length(reference_name)){
    matrix_ref = dataFrameReader_T(reference_name[k])
    matrix_different_PPV = read.csv(paste("Fitness_Q_value_PPV_threshold", as.character(k), "data.csv", sep = "_"))
    PPV_unique = unique(matrix_different_PPV[,1])
    matrix_fit_q_value = matrix(0,length(PPV_unique), 11)
    matrix_fit_q_value[,1] = rep(paste("data", as.character(k), sep = "_"), length(PPV_unique))
    matrix_fit_q_value[,2] = split_string_vector(PPV_unique)[,2]
    for (j in 1:length(PPV_unique)){
      PPV_chosen = as.numeric(matrix_fit_q_value[j,2])
      fit = as.numeric(coeff_matrix[which(as.numeric(coeff_matrix[,1]) == PPV_chosen), 2:4]) 
      ## use the fitted parameter generated from all 50 reference sets
      matrix_fit_q_value[j,3:5] = fit
      c = rate_calculate_PPV(matrix_ref, fit, 0)
      matrix_fit_q_value[j,6] = c[1]
      matrix_fit_q_value[j,7] = c[2]
      matrix_fit_q_value[j,8] = c[3]
      matrix_fit_q_value[j,9] = c[4]
      matrix_fit_q_value[j,10] = c[5]
      matrix_fit_q_value[j,11] = MCC_calculate(matrix_ref, fit)  
    }
    matrix_fit_q_value = na.omit(matrix_fit_q_value)
    matrix_final = rbind(matrix_final, matrix_fit_q_value)
  }
  matrix_final = matrix_final[2:nrow(matrix_final),]
  colnames(matrix_final) = c("Reference","PPV", "asym", "xmid", "scal", "FPR",
                             "TPR", "Real_PPV", "True_pos", "False_pos", "MCC")
  csvWriter(matrix_final, "PPV_threshold_all_data.csv")
  
  
  ###### Take the mean TPR under different specific FPRs so that we can choose an optimum combination of FPR and TPR
  matrix_final = dataFrameReader_T("PPV_threshold_all_data.csv")
  PPV_unique = unique(matrix_final$PPV)
  FPR_mean = rep(0, length(PPV_unique))
  TPR_mean = rep(0, length(PPV_unique))
  Real_PPV_mean = rep(0, length(PPV_unique))
  true_pos_mean = rep(0, length(PPV_unique))
  false_pos_mean = rep(0, length(PPV_unique))
  MCC_mean = rep(0, length(PPV_unique))
  for(i in 1:length(PPV_unique)){
    a = which(matrix_final$PPV == PPV_unique[i])
    FPR_mean[i] = mean(matrix_final[a,6])
    TPR_mean[i] = mean(matrix_final[a,7])
    Real_PPV_mean[i] = mean(matrix_final[a,8])
    true_pos_mean[i] = mean(matrix_final[a,9])
    false_pos_mean[i] = mean(matrix_final[a,10])
    MCC_mean[i] = mean(matrix_final[a,11])
  }
  matrix_average = cbind(PPV_unique, FPR_mean, TPR_mean, Real_PPV_mean, true_pos_mean, false_pos_mean, MCC_mean)
  colnames(matrix_average) = c("Combined_PPV", "FPR", "TPR", "PPV", "True_pos", "False_pos", "MCC")
  csvWriter(matrix_average, "Mean_FPR_TPR_PPV_1.csv")
  
  mean_FPR_TPR_PPV = dataFrameReader_T("Mean_FPR_TPR_PPV_1.csv")
  number_pos_PPI = dataFrameReader_T("Positive_PPI_number_threshold_1.csv")
  specificity = 1- mean_FPR_TPR_PPV$FPR
  balanced_accuracy = (mean_FPR_TPR_PPV$TPR + specificity)/2
  F1_score = rep(0, nrow(mean_FPR_TPR_PPV))
  for(i in 1:length(F1_score)){
    a = mean_FPR_TPR_PPV[i, 3:4]
    F1_score[i] = 1/mean(as.numeric(1/a))
  }
  matrix_threshold_final = cbind(mean_FPR_TPR_PPV, balanced_accuracy, F1_score, number_pos_PPI[,2:ncol(number_pos_PPI)])
  
  csvWriter(matrix_threshold_final, "Threshold_metrics_final.csv")
  specificity = 1- mean_FPR_TPR_PPV$FPR
  balanced_accuracy = (mean_FPR_TPR_PPV$TPR + specificity)/2
  F1_score = rep(0, nrow(mean_FPR_TPR_PPV))
  for(i in 1:length(F1_score)){
    a = mean_FPR_TPR_PPV[i, 3:4]
    F1_score[i] = 1/mean(as.numeric(1/a))
  }
  matrix_threshold_final = cbind(mean_FPR_TPR_PPV, balanced_accuracy, F1_score, number_pos_PPI[,2:ncol(number_pos_PPI)])
  
  csvWriter(matrix_threshold_final, "Threshold_metrics_final.csv")
  
}


