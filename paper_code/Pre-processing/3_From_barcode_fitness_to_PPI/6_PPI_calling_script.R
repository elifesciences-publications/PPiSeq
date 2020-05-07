source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

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

### DMSO
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.7

### DMSO_2
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_2/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_2/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_2/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.73

### H2O2
setwd("/Volumes/zmliu_02/PPiseq_03/H2O2/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/H2O2/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/H2O2/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.73

### HU
setwd("/Volumes/zmliu_02/PPiseq_03/HU/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/HU/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/HU/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.74

### Forskolin
setwd("/Volumes/zmliu_02/PPiseq_03/Forskolin/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Forskolin/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Forskolin/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.64

### Dox
setwd("/Volumes/zmliu_02/PPiseq_03/Dox/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Dox/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Dox/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.77

### NaCl
setwd("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.73

### 16C
setwd("/Volumes/zmliu_02/PPiseq_03/16C/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/16C/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/16C/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.41

### Raffinose
setwd("/Volumes/zmliu_02/PPiseq_03/Raffinose/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Raffinose/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Raffinose/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.48

### FK506
setwd("/Volumes/zmliu_02/PPiseq_03/FK506/counts/")
PPI_multiple = dataFrameReader_T("PPI_multiple_p.values.csv")
coeff_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/FK506/reference_set/p_value/coefficient_matrix.csv")
metrics_matrix = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/FK506/reference_set/p_value/Threshold_metrics_final.csv")
p_threshold = -4
p_loc = 6
PPI_filtering(coeff_matrix, metrics_matrix, p_threshold, p_loc) # 0.72