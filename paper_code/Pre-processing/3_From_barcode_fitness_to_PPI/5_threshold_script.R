source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

###############################################################################


# DMSO
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.9, by = 0.02))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

#DMSO_02
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_2/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.9, by = 0.02))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_2/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

#H2O2
setwd("/Volumes/zmliu_02/PPiseq_03/H2O2/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.9, by = 0.02))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/H2O2/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

#HU
setwd("/Volumes/zmliu_02/PPiseq_03/HU/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.9, by = 0.02))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/HU/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

#Forskolin
setwd("/Volumes/zmliu_02/PPiseq_03/Forskolin/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.9, by = 0.02))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Forskolin/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

#Dox
setwd("/Volumes/zmliu_02/PPiseq_03/Dox/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.9, by = 0.02))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Dox/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

#NaCl
setwd("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.9, by = 0.02))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

# 16C

setwd("/Volumes/zmliu_02/PPiseq_03/16C/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.2, 0.4, by= 0.01), seq(0.41, 0.6, by= 0.01))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/16C/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

# FK506
setwd("/Volumes/zmliu_02/PPiseq_03/FK506/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.84, by = 0.02))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/FK506/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)

# Raffinose
setwd("/Volumes/zmliu_02/PPiseq_03/Raffinose/reference_set/p_value/")
Neg_number_PPI = 4.5e4
Neg_ref_number = 50
p_values = seq(-4, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)
p_threshold = -4
p_loc = 6 # the column in the reference set as the threshold
specific_PPV = c(seq(0.3, 0.4, by= 0.02), seq(0.41, 0.58, by= 0.01))
PPI_multiple = dataFrameReader_T("/Volumes/zmliu_02/PPiseq_03/Raffinose/counts/PPI_multiple_p.values.csv")
PPI_calling_sigmoid(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold)


####################################################################################
##### sigmoid function to call PPIs

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
  PPI_paper= read.delim(file = "~/Dropbox/PPiSeq_02/Paper_data/Outside_datasets/PPI_set_PCA_science.txt", sep=" ")
  PPI_paper_Michnick= paste(PPI_paper[,1], PPI_paper[,3], sep="_") # 2770
  
  PPI_summary = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/Useful_datasets/multiple_validated_PPI.csv")
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

####################### Polynomial regerssion to fit the line
PPI_calling = function(Neg_number_PPI, Neg_ref_number, p_value, fitness, p_loc, specific_PPV, p_threshold){
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
                # Fit the data with bionomial 
                fitness_threshold = as.numeric(PPV_threshold_all[,3])
                Q_value_threshold = as.numeric(PPV_threshold_all[,2])
                coef= lm(fitness_threshold~poly(Q_value_threshold,2, raw = TRUE))$coefficients
                simulated_Q = unique(Q_value_threshold)
                simulated_fit = rep(0, length(simulated_Q))
                for(i in 1:length(simulated_Q)){
                        m = simulated_Q[i]
                        n = c(1, m, m^2)
                        simulated_fit[i] = sum(n * coef)
                }
                
                simulated_data = cbind(rep("Fitted line", length(simulated_Q)), rep("Simulation", length(simulated_Q)), simulated_Q, simulated_fit)
                colnames(simulated_data) = c("Data_type", "Lineage", "Q_value", "Fitness")
                csvWriter(simulated_data, paste("PPV", as.character(specific_PPV), "simulated_threshold.csv", sep= "_"))
                
                real_data = cbind(rep("Threshold", nrow(PPV_threshold_all)), PPV_threshold_all)
                colnames(real_data) = c("Data_type", "Lineage", "Q_value", "Fitness")
                csvWriter(real_data, paste("PPV", as.character(specific_PPV), "threshold_plot.csv", sep= "_"))
                
                real_data = dataFrameReader_T(paste("PPV", as.character(specific_PPV), "threshold_plot.csv", sep= "_"))
                simulated_data = dataFrameReader_T(paste("PPV", as.character(specific_PPV), "simulated_threshold.csv", sep= "_"))
                #######
                
                ggplot() +
                        geom_point(aes(Q_value, Fitness, color = Data_type), real_data, alpha = 0.5)+
                        geom_line(aes(Q_value, Fitness, color = Data_type), simulated_data) +
                        scale_color_manual('', breaks = c("Threshold", "Fitted line"),
                                           values = apple_colors[c(7,8)], 
                                           guide = guide_legend(override.aes = list(linetype = c(NA,1),
                                                                                    shape = c(16, NA)))) +
                        scale_x_continuous(name= "Log10(Q-value)",
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
                        coeff_matrix[i, 2: ncol(coeff_matrix)] = lm(fitness_threshold~poly(Q_value_threshold,2, raw = TRUE))$coefficients
                }
                colnames(coeff_matrix) = c("threshold", "Intercept", "X", "X2")
                csvWriter(coeff_matrix, "coefficient_matrix.csv")
                return(coeff_matrix)
        }
        ## write coefficients into a matrix
        coeff_matrix = coeff_calling(PPV_file_name, specific_PPV)
        
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
        #for(i in 1:10){
                #i = 10
                coeff_select = coeff_matrix[i, 2:ncol(coeff_matrix)]
                pos_PPI = rep(0, nrow(PPI_multiple_pos_fit))
                a = log10_Q_value[index_1]
                b = a^2
                d = rep(1, length(a))
                Q_matrix = cbind.data.frame(d, a, b)
                Q_matrix_cal = as.matrix(Q_matrix) %*% coeff_select
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
                index_1 = which((fitness >= (coeff[1] + coeff[2] * Q_value + coeff[3] * Q_value^2)) & Q_value >= p_threshold)
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
                index_1 = which((fitness >= (coeff[1] + coeff[2] * Q_value + coeff[3] * Q_value^2)) & Q_value >= p_threshold)
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
                        index = which(matrix_different_PPV[,1] == PPV_unique[j])
                        a = matrix_different_PPV[index,2]
                        b = matrix_different_PPV[index,3]
                        matrix_a_b = cbind(a,b)
                        matrix_a_b = na.omit(matrix_a_b)
                        if (nrow(matrix_a_b) > 3){
                                fit= lm(as.numeric(matrix_a_b[,2])~poly(as.numeric(matrix_a_b[,1]),2,raw=TRUE))$coefficients
                                matrix_fit_q_value[j,3:5] = fit
                                c = rate_calculate_PPV(matrix_ref, fit, 0)
                                matrix_fit_q_value[j,6] = c[1]
                                matrix_fit_q_value[j,7] = c[2]
                                matrix_fit_q_value[j,8] = c[3]
                                matrix_fit_q_value[j,9] = c[4]
                                matrix_fit_q_value[j,10] = c[5]
                                matrix_fit_q_value[j,11] = MCC_calculate(matrix_ref, fit)  
                        }else {matrix_fit_q_value[j, 3:11] = NA}
                        
                }
                matrix_fit_q_value = na.omit(matrix_fit_q_value)
                matrix_final = rbind(matrix_final, matrix_fit_q_value)
        }
        matrix_final = matrix_final[2:nrow(matrix_final),]
        colnames(matrix_final) = c("Reference","PPV", "Intercetp", "X", "X2", "FPR",
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


