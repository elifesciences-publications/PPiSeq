#This script will combine the lineage counts and their fitness and filter out bad fitness measurements
source("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
library(ggplot2)
library(scales)

distance_calculate = function(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment){
        real_lineage_counts = real_lineage[,4:(4 + time_point -1)]
        real_lineage_counts[real_lineage_counts == 0] = 1
        distance = (rowSums(((predict_lineage - real_lineage_counts)/real_lineage_counts)^2))^0.5
        real_lineage_counts = real_lineage[,4:(4 + time_point -1)]
        lineage_sum = rowSums(real_lineage_counts) 
        data_threshold = data.frame(lineage_sum, distance)
        ### Make a hexplot to show distance versus lineage_sum
        ggplot() +
          geom_hex(aes(x= lineage_sum, y= distance, fill = log10(..count..)), data_threshold, bins = 60)+
          geom_hline(yintercept=distance_threshold, linetype="dashed", color = apple_colors[5])+
          scale_fill_gradientn(colours = apple_colors[c(8, 3, 7)])+
          #geom_hline(yintercept = 26, color = apple_colors[5])+
          #scale_y_continuous(name = expression(italic(d)),
          scale_y_continuous(name = expression(italic(d)),
                             limits=c(1, 2^8),
                             trans = log2_trans(),
                             breaks = trans_breaks("log2", function(x) 2^x),
                             labels = trans_format("log2", math_format(2^.x))) +
          scale_x_continuous(name = "Sum of lineage counts",
                             limits = c(9, 1e4),
                             trans = log10_trans(),
                             breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x)))+
          labs(fill = expression('Log'[10]* '(count)')) +     
          theme(legend.position =c(0.9,0.8), legend.key=element_blank(), legend.text=element_text(size=10)) +
          #guides(fill=guide_legend(title="Log10(Count)")) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          theme(axis.text.x = element_text(size = 10, color = "black"),
                axis.text.y.left = element_text(size = 10, color = "black"))
        
        ggsave("distance_threshold/distance_sum_count.pdf", height =5, width =5)
        
        ### Make a histogram of distance
        library(ggplot2)
        library(scales)
        ggplot() + aes(distance) +
          geom_histogram(breaks = seq(0,380,by=1), color = apple_colors[11], fill = apple_colors[3]) +
          geom_vline(xintercept=distance_threshold, linetype="dashed", color = apple_colors[5])+
          scale_y_continuous(name = "Frequency",
                             trans = log10_trans(),
                             breaks = trans_breaks("log10", function(x) 10^x),
                             labels = trans_format("log10", math_format(10^.x)))+
          scale_x_continuous(name = "Distance",
                             limits=c(0, 100),
                             breaks = seq(0,100, by =5),
                             labels = seq(0,100, by= 5)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          theme(axis.text.x = element_text(size = 10, color = "black"),
                axis.text.y.left = element_text(size = 10, color = "black"))
        ggsave("distance_threshold/Histogram_distance.pdf", width = 6, height = 5)
        
        ### filter out lineages that have bad fitness measurement
        real_lineage_filter = real_lineage[which(distance < distance_threshold),]
        bad_lineage_count = length(which(distance >= distance_threshold)) 
        
        fitness_filter = fitness[which(distance < distance_threshold)]
        distance_filter = distance[which(distance < distance_threshold)]
        known_PPI_filtered = as.matrix(cbind(real_lineage_filter[,1:3], rep(environment, length(fitness_filter)), 
                                          fitness_filter, distance_filter, real_lineage_filter[,4:(4 + time_point -1)]))
        colnames(known_PPI_filtered) = c(colnames(real_lineage_filter)[1:3], "Environment", 
                                                "Fitness", "Distance", colnames(real_lineage_filter)[4:(4 + time_point -1)])
        
        ### Split PRS and RRS into two halves (Consider the same PPI with different directions as different PPIs)
        PRS_PPI = csvReader_T("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/PPI_barcodes/Pos_PPI_name_split.csv")
        RRS_PPI = csvReader_T("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/PPI_barcodes/Neg_PPI_name_split.csv")
        PRS_RRS = rbind(PRS_PPI, RRS_PPI) 
        known_PPI_filtered_PRS_RRS = known_PPI_filtered[which(known_PPI_filtered[,3] %in% PRS_RRS[,2]),] 
        known_PPI_filtered_remaining = known_PPI_filtered[which(!known_PPI_filtered[,3] %in% PRS_RRS[,2]),] 
        for(i in 1:nrow(PRS_RRS)){
          index = which(known_PPI_filtered_PRS_RRS[,3] == PRS_RRS[i,2])
          known_PPI_filtered_PRS_RRS[index,1] = PRS_RRS[i,3]
          known_PPI_filtered_PRS_RRS[index,2] = PRS_RRS[i,3]
        }
        
        known_PPI_filtered_rename = rbind(known_PPI_filtered_PRS_RRS, known_PPI_filtered_remaining)
        csvWriter(known_PPI_filtered_rename, "known_PPI_fitness_barcodes_filtered_rename.csv")
        return (bad_lineage_count)
}

### DMSO
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 5
fit_seq = dataFrameReader_T("SD1_FitSeq.csv")
predict_lineage = dataFrameReader_F("SD1_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "DMSO"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment) #336


### DMSO_2
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_2/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 5
fit_seq = dataFrameReader_T("SD2_FitSeq.csv")
predict_lineage = dataFrameReader_F("SD2_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "DMSO2"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment) # 620

### H2O2
setwd("/Volumes/zmliu_02/PPiseq_03/H2O2/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 5
fit_seq = dataFrameReader_T("H2O2_FitSeq.csv")
predict_lineage = dataFrameReader_F("H2O2_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "H2O2"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment) #518

### HU
setwd("/Volumes/zmliu_02/PPiseq_03/HU/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 4
fit_seq = dataFrameReader_T("HU_FitSeq.csv")
predict_lineage = dataFrameReader_F("HU_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "HU"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment) # 289

### Dox
setwd("/Volumes/zmliu_02/PPiseq_03/Dox/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 4
fit_seq = dataFrameReader_T("Dox_FitSeq.csv")
predict_lineage = dataFrameReader_F("Dox_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "Doxorubicin"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment) # 448

### FK506
setwd("/Volumes/zmliu_02/PPiseq_03/FK506/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 4
fit_seq = dataFrameReader_T("FK506_FitSeq.csv")
predict_lineage = dataFrameReader_F("FK506_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "FK506"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment) # 396

### Raffinose
setwd("/Volumes/zmliu_02/PPiseq_03/Raffinose/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 4
fit_seq = dataFrameReader_T("Raffinose_FitSeq.csv")
predict_lineage = dataFrameReader_F("Raffinose_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "Raffinose"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment) # 346


### Forskolin
setwd("/Volumes/zmliu_02/PPiseq_03/Forskolin/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 4
fit_seq = dataFrameReader_T("Forskolin_FitSeq.csv")
predict_lineage = dataFrameReader_F("Forskolin_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "Forskolin"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment)#198


### NaCl
setwd("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 4
fit_seq = dataFrameReader_T("NaCl_FitSeq.csv")
predict_lineage = dataFrameReader_F("NaCl_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "NaCl"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment)#1312

### 16C
setwd("/Volumes/zmliu_02/PPiseq_03/16C/counts/")
real_lineage = dataFrameReader_T("Good_lineage_trajectories_add_sum_new.csv")
time_point = 4
fit_seq = dataFrameReader_T("16C_FitSeq.csv")
predict_lineage = dataFrameReader_F("16C_predict.csv")
fitness = fit_seq$Estimated_Fitness
distance_threshold = 19
environment = "16C"
distance_calculate(real_lineage, time_point, predict_lineage, fitness, distance_threshold, environment) #2232

