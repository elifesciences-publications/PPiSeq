# This scirpt is to calculate a p-value for each PPI compare to 
# a negative control group (fragment strains (ORF-DHFR[1,2] or DHFR[3] X ORF)
source("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
library(ggplot2)
library(scales)

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

## DMSO
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"

p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## DMSO_2
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_2/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"

p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## H2O2
setwd("/Volumes/zmliu_02/PPiseq_03/H2O2/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## HU
setwd("/Volumes/zmliu_02/PPiseq_03/HU/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## Forskolin
setwd("/Volumes/zmliu_02/PPiseq_03/Forskolin/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## Doxorubicin
setwd("/Volumes/zmliu_02/PPiseq_03/Dox/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## FK506
setwd("/Volumes/zmliu_02/PPiseq_03/FK506/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## Raffinose
setwd("/Volumes/zmliu_02/PPiseq_03/Raffinose/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## NaCl
setwd("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

## 16C
setwd("/Volumes/zmliu_02/PPiseq_03/16C/counts/")
PPI_fitness= "known_PPI_fitness_barcodes_sorted.csv"
P_value_output_01= "PPI_1_p.values.csv"
P_value_output_02= "PPI_multiple_p.values.csv"
p_value_calculate(PPI_fitness, P_value_output_01, P_value_output_02)

