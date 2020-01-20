source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
## (2-2-2) Make a figure to show the above method is better than arbitrarily choose a specific Q_value and fitness
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/reference_set/p_value/")
ROC_matrix = function(matrix_ref, selected_fitness, selected_Q_value, pos_neg_ratio, index, Fitness, Q_values){
        number_Pos = length(which(matrix_ref[,1] == "Positive")) 
        number_Neg = length(which(matrix_ref[,1] == "Negative")) 
        l_f = length(Fitness)
        l_q = length(Q_values)
        row_number = l_f * l_q
        ROC_matrix = matrix(0, row_number, 7)
        colnames(ROC_matrix) = c("Index", "Ratio", "Fitness", "Q_value", "TPR", "FPR", "PPV")
        if (selected_fitness != "empty" & selected_Q_value == "empty"){
                for(i in 1:l_f){
                        ROC_matrix[(i *l_q - (l_q-1)): (i*l_q),3] = rep(Fitness[i], l_q)
                        ROC_matrix[(i *l_q - (l_q-1)): (i*l_q),4] = Q_values
                }
        }else if (selected_fitness == "empty" & selected_Q_value != "empty"){
                for(i in 1:l_q){
                        ROC_matrix[(i *l_f - (l_f-1)): (i*l_f),3] = Fitness
                        ROC_matrix[(i *l_f - (l_f-1)): (i*l_f),4] = rep(Q_values[i], l_f)
                }
        }
        
        for(i in 1:nrow(ROC_matrix)){
                positive_threshold = matrix_ref[which(matrix_ref[,3] >= ROC_matrix[i,3] &
                                                              log10(matrix_ref[,6]) <= ROC_matrix[i,4]),1]
                true_positive = length(which(positive_threshold == "Positive"))
                false_positive = length(which(positive_threshold == "Negative"))
                ROC_matrix[i, 5] = true_positive/number_Pos
                specificity = (number_Neg - false_positive)/number_Neg
                ROC_matrix[i, 6] = 1 - specificity
                ROC_matrix[i, 7] = true_positive/(true_positive + false_positive)
                
        }
        if (selected_fitness != "empty" & selected_Q_value == "empty"){
                ROC_matrix_selected = ROC_matrix[which(ROC_matrix[,3] == selected_fitness),]
        }else if (selected_fitness == "empty" & selected_Q_value != "empty"){
                ROC_matrix_selected = ROC_matrix[which(ROC_matrix[,4] == selected_Q_value),]
        }
        
        ROC_matrix_selected[,1] = index
        ROC_matrix_selected[,2] = pos_neg_ratio
        return (ROC_matrix_selected)
}

number_PPI = 6e4
data_set_number = 50
pos_name = "PPI_pos_3_assay"
neg_name = vector("character", data_set_number)
for (i in 1:length(neg_name)){
        neg_name[i] = paste("neg", as.character(number_PPI), as.character(i), sep = "_")
}

reference_name = vector("character", length(pos_name)*length(neg_name))
length(reference_name) #50
start_count = 1
for (i in 1:length(pos_name)){
        for (j in 1:length(neg_name)){
                reference_name[start_count] = paste(pos_name[i], neg_name[j], "reference.csv", sep= "_")
                start_count = start_count + 1
        }
}

# Combine 50 datasets into one matrix
matrix_fit_final = rep(0,7)
matrix_Q_value_final = rep(0,7)
for(k in 1:length(reference_name)){
        matrix_ref = dataFrameReader_T(reference_name[k])
        Q_values = seq(-4, 0, by = 0.1)
        Fitness = seq(0, 1, by = 0.05)
        Q_values_select = seq(-4, 0, by = 0.1)
        Fitness_select = seq(0, 0.5, by = 0.05)
        matrix_fit = rep(0,7)
        for (i in 1:length(Fitness_select)){
                pos_count = length(which(matrix_ref[,1] == "Positive"))
                neg_count = length(which(matrix_ref[,1] == "Negative"))
                pos_neg_ratio = pos_count/neg_count
                matrix_ref_select = ROC_matrix(matrix_ref, Fitness_select[i], "empty", pos_neg_ratio, 
                                               paste(as.character(k), as.character(i), sep= "_"),
                                               Fitness, Q_values)
                matrix_fit = rbind(matrix_fit, matrix_ref_select)
        }
        matrix_fit = matrix_fit[2:nrow(matrix_fit),]
        matrix_fit_final = rbind(matrix_fit_final, matrix_fit)
        
        matrix_Q_value = rep(0,7)
        for (i in 1:length(Q_values_select)){
                pos_count = length(which(matrix_ref[,1] == "Positive"))
                neg_count = length(which(matrix_ref[,1] == "Negative"))
                pos_neg_ratio = pos_count/neg_count
                matrix_ref_select = ROC_matrix(matrix_ref, "empty", Q_values_select[i], pos_neg_ratio, 
                                               paste(as.character(k), as.character(i), sep= "_"),
                                               Fitness, Q_values)
                matrix_Q_value = rbind(matrix_Q_value, matrix_ref_select)
        }
        matrix_Q_value = matrix_Q_value[2:nrow(matrix_Q_value),]
        matrix_Q_value_final = rbind(matrix_Q_value_final, matrix_Q_value)
}
matrix_fit_final = matrix_fit_final[2:nrow(matrix_fit_final),]
matrix_Q_value_final = matrix_Q_value_final[2:nrow(matrix_Q_value_final),]

csvWriter(matrix_fit_final, "Specific_fitness_different_Q_value_ROC_curve.csv")
csvWriter(matrix_Q_value_final, "Specific_Q_value_different_fitness_ROC_curve.csv")

# Scatter plot to show the TPR value of 50 datasets with FPR = 5e-4 with three data sets
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/reference_set/p_value/")
Con_fit_var_Q = dataFrameReader_T("Specific_fitness_different_Q_value_ROC_curve.csv")
Con_Q_var_fit = dataFrameReader_T("Specific_Q_value_different_fitness_ROC_curve.csv")
matrix_FPR_TPR = dataFrameReader_T("PPV_threshold_all_data.csv")
metrics_matrix = dataFrameReader_T("Threshold_metrics_final.csv")
best_metrics_threshold = metrics_matrix[which(metrics_matrix$MCC == max(metrics_matrix$MCC)),1] # 0.71
matrix_FPR_TPR_chosen = matrix_FPR_TPR[which(matrix_FPR_TPR$PPV == best_metrics_threshold),] # 50
min_FPR = min(matrix_FPR_TPR_chosen$FPR) # 0.001832401
max_FPR = max(matrix_FPR_TPR_chosen$FPR) # 0.003860438

matrix_FPR_TPR_select = matrix_FPR_TPR[which(matrix_FPR_TPR$FPR >= min_FPR & matrix_FPR_TPR$FPR <= max_FPR),] # 266
Con_fit_var_Q_select = Con_fit_var_Q[which(Con_fit_var_Q$FPR >= min_FPR & Con_fit_var_Q$FPR <= max_FPR),] # 668
#csvWriter(Con_fit_var_Q_select, "Constant combination matrix.csv")
Con_Q_var_fit_select = Con_Q_var_fit[which(Con_Q_var_fit$FPR >= min_FPR & Con_Q_var_fit$FPR <= max_FPR),] # 668

Con_fit_var_Q_select_matrix = Con_fit_var_Q_select[, c(1, 6, 5,7)]
Con_Q_var_fit_select_matrix = Con_Q_var_fit_select[, c(1, 6, 5,7)]
matrix_FPR_TPR_final = matrix_FPR_TPR_select[, c(1,6,7,8)]
matrix_FPR_TPR_final[,1] = "Dynamic combination"
colnames(matrix_FPR_TPR_final) = c("Method", "FPR", "TPR", "PPV")
colnames(Con_fit_var_Q_select_matrix) = c("Method", "FPR", "TPR", "PPV")
Con_fit_var_Q_select_matrix[,1] = "Constant combination"

matrix_all = rbind(Con_fit_var_Q_select_matrix, matrix_FPR_TPR_final)
t.test(matrix_FPR_TPR_final[,3], Con_fit_var_Q_select_matrix[,3],alternative = "g") # p-value < 2.2e-16
t.test(matrix_FPR_TPR_final[,4], Con_fit_var_Q_select_matrix[,4],alternative = "g") # p-value < 2.2e-16

library(ggplot2)
library(gridExtra)
p1 <- ggplot(matrix_all, aes(x = Method, y = TPR, group = Method, col = Method))+
        geom_violin() +
        #geom_boxplot(width = 0.1, col ="black", alpha = 0.3)+
        geom_jitter(position=position_jitter(0.08), size = 0.6, alpha = 0.3) + 
        stat_summary(fun.y="median", geom="point", col = apple_colors[11], shape = 23, size = 2)+
        scale_color_manual(name = "", breaks = c('Constant combination', "Dynamic combination"),
                           values  = apple_colors[c(5,7)]) +theme(legend.key=element_blank()) +
        scale_y_continuous(name = "True positive rate", 
                           limits=c(0.05, 0.5),
                           breaks = seq(0.05,0.5, by =0.05),
                           labels = seq(0.05,0.5, by= 0.05))+
        theme(legend.position =c(0.7, 0.25), legend.key=element_blank(), legend.text=element_text(size=14)) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_blank(), axis.title.x= element_blank(),axis.line.x = element_blank(),
              axis.ticks.x = element_blank(), axis.text.y.left = element_text(size = 14, color = "black")) +
        theme(text = element_text(size=14))


p2 <- ggplot(matrix_all, aes(x = Method, y = PPV, group = Method, col = Method))+
        geom_violin() +
        geom_jitter(position=position_jitter(0.08), size = 0.6, alpha = 0.3) + 
        stat_summary(fun.y="median", geom="point", col = apple_colors[11], shape = 23, size = 2)+
        scale_color_manual(name = "", breaks = c('Constant combination', "Dynamic combination"),
                           values  = apple_colors[c(5,7)]) +theme(legend.key=element_blank()) +
        scale_y_continuous(name = "Positive predictive value", 
                           limits=c(0.3, 0.8),
                           breaks = seq(0.3,0.8, by =0.1),
                           labels = seq(0.3,0.8, by= 0.1))+
        theme(legend.position =c(0.7, 0.25), legend.key=element_blank(), legend.text=element_text(size=14))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_blank(), axis.title.x= element_blank(),axis.line.x = element_blank(),
              axis.ticks.x = element_blank(), axis.text.y.left = element_text(size = 14, color = "black"))+
        theme(text = element_text(size=14))

g <- arrangeGrob(p1, p2, nrow=1)
#ggsave("TPR and PPV of bestPPV threshold comparison.pdf", plot = g, width = 10, height = 5)
ggsave("~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/SD_comparison.pdf", 
       plot = g, width = 10, height = 5)

##### Write a function to make this plot for each environment

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

Create_matrix_plot_dynamic = function(Neg_number_PPI, Neg_ref_number, output_plot){
        reference_name = reference_name_generation(Neg_number_PPI, Neg_ref_number)
        # Combine 50 datasets into one matrix
        matrix_fit_final = rep(0,7)
        matrix_Q_value_final = rep(0,7)
        for(k in 1:length(reference_name)){
                matrix_ref = dataFrameReader_T(reference_name[k])
                Q_values = seq(-4, 0, by = 0.1)
                Fitness = seq(0, 1, by = 0.05)
                Q_values_select = seq(-4, 0, by = 0.1)
                Fitness_select = seq(0, 0.5, by = 0.05)
                matrix_fit = rep(0,7)
                for (i in 1:length(Fitness_select)){
                        pos_count = length(which(matrix_ref[,1] == "Positive"))
                        neg_count = length(which(matrix_ref[,1] == "Negative"))
                        pos_neg_ratio = pos_count/neg_count
                        matrix_ref_select = ROC_matrix(matrix_ref, Fitness_select[i], "empty", pos_neg_ratio, 
                                                       paste(as.character(k), as.character(i), sep= "_"),
                                                       Fitness, Q_values)
                        matrix_fit = rbind(matrix_fit, matrix_ref_select)
                }
                matrix_fit = matrix_fit[2:nrow(matrix_fit),]
                matrix_fit_final = rbind(matrix_fit_final, matrix_fit)
                
                matrix_Q_value = rep(0,7)
                for (i in 1:length(Q_values_select)){
                        pos_count = length(which(matrix_ref[,1] == "Positive"))
                        neg_count = length(which(matrix_ref[,1] == "Negative"))
                        pos_neg_ratio = pos_count/neg_count
                        matrix_ref_select = ROC_matrix(matrix_ref, "empty", Q_values_select[i], pos_neg_ratio, 
                                                       paste(as.character(k), as.character(i), sep= "_"),
                                                       Fitness, Q_values)
                        matrix_Q_value = rbind(matrix_Q_value, matrix_ref_select)
                }
                matrix_Q_value = matrix_Q_value[2:nrow(matrix_Q_value),]
                matrix_Q_value_final = rbind(matrix_Q_value_final, matrix_Q_value)
        }
        matrix_fit_final = matrix_fit_final[2:nrow(matrix_fit_final),]
        matrix_Q_value_final = matrix_Q_value_final[2:nrow(matrix_Q_value_final),]
        csvWriter(matrix_fit_final, "Specific_fitness_different_Q_value_ROC_curve.csv")
        csvWriter(matrix_Q_value_final, "Specific_Q_value_different_fitness_ROC_curve.csv")
        ### make the plot
        Con_fit_var_Q = dataFrameReader_T("Specific_fitness_different_Q_value_ROC_curve.csv")
        Con_Q_var_fit = dataFrameReader_T("Specific_Q_value_different_fitness_ROC_curve.csv")
        matrix_FPR_TPR = dataFrameReader_T("PPV_threshold_all_data.csv")
        metrics_matrix = dataFrameReader_T("Threshold_metrics_final.csv")
        best_metrics_threshold = metrics_matrix[which(metrics_matrix$MCC == max(metrics_matrix$MCC)),1] 
        matrix_FPR_TPR_chosen = matrix_FPR_TPR[which(matrix_FPR_TPR$PPV == best_metrics_threshold),] # 50
        min_FPR = min(matrix_FPR_TPR_chosen$FPR) # 0.001832401
        max_FPR = max(matrix_FPR_TPR_chosen$FPR) # 0.003860438
        
        matrix_FPR_TPR_select = matrix_FPR_TPR[which(matrix_FPR_TPR$FPR >= min_FPR & matrix_FPR_TPR$FPR <= max_FPR),] # 266
        Con_fit_var_Q_select = Con_fit_var_Q[which(Con_fit_var_Q$FPR >= min_FPR & Con_fit_var_Q$FPR <= max_FPR),] # 668
        #csvWriter(Con_fit_var_Q_select, "Constant combination matrix.csv")
        Con_Q_var_fit_select = Con_Q_var_fit[which(Con_Q_var_fit$FPR >= min_FPR & Con_Q_var_fit$FPR <= max_FPR),] # 668
        
        Con_fit_var_Q_select_matrix = Con_fit_var_Q_select[, c(1, 6, 5,7)]
        Con_Q_var_fit_select_matrix = Con_Q_var_fit_select[, c(1, 6, 5,7)]
        matrix_FPR_TPR_final = matrix_FPR_TPR_select[, c(1,6,7,8)]
        matrix_FPR_TPR_final[,1] = "Dynamic combination"
        colnames(matrix_FPR_TPR_final) = c("Method", "FPR", "TPR", "PPV")
        colnames(Con_fit_var_Q_select_matrix) = c("Method", "FPR", "TPR", "PPV")
        Con_fit_var_Q_select_matrix[,1] = "Constant combination"
        
        matrix_all = rbind(Con_fit_var_Q_select_matrix, matrix_FPR_TPR_final)
        t.test(matrix_FPR_TPR_final[,3], Con_fit_var_Q_select_matrix[,3],alternative = "g") # p-value < 2.2e-16
        t.test(matrix_FPR_TPR_final[,4], Con_fit_var_Q_select_matrix[,4],alternative = "g") # p-value < 2.2e-16
        
        library(ggplot2)
        library(gridExtra)
        p1 <- ggplot(matrix_all, aes(x = Method, y = TPR, group = Method, col = Method))+
                geom_violin() +
                #geom_boxplot(width = 0.1, col ="black", alpha = 0.3)+
                geom_jitter(position=position_jitter(0.08), size = 0.6, alpha = 0.3) + 
                stat_summary(fun.y="median", geom="point", col = apple_colors[11], shape = 23, size = 2)+
                scale_color_manual(name = "", breaks = c('Constant combination', "Dynamic combination"),
                                   values  = apple_colors[c(5,7)]) +theme(legend.key=element_blank()) +
                scale_y_continuous(name = "True positive rate", 
                                   limits=c(0.05, 0.5),
                                   breaks = seq(0.05,0.5, by =0.05),
                                   labels = seq(0.05,0.5, by= 0.05))+
                theme(legend.position =c(0.7, 0.75), legend.key=element_blank(), legend.text=element_text(size=14)) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                theme(axis.text.x = element_blank(), axis.title.x= element_blank(),axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(), axis.text.y.left = element_text(size = 14, color = "black")) +
                theme(text = element_text(size=14))
        
        
        p2 <- ggplot(matrix_all, aes(x = Method, y = PPV, group = Method, col = Method))+
                geom_violin() +
                geom_jitter(position=position_jitter(0.08), size = 0.6, alpha = 0.3) + 
                stat_summary(fun.y="median", geom="point", col = apple_colors[11], shape = 23, size = 2)+
                scale_color_manual(name = "", breaks = c('Constant combination', "Dynamic combination"),
                                   values  = apple_colors[c(5,7)]) +theme(legend.key=element_blank()) +
                scale_y_continuous(name = "Positive predictive value", 
                                   limits=c(0.1, 0.5),
                                   breaks = seq(0.1,0.5, by =0.1),
                                   labels = seq(0.1,0.5, by= 0.1))+
                theme(legend.position =c(0.7, 0.25), legend.key=element_blank(), legend.text=element_text(size=14))+
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                theme(axis.text.x = element_blank(), axis.title.x= element_blank(),axis.line.x = element_blank(),
                      axis.ticks.x = element_blank(), axis.text.y.left = element_text(size = 14, color = "black"))+
                theme(text = element_text(size=14))
        
        g <- arrangeGrob(p1, p2, nrow=1)
        #ggsave("TPR and PPV of bestPPV threshold comparison.pdf", plot = g, width = 10, height = 5)
        ggsave(output_plot, plot = g, width = 10, height = 5)
}

### SD2
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_2/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/SD2_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### SD_merge
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO_merge/reference_set/p_value/")
Neg_number_PPI = 6.4e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/SD_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### H2O2
setwd("/Volumes/zmliu_02/PPiseq_03/H2O2/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/H2O2_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### HU
setwd("/Volumes/zmliu_02/PPiseq_03/HU/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/HU_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### Forskolin
setwd("/Volumes/zmliu_02/PPiseq_03/Forskolin/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/Forskolin_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### Dox
setwd("/Volumes/zmliu_02/PPiseq_03/Dox/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/Dox_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### NaCl
setwd("/Volumes/zmliu_02/PPiseq_03/NaCl_0.4M/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/NaCl_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### 16C
setwd("/Volumes/zmliu_02/PPiseq_03/16C/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/16C_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### FK506
setwd("/Volumes/zmliu_02/PPiseq_03/FK506/reference_set/p_value/")
Neg_number_PPI = 6e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/FK506_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)

### Raffinose
setwd("/Volumes/zmliu_02/PPiseq_03/Raffinose/reference_set/p_value/")
Neg_number_PPI = 4.5e4
Neg_ref_number = 50
output_plot = "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM5_TPR_PPV_Best_comparison/Raffinose_merge_comparison.pdf"
Create_matrix_plot_dynamic(Neg_number_PPI, Neg_ref_number, output_plot)







