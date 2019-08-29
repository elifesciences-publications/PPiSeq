source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
# Scatter plot to show the TPR value of 50 datasets with FPR = 5e-4 with three data sets
setwd("/Volumes/zmliu_02/PPiseq/DMSO/reference_set/")
Con_fit_var_Q = dataFrameReader_T("Specific_fitness_different_Q_value_ROC_curve_new.csv")
Con_Q_var_fit = dataFrameReader_T("Specific_Q_value_different_fitness_ROC_curve_new.csv")
matrix_FPR_TPR = dataFrameReader_T("PPV_threshold_all_data_new.csv")
metrics_matrix = dataFrameReader_T("Threshold_metrics_final.csv")
best_metrics_threshold = metrics_matrix[which(metrics_matrix$MCC == max(metrics_matrix$MCC)),1] # 0.71
matrix_FPR_TPR_chosen = matrix_FPR_TPR[which(matrix_FPR_TPR$PPV == best_metrics_threshold),] # 50
min_FPR = min(matrix_FPR_TPR_chosen$FPR) # 0.001832401
max_FPR = max(matrix_FPR_TPR_chosen$FPR) # 0.003860438

matrix_FPR_TPR_select = matrix_FPR_TPR[which(matrix_FPR_TPR$FPR >= min_FPR & matrix_FPR_TPR$FPR <= max_FPR),] # 266
Con_fit_var_Q_select = Con_fit_var_Q[which(Con_fit_var_Q$FPR >= min_FPR & Con_fit_var_Q$FPR <= max_FPR),] # 668
#csvWriter(Con_fit_var_Q_select, "Constant combination matrix.csv")
Con_Q_var_fit_select = Con_Q_var_fit[which(Con_Q_var_fit$FPR >= min_FPR & Con_Q_var_fit$FPR <= max_FPR),] # 668

Con_fit_var_Q_select_matrix = Con_fit_var_Q_select[, c(1, 7, 6,8)]
Con_Q_var_fit_select_matrix = Con_Q_var_fit_select[, c(1, 7, 6,8)]
matrix_FPR_TPR_final = matrix_FPR_TPR_select[, c(1,8,9,10)]
matrix_FPR_TPR_final[,1] = "Dynamic combination"
colnames(matrix_FPR_TPR_final) = c("Method", "FPR", "TPR", "PPV")
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
ggsave("~/Dropbox/PPiSeq_02/Working_figure/SFigures/SFigureSM5_TPR_PPV_Best_comparison.pdf", 
       plot = g, width = 10, height = 5)