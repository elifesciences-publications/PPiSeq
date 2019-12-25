# Source some basic functions froma function.R in Github repository
source_https <- function(u, unlink.tmp.certs = FALSE) {
        # load package
        require(RCurl)
        # read script lines from website using a security certificate
        if(!file.exists("cacert.pem")) download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.pem")
        script <- getURL(u, followlocation = TRUE, cainfo = "cacert.pem")
        if(unlink.tmp.certs) unlink("cacert.pem")
        
        # parase lines and evealuate in the global environement
        eval(parse(text = script), envir= .GlobalEnv)
}
source_https("https://raw.githubusercontent.com/sashaflevy/PPiSeq/master/working_code/function.R", unlink.tmp.certs = TRUE)

#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/reference_set_large_Q_range/")
## First get the thresholds for larger Q-value range (-6, 0) instead of (-4, 0)
## Use the same reference sets generated before

Q_values = seq(-6, 0, by = 0.1)
Fitness = seq(0, 1.0, by = 0.01)

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
                                                                      log10(matrix_ref[,6]) <= Q_values[j]),1]
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

for (k in 1: length(reference_name)){
        matrix_ref = dataFrameReader_T(reference_name[k])
        PPV_matrix = ROC_line_matrix(matrix_ref, Fitness, Q_values, "PPV")
        csvWriter_rownames(PPV_matrix, paste("Fitness_Q_value_PPV", as.character(k), "data.csv", sep = "_"))
        
        PPV_bin = c(seq(0.5, 0.58, by= 0.02), seq(0.6, 0.8, by = 0.01), seq(0.82, 0.9, by= 0.02))
        matrix_different_PPV = rep(0, 3)
        for (m in 1:length(PPV_bin)){
                matrix_PPV_threshold = Q_fit_PPV_threshold(PPV_matrix, PPV_bin[m])
                matrix_different_PPV = rbind(matrix_different_PPV, matrix_PPV_threshold)
        }
        matrix_different_PPV = matrix_different_PPV[2:nrow(matrix_different_PPV),]
        csvWriter(matrix_different_PPV, paste("Fitness_Q_value_PPV_threshold", as.character(k), "data.csv", sep = "_"))
}


#### Figure SM3 For the same reference set, different PPV thresholds
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/reference_set_large_Q_range/")
PPV_threshold = dataFrameReader_T("Fitness_Q_value_PPV_threshold_1_data.csv")
PPV_chosen = c("PPV_0.5", "PPV_0.56", "PPV_0.6", "PPV_0.64", "PPV_0.71", "PPV_0.75", "PPV_0.8", "PPV_0.9")
PPV_threshold_chosen = PPV_threshold[which(PPV_threshold$PPV %in% PPV_chosen),]
PPV_threshold_chosen$PPV = factor(PPV_threshold_chosen$PPV, levels = rev(levels(PPV_threshold_chosen$PPV)))
col_chosen = rev(c("#3288bd", "#66c2a5", "#abdda4", "#e6f598", "#fee08b", "#fdae61",
               "#f46d43", "#f46d43"))
library("ggplot2")
ggplot()+
        geom_line(aes(x = Q_value, y = Fitness, group = PPV, col = PPV), PPV_threshold_chosen) +
        scale_color_manual("", values = col_chosen )+
        scale_y_continuous(name = expression(italic(f)),
                           limits=c(0, 0.7),
                           breaks=seq(0,0.7, by =0.1),
                           labels = seq(0, 0.7, by= 0.1)) +
        scale_x_continuous(name = expression(italic(p)), 
                           limits=c(-6, 0),
                           breaks=seq(-6, 0, by =1),
                           labels = seq(-6, 0, by= 1))+
        theme(legend.key=element_blank()) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 10, color = "black"),
              axis.title=element_text(size=14,face="bold"),
              axis.text.y.left = element_text(size = 10, color = "black")) + 
        theme(text = element_text(size=12))
ggsave("~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/FigureSM3A_Method_Different_dynamic_PPV_threshold.pdf", width = 5, height = 4)
#ggsave("~/Desktop/FigureS2_Method_Different_dynamic_PPV_threshold.pdf", width = 6, height = 5)


##################### Figure SM3B plot 50 lines for the same PPV = 0.71 from 50 rereference sets

PPV_threshold_files = rep("Name", 50)
for (i in 1:50){
        PPV_threshold_files[i] = paste("Fitness_Q_value_PPV_threshold", as.character(i), "data.csv", sep = "_")
}
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
        fitmodel = nls(fitness_threshold  ~ SSlogis( Q_value_threshold, Asym, xmid, scal))
        coeff = coef(fitmodel)
        simulated_Q = unique(Q_value_threshold)
        simulated_fit = rep(0, length(simulated_Q))
        for(i in 1:length(simulated_Q)){
                m = simulated_Q[i]
                simulated_fit[i] = coeff[1]/(1+exp((coeff[2] - m)/coeff[3]))
        }
        
        simulated_data = cbind(rep("Dynamic_combination", length(simulated_Q)), rep("Simulation", length(simulated_Q)), simulated_Q, simulated_fit)
        colnames(simulated_data) = c("Data_type", "Lineage", "Q_value", "Fitness")
        csvWriter(simulated_data, paste("PPV", as.character(specific_PPV), "simulated_threshold.csv", sep= "_"))
        
        real_data = cbind(rep("Discrete_combination", nrow(PPV_threshold_all)), PPV_threshold_all)
        colnames(real_data) = c("Data_type", "Lineage", "Q_value", "Fitness")
        csvWriter(real_data, paste("PPV", as.character(specific_PPV), "threshold_plot.csv", sep= "_"))
        
        real_data = dataFrameReader_T(paste("PPV", as.character(specific_PPV), "threshold_plot.csv", sep= "_"))
        simulated_data = dataFrameReader_T(paste("PPV", as.character(specific_PPV), "simulated_threshold.csv", sep= "_"))
        #######
        
        ggplot() +
          #geom_point(aes(Q_value, Fitness, color = Data_type), real_data, alpha = 0.5)+
          #geom_line(aes(Q_value, Fitness, color = Data_type), real_data, alpha = 0.5)+
          geom_hex(aes(x= Q_value, y= Fitness), real_data, bins = 50)+
          scale_fill_gradient(low= apple_colors[10], high = apple_colors[7])+
          #scale_fill_gradient(low= "white", high = apple_colors[7])+
          geom_line(aes(Q_value, Fitness), col = apple_colors[11], simulated_data) +
          #scale_color_manual('', breaks = c("Discrete_combination", "Dynamic_combination"),
          #values = apple_colors[c(8,7)], 
          #guide = guide_legend(override.aes = list(linetype = c(NA,1),shape = c(16, NA)))) +
          scale_x_continuous(name= expression(italic(p)),
                             limits = c(-4, 0),
                             breaks = seq(-4, 0, by = 0.5),
                             labels = seq(-4, 0, by = 0.5)) +
          ylab(expression(italic(f)))+
          #scale_y_continuous(name = "Fitness",
          #limits = c(0, 0.6),
          #breaks = seq(0, 0.6, by = 0.1),
          #labels = seq(0, 0.6, by = 0.1)) +
          
          theme(legend.key=element_blank(), legend.position = c(0.9,0.2)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")) +
          theme(axis.text.x = element_text(size = 10, color = "black"), 
                axis.title=element_text(size=14,face="bold"),
                axis.text.y.left = element_text(size = 10, color = "black")) + 
          theme(text = element_text(size=12))
        
        ggsave(output, width = 5, height =5, device = "pdf")
        
        
        
}
# Check the linear regression effect of different FPR
plot_threshold("0.7", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/FigureSM3B_Method_Different_references_PPV_0.7_threshold_dot.pdf")



