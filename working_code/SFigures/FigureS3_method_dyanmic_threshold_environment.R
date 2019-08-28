source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")


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
        coef= lm(fitness_threshold~poly(Q_value_threshold,3, raw = TRUE))$coefficients
        simulated_Q = unique(Q_value_threshold)
        simulated_fit = rep(0, length(simulated_Q))
        for(i in 1:length(simulated_Q)){
                m = simulated_Q[i]
                n = c(1, m, m^2, m^3)
                simulated_fit[i] = sum(n * coef)
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
                geom_line(aes(Q_value, Fitness), col = apple_colors[11],simulated_data) +
                #scale_color_manual('', breaks = c("Discrete_combination", "Dynamic_combination"),
                                   #values = apple_colors[c(8,7)], 
                                   #guide = guide_legend(override.aes = list(linetype = c(NA,1), shape = c(16, NA)))) +
                #scale_x_continuous(name= expression(paste("Log10(", italic(q), "): ", "Q",sep = "")),
                scale_x_continuous(name= expression(italic(q)),
                                   limits = c(-1.4, 0),
                                   breaks = seq(-1.4, 0, by = 0.2),
                                   labels = seq(-1.4, 0, by = 0.2)) +
                ylab(expression(italic(f))) +
                
                theme(legend.key=element_blank(), legend.position = c(0.2,0.7)) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                theme(axis.text.x = element_text(size = 10, color = "black"), 
                      axis.title=element_text(size=14,face="bold"),
                      axis.text.y.left = element_text(size = 10, color = "black")) + 
                theme(text = element_text(size=10))
        
        ggsave(output, width = 3, height =3, device = "pdf")
        
        
}
# Check the linear regression effect of different FPR
PPV_threshold_files = rep("Name", 50)
for (i in 1:50){
        PPV_threshold_files[i] = paste("Fitness_Q_value_PPV_threshold", as.character(i), "data.csv", sep = "_")
}

setwd("/Volumes/zmliu_02/PPiseq/DMSO/reference_set/")
plot_threshold("0.71", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/SD_threshold.pdf")

#setwd("/Volumes/zmliu_02/PPiseq/DMSO_2/reference_set/")
#plot_threshold("0.72", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/SD_2_threshold.pdf")

setwd("/Volumes/zmliu_02/PPiseq/Forskolin/reference_set/one_binary/PPV_threshold/")
plot_threshold("0.68", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/Forskolin_threshold.pdf")

setwd("/Volumes/zmliu_02/PPiseq/FK506/reference_set/one_binary/PPV_threshold/")
plot_threshold("0.7", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/FK506_threshold.pdf")

#NaCl Bionomial
#setwd("/Volumes/zmliu_02/PPiseq/NaCl_0.4M/reference_set/one_binary/PPV_threshold/")
#plot_threshold("0.7", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/NaCl_threshold.pdf")

#Raffinose
setwd("/Volumes/zmliu_02/PPiseq/Raffinose/reference_set/one_binary/PPV_threshold/")
plot_threshold("0.44", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/Raffinose_threshold.pdf")

## HU
setwd("/Volumes/zmliu_02/PPiseq/HU/reference_set/one_binary/PPV_threshold/")
plot_threshold("0.7", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/HU_threshold.pdf")

## H2O2
setwd("/Volumes/zmliu_02/PPiseq/H2O2/reference_set/one_binary/PPV_threshold/")
plot_threshold("0.7", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/H2O2_threshold.pdf")

## Dox
setwd("/Volumes/zmliu_02/PPiseq/Dox/reference_set/one_binary/PPV_threshold/")
plot_threshold("0.74", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/Dox_threshold.pdf")

## 16C Binomial
#setwd("/Volumes/zmliu_02/PPiseq/16C/reference_set/one_binary/PPV_threshold/")
#plot_threshold("0.41", "~/Dropbox/PPiSeq_02/Working_figure/SFigures/FigureS3_Method_threshold_do_each_environment/16C_threshold.pdf")
