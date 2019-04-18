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

#### Caculate an variation score for each PPI based on their fitness values across
#### different environments. Use the normalized value and Keep the primary normalized fitess even a PPI is called negative
setwd("~/Dropbox/PPiSeq_02/")

PPI_fit_norm = csvReader_T("Working_data/Positive_PPI_environment/Pos_PPI_normalized_fit_primary.csv") # 14564
PPI_dup = mark_duplicates_fast(PPI_fit_norm[,1]) #13430
matrix = matrix(NA, nrow(PPI_dup), 2* ncol(PPI_fit_norm))
for(i in 1:nrow(PPI_dup)){
        if (PPI_dup[i,2] != "0"){
                matrix[i,1:10] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
                matrix[i,11:20] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,2]),]
        }else{
                matrix[i,1:10] = PPI_fit_norm[which(PPI_fit_norm[,1] == PPI_dup[i,1]),]
        }
}
colnames(matrix) = c("PPI", "DMSO", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", 
                     "NaCl", "16C", "FK506", "PPI_opposite", "DMSO", "H2O2", "HU",
                     "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
mean_fitness = rep(0, nrow(matrix))
for(i in 1:length(mean_fitness)){
        mean_fitness[i] = mean(na.omit(as.numeric(matrix[i, c(2:10, 12:20)])))
}

matrix_final = cbind(matrix[,1], mean_fitness, matrix[,2:ncol(matrix)])
csvWriter(matrix_final, "Working_data/Positive_PPI_environment/Normalized_fitness_PPI_all_primary.csv")

PPI_norm = csvReader_T("Working_data/Positive_PPI_environment/Normalized_fitness_PPI_all_primary.csv")
PPI_norm_matrix = matrix(0, nrow(PPI_norm), 10)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
        for(j in 3:11){
                mean_fit = mean(as.numeric(na.omit(PPI_norm[i,c(j, j + 10)])))
                if(is.na(mean_fit)){
                        PPI_norm_matrix[i, j-1] = 0
                }
                else{
                        PPI_norm_matrix[i,j-1] = mean_fit
                }
                
        }  
}
variation_score = rep(0, nrow(PPI_norm_matrix))
environment_number = rep(0, nrow(PPI_norm_matrix))
for(i in 1:length(variation_score)){
        fitness = as.numeric(PPI_norm_matrix[i,2:10])
        variation_score[i] = sd(fitness)/mean(fitness)
        environment_number[i] = length(which(fitness != 0))
}
# Include the variation score into the matrix, and output the matrix
PPI_norm_matrix_final = cbind(PPI_norm_matrix[,1],environment_number, variation_score,
                              PPI_norm_matrix[,2:10])
colnames(PPI_norm_matrix_final) = c("PPI", "Environment_number", "Variation_score",
                                    "DMSO", "H2O2", "HU", "Dox", "Forskolin", "Raffinose",
                                    "NaCl", "16C", "FK506")
csvWriter(PPI_norm_matrix_final, "Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
cor(environment_number, variation_score, method="spearman") # -0.7742618
# Check the correlation between variation score and environments_number
pdf("Working_figure/Figure6/SFigure6_1_correlation_variation_score_environment_number.pdf", width =5, height =5)
plot(environment_number, variation_score, type = "p", ylim = c(0, 10),
     xlab = "Number of environments in which a PPI is observed",
     ylab = "Variation score", bty = "n", col = rgb(0,0,1,alpha = 0.3))
axis(1, at= 1:9, labels = as.character(1:9))
text(3, 0.1, expression(paste("Spearman's ", italic(r), " = -0.77")))
dev.off()

pdf("Working_figure/Figure6/SFigure6_2_histogram_variation_score_all.pdf", width =5, height =5)
hist(variation_score, breaks = seq(-1400, 5200, by = 0.02), xlim = c(0,4),xlab = "Varation Score", main = NA, col = apple_colors[5])
dev.off()

pdf("Working_figure/Figure6/SFigure6_2_histogram_variation_score_multiple_environment.pdf", width =5, height =5)
hist(variation_score[which(environment_number !=1)], breaks = seq(-1400, 5200, by = 0.02), 
     xlab = "Varation Score", main = NA, col = apple_colors[5], xlim = c(0,2.5))
dev.off()

#Check correlation between variation score and gene features
setwd("~/Dropbox/PPiSeq_02/")
variation_score = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
PP_pair = split_string_vector(variation_score[,1])
protein_unique = unique(as.vector(PP_pair)) # 2111
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
        index_protein = unique(c(which(PP_pair[,1] == protein_unique[i]),
                                 which(PP_pair[,2] == protein_unique[i])))
        vScore_protein[i,2] = mean(as.numeric(variation_score[index_protein,3]))
}
gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39

gene_feature_matched = gene_feature[match(vScore_protein[,1], gene_feature[,1]),]
features = colnames(gene_feature)
features_chosen = features[2:(length(features)-7)]
vScore = vScore_protein[,2]
co_efficient = rep(0, length(features_chosen))
for(i in 1:length(features_chosen)){
        feature_column = which(features == features_chosen[i])
        feature_data = gene_feature_matched[, feature_column]
        cor_matrix = data.frame(vScore, feature_data) # 2111
        cor_matrix_filter = na.omit(cor_matrix) # 2109
        co_efficient[i] = cor(as.numeric(cor_matrix_filter[,1]), 
                              as.numeric(cor_matrix_filter[,2]), method = 'spearman')
        
}
coefficient_matrix = cbind(features_chosen, co_efficient)
coefficient_matrix_order = coefficient_matrix[order(as.numeric(coefficient_matrix[,2]), decreasing = F),]
coefficient_order = as.numeric(coefficient_matrix_order[,2])
features_order = coefficient_matrix_order[,1]
pdf("Working_figure/Figure6/Figure6A_Correlation_variation_score_gene_features.pdf", width = 5, height = 5)
par(mar= c(5,5,1,1))
barCenter = barplot(coefficient_order, horiz=T, beside=F, xlim=c(-0.4, 0.2), 
                    space= rep(0.4, 31),
                    xlab="Spearman correlation with protein variation score",
                    axisnames=F, border=NA, cex.axis=0.7, col = apple_colors[5])
text(-0.4, barCenter, labels= features_order, cex = 0.5, xpd = TRUE)
dev.off()

### Check the distribution of variation scores for binary features
setwd("~/Dropbox/PPiSeq_02/")
variation_score = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment.csv")
PP_pair = split_string_vector(variation_score[,1])
protein_unique = unique(as.vector(PP_pair)) # 2111
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
        index_protein = unique(c(which(PP_pair[,1] == protein_unique[i]),
                                 which(PP_pair[,2] == protein_unique[i])))
        vScore_protein[i,2] = mean(as.numeric(variation_score[index_protein,3]))
}

gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39
gene_feature_matched = gene_feature[match(vScore_protein[,1], gene_feature[,1]),]
features = colnames(gene_feature)
features_chosen = features[2:length(features)]
gene_feature_binary = features_chosen[(length(features_chosen)-6):length(features_chosen)]

library(ggplot2)
binary_feature_split = function(binary_feature_chosen, features,
                                gene_feature_matched, vScore_protein){
        feature_column = which(features == binary_feature_chosen)
        binary_variation = na.omit(cbind(vScore_protein, gene_feature_matched[,feature_column]))
        PPI = binary_variation[,1]
        binary = as.numeric(binary_variation[,3])
        vScore = as.numeric(binary_variation[,2])
        binary_variation_data = data.frame(PPI, binary, vScore)
        #binary_number = as.numeric(unique(binary))
        vScore_1 = binary_variation_data[which(binary_variation_data$binary == 0),3]
        vScore_2 = binary_variation_data[which(binary_variation_data$binary == 1),3]
        p_value= t.test(vScore_1, vScore_2)$p.value
        binary_variation_data$binary = factor(binary_variation_data$binary, levels = c("0", "1"))
        ggplot()+
                geom_violin(aes(x = binary, y = vScore), binary_variation_data, 
                            draw_quantiles = c(0.25, 0.5, 0.75), col = apple_colors[5],
                            show.legend = FALSE)+
                stat_summary(aes(x = binary, y = vScore),binary_variation_data,
                             fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
                xlab(binary_feature_chosen) +
                ylab("Variation score") +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                theme(axis.text.x = element_text(size = 10, color = "black"), 
                      axis.text.y.left = element_text(size = 10, color = "black"))
        ggsave(paste(binary_feature_chosen, "violin_plot.pdf", sep= "_"), width = 5, height = 5)
        return(p_value)
}
setwd("~/Dropbox/PPiSeq_02/Working_figure/Figure6/binary_feature/")
p.value = rep(0, length(gene_feature_binary))
for(i in 1:length(gene_feature_binary)){
        p.value[i] = binary_feature_split(gene_feature_binary[i], features,gene_feature_matched, vScore_protein)
}
p.value.matrix = cbind(gene_feature_binary, p.value)
csvWriter(p.value.matrix, "P.value_matrix_difference.csv")

#################################################################################################################
#### Caculate an variation score for each PPI based on their fitness values across
#### different environments. Use the normalized value and consider the negative PPI as 0
setwd("~/Dropbox/PPiSeq_02/")
PPI_norm = csvReader_T("Working_data/Positive_PPI_environment/Normalized_fitness_PPI_all.csv")
PPI_norm_matrix = matrix(0, nrow(PPI_norm), 10)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
        for(j in 3:11){
                mean_fit = mean(as.numeric(na.omit(PPI_norm[i,c(j, j + 10)])))
                if(is.na(mean_fit)){
                        PPI_norm_matrix[i, j-1] = 0
                }
                else{
                        PPI_norm_matrix[i,j-1] = mean_fit
                }
                 
        }  
}
variation_score = rep(0, nrow(PPI_norm_matrix))
environment_number = rep(0, nrow(PPI_norm_matrix))
for(i in 1:length(variation_score)){
        fitness = as.numeric(PPI_norm_matrix[i,2:10])
        variation_score[i] = sd(fitness)/mean(fitness)
        environment_number[i] = length(which(fitness != 0))
}
# Include the variation score into the matrix, and output the matrix
PPI_norm_matrix_final = cbind(PPI_norm_matrix[,1],environment_number, variation_score,
                              PPI_norm_matrix[,2:10])
colnames(PPI_norm_matrix_final) = c("PPI", "Environment_number", "Variation_score",
                                    "DMSO", "H2O2", "HU", "Dox", "Forskolin", "Raffinose",
                                    "NaCl", "16C", "FK506")
csvWriter(PPI_norm_matrix_final, "Working_data/Positive_PPI_environment/Variation_score_PPI_environment.csv")
cor(environment_number, variation_score, method="spearman") # -0.9254962
# Check the correlation between variation score and environments_number
pdf("Working_figure/Figure6/SFigure6_1_correlation_variation_score_environment_number.pdf", width =5, height =5)
plot(environment_number, variation_score, type = "p", 
     xlab = "Number of environments in which a PPI is observed",
     ylab = "Variation score", bty = "n", col = rgb(0,0,1,alpha = 0.3))
axis(1, at= 1:9, labels = as.character(1:9))
text(7, 2.8, expression(paste("Spearman's ", italic(r), " = -0.92")))
dev.off()

pdf("Working_figure/Figure6/SFigure6_2_histogram_variation_score_all.pdf", width =5, height =5)
hist(variation_score, breaks = seq(0, 3, by = 0.02), xlab = "Varation Score", main = NA, col = apple_colors[5])
dev.off()

pdf("Working_figure/Figure6/SFigure6_2_histogram_variation_score_multiple_environment.pdf", width =5, height =5)
hist(variation_score[which(environment_number !=1)], breaks = seq(0, 3, by = 0.02), 
     xlab = "Varation Score", main = NA, col = apple_colors[5], xlim = c(0,2.5))
dev.off()

#Check correlation between variation score and gene features
setwd("~/Dropbox/PPiSeq_02/")
variation_score = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment.csv")
PP_pair = split_string_vector(variation_score[,1])
protein_unique = unique(as.vector(PP_pair)) # 2111
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
        index_protein = unique(c(which(PP_pair[,1] == protein_unique[i]),
                                 which(PP_pair[,2] == protein_unique[i])))
        vScore_protein[i,2] = mean(as.numeric(variation_score[index_protein,3]))
}
gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39

gene_feature_matched = gene_feature[match(vScore_protein[,1], gene_feature[,1]),]
features = colnames(gene_feature)
features_chosen = features[2:(length(features)-7)]
vScore = vScore_protein[,2]
co_efficient = rep(0, length(features_chosen))
for(i in 1:length(features_chosen)){
        feature_column = which(features == features_chosen[i])
        feature_data = gene_feature_matched[, feature_column]
        cor_matrix = data.frame(vScore, feature_data) # 2111
        cor_matrix_filter = na.omit(cor_matrix) # 2109
        co_efficient[i] = cor(as.numeric(cor_matrix_filter[,1]), 
                              as.numeric(cor_matrix_filter[,2]), method = 'spearman')
        
}
coefficient_matrix = cbind(features_chosen, co_efficient)
coefficient_matrix_order = coefficient_matrix[order(as.numeric(coefficient_matrix[,2]), decreasing = F),]
coefficient_order = as.numeric(coefficient_matrix_order[,2])
features_order = coefficient_matrix_order[,1]
pdf("Working_figure/Figure6/Figure6A_Correlation_variation_score_gene_features.pdf", width = 5, height = 5)
par(mar= c(5,5,1,1))
barCenter = barplot(coefficient_order, horiz=T, beside=F, xlim=c(-0.6, 0.3), 
                    space= rep(0.4, 31),
                    xlab="Spearman correlation with protein variation score",
                    axisnames=F, border=NA, cex.axis=0.7, col = apple_colors[5])
text(-0.6, barCenter, labels= features_order, cex = 0.5, xpd = TRUE)
dev.off()

### Check the distribution of variation scores for binary features
setwd("~/Dropbox/PPiSeq_02/")
variation_score = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment.csv")
PP_pair = split_string_vector(variation_score[,1])
protein_unique = unique(as.vector(PP_pair)) # 2111
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
        index_protein = unique(c(which(PP_pair[,1] == protein_unique[i]),
                                 which(PP_pair[,2] == protein_unique[i])))
        vScore_protein[i,2] = mean(as.numeric(variation_score[index_protein,3]))
}

gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39
gene_feature_matched = gene_feature[match(vScore_protein[,1], gene_feature[,1]),]
features = colnames(gene_feature)
features_chosen = features[2:length(features)]
gene_feature_binary = features_chosen[(length(features_chosen)-6):length(features_chosen)]

library(ggplot2)
binary_feature_split = function(binary_feature_chosen, features,
                                gene_feature_matched, vScore_protein){
        feature_column = which(features == binary_feature_chosen)
        binary_variation = na.omit(cbind(vScore_protein, gene_feature_matched[,feature_column]))
        PPI = binary_variation[,1]
        binary = as.numeric(binary_variation[,3])
        vScore = as.numeric(binary_variation[,2])
        binary_variation_data = data.frame(PPI, binary, vScore)
        #binary_number = as.numeric(unique(binary))
        vScore_1 = binary_variation_data[which(binary_variation_data$binary == 0),3]
        vScore_2 = binary_variation_data[which(binary_variation_data$binary == 1),3]
        p_value= t.test(vScore_1, vScore_2)$p.value
        binary_variation_data$binary = factor(binary_variation_data$binary, levels = c("0", "1"))
        ggplot()+
                geom_violin(aes(x = binary, y = vScore), binary_variation_data, 
                            draw_quantiles = c(0.25, 0.5, 0.75), col = apple_colors[5],
                            show.legend = FALSE)+
                stat_summary(aes(x = binary, y = vScore),binary_variation_data,
                             fun.y="mean", geom="point", col = apple_colors[11], shape = 23, size = 1)+
                xlab(binary_feature_chosen) +
                ylab("Variation score") +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                theme(axis.text.x = element_text(size = 10, color = "black"), 
                      axis.text.y.left = element_text(size = 10, color = "black"))
        ggsave(paste(binary_feature_chosen, "violin_plot.pdf", sep= "_"), width = 5, height = 5)
        return(p_value)
}
setwd("~/Dropbox/PPiSeq_02/Working_figure/Figure6/binary_feature/")
p.value = rep(0, length(gene_feature_binary))
for(i in 1:length(gene_feature_binary)){
        p.value[i] = binary_feature_split(gene_feature_binary[i], features,gene_feature_matched, vScore_protein)
}
p.value.matrix = cbind(gene_feature_binary, p.value)
csvWriter(p.value.matrix, "P.value_matrix_difference.csv")