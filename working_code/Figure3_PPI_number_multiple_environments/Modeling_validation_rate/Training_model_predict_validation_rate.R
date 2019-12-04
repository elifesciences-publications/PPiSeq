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

setwd("~/Dropbox/PPiseq_02")
count_summary = csvReader_T("Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")
protein_abundance = csvReader_T("Paper_data/Outside_datasets/protein_abundance/table_S4.csv")
protein_pair = split_string_vector(count_summary[,1])
protein_1_abun = as.numeric(protein_abundance[match(protein_pair[,1], protein_abundance[,1]),5])
protein_2_abun = as.numeric(protein_abundance[match(protein_pair[,2], protein_abundance[,1]),5])
protein_12_abun = data.frame(protein_1_abun, protein_2_abun)
mean_abun = rep(0, nrow(protein_12_abun))
#mean_abun = apply()
mean_abun = rowMeans(protein_12_abun, na.rm = T)
mean_abun_count = data.frame(count_summary[,1], mean_abun, count_summary[,2])
colnames(mean_abun_count) = c("PPI", "Mean_abundance", "Env_count")

#col_purple = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
library(ggplot2)
ggplot() +
  #geom_jitter(aes(x = Env_count, y = Mean_abundance, group = Env_count, fill = Env_count, col = Env_count), mean_abun_count, 
              #width = 0.2, alpha = 0.5)+
  geom_boxplot(aes(x = Env_count, y = Mean_abundance), mean_abun_count, outlier.shape=NA) +
  geom_dotplot(aes(x = Env_count, y = Mean_abundance), mean_abun_count, 
               binaxis="y",stackdir="center",binwidth=80, alpha=0.2, col = apple_colors[8]) +
  #geom_violin(aes(x = Env_count, y = Mean_abundance), mean_abun_count,
              #draw_quantiles = c(0.25, 0.5, 0.75)) +
  #scale_x_discrete("Number of environments in which a PPI is identified")+
  #scale_color_manual(name = "", values  = col_purple) +
  xlab("Number of environments in which a PPI is identified") +
  ylab("Average molecular number per cell of a protein pair") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black")) 
  #theme(plot.margin = unit(c(1,1,2,1), "cm"))
ggsave("~/Desktop/Average_protein_abundance_PPI_environment_number.pdf", width=5, height =5)


### Check the normalized fitness values for these PPIs
vScore = dataFrameReader_T("Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
vScore_fit = vScore[,4:ncol(vScore)]


### Take the mean fitness across all environments for all the positive PPIs
mean_fitness_pos = rowMeans(vScore_fit)
min(mean_fitness_pos) # 0.01263292
fitness_count = data.frame(mean_fitness_pos, count_summary[,2])
colnames(fitness_count) = c("Mean_fitness", "Env_count")
ggplot() +
  
  geom_boxplot(aes(x = Env_count, y = Mean_fitness), fitness_count, outlier.shape=NA) +
  geom_dotplot(aes(x = Env_count, y = Mean_fitness), fitness_count, 
               binaxis="y",stackdir="center",binwidth=0.0015, alpha=0.2, col = apple_colors[8]) +
  xlab("Number of environments in which a PPI is identified") +
  ylab("Mean fitness of a PPI across different environments") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black")) 
#theme(plot.margin = unit(c(1,1,2,1), "cm"))
ggsave("~/Desktop/All_Mean_fitness_PPI_environment_number.pdf", width=5, height =5)


### Take the mean fitness for all the positive PPIs
mean_fitness_pos = rep(0, nrow(vScore))
for(i in 1:length(mean_fitness_pos)){
  pos_env_index = which(count_summary[i,3:ncol(count_summary)] == "1")
  pos_env_index_fit = pos_env_index + 3
  mean_fitness_pos[i] = mean(as.numeric(vScore[i, pos_env_index_fit]))
}
min(mean_fitness_pos) # 0.1121392
fitness_count = data.frame(mean_fitness_pos, count_summary[,2])
colnames(fitness_count) = c("Mean_fitness", "Env_count")
ggplot() +
  
  geom_boxplot(aes(x = Env_count, y = Mean_fitness), fitness_count, outlier.shape=NA) +
  geom_dotplot(aes(x = Env_count, y = Mean_fitness), fitness_count, 
               binaxis="y",stackdir="center",binwidth=0.002, alpha=0.2, col = apple_colors[8]) +
  xlab("Number of environments in which a PPI is identified") +
  ylab("Mean fitness of a PPI across different environments") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black")) 
#theme(plot.margin = unit(c(1,1,2,1), "cm"))
ggsave("~/Desktop/Mean_fitness_PPI_environment_number.pdf", width=5, height =5)

### Take all the positive fitness values 
matrix_pos_PPI = matrix(0, sum(as.numeric(count_summary[,2])), 3)
b = 1
for(i in 1:length(mean_fitness_pos)){
  a = as.numeric(count_summary[i,2])
  index = b: (b + a-1)
  pos_env_index = which(count_summary[i,3:ncol(count_summary)] == "1")
  pos_env_index_fit = pos_env_index + 3
  matrix_pos_PPI[index, 1] = count_summary[i,1]
  matrix_pos_PPI[index, 2] = as.numeric(vScore[i, pos_env_index_fit])
  matrix_pos_PPI[index, 3] = count_summary[i,2]
  b = b + a
}

matrix_pos_PPI = data.frame(matrix_pos_PPI[,1], as.numeric(matrix_pos_PPI[,2]), matrix_pos_PPI[,3])
colnames(matrix_pos_PPI) = c("PPI", "Fitness", "Env_count")


ggplot() +
  
  geom_boxplot(aes(x = Env_count, y = Fitness), matrix_pos_PPI, outlier.shape=NA) +
  geom_dotplot(aes(x = Env_count, y = Fitness), matrix_pos_PPI, 
               binaxis="y",stackdir="center",binwidth=0.002, alpha=0.2, col = apple_colors[8]) +
  xlab("Number of environments in which a PPI is identified") +
  ylab("Fitness of a PPI") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black")) 
#theme(plot.margin = unit(c(1,1,2,1), "cm"))
ggsave("~/Desktop/Fitness_PPI_environment_number.pdf", width=5, height =5)


#######################################
#### Check the fitness of validated PPIs
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/")
d_1_1 = csvReader_T("Diploid_01_01_both_T5.csv")
d_1_2 = csvReader_T("Diploid_01_02_both_T14.csv")
d_1_3 = csvReader_T("Diploid_01_03_both_T15.csv")
d_1_4 = csvReader_T("Diploid_01_04_both_T2.csv")
d_1_5 = csvReader_T("Diploid_01_05_both_T5.csv")
d_1_6 = csvReader_T("Diploid_01_06_both_T7.csv")
d_1_7 = csvReader_T("Diploid_01_07_both_T14.csv")

d_2_1 = csvReader_T("Diploid_02_03_04_01_both_T15.csv")
d_2_2 = csvReader_T("Diploid_02_03_04_02_T7_T2.csv")
d_2_3 = csvReader_T("Diploid_02_03_04_03_both_T5.csv")
d_2_4 = csvReader_T("Diploid_02_03_04_04_both_T7.csv")
d_2_5 = csvReader_T("Diploid_02_03_04_05_both_T14.csv")
d_2_6 = csvReader_T("Diploid_02_03_04_06_both_T15.csv")
d_2_7 = csvReader_T("Diploid_02_03_04_07_both_T2.csv")
d_2_8 = csvReader_T("Diploid_02_03_04_08_both_T5.csv")
d_2_9 = csvReader_T("Diploid_02_03_04_09_both_T7.csv")
d_2_10 = csvReader_T("Diploid_02_03_04_10_both_T14.csv")

d_5_1 = csvReader_T("Diploid_05_01_both_T15.csv")
d_5_2 = csvReader_T("Diploid_05_02_both_T2.csv")
d_5_3 = csvReader_T("Diploid_05_03_both_T5.csv")
d_5_4 = csvReader_T("Diploid_05_04_both_T7.csv")
d_5_5 = csvReader_T("Diploid_05_05_both_T14.csv")
d_5_6 = csvReader_T("Diploid_05_06_both_T15.csv")
d_5_7 = csvReader_T("Diploid_05_07_both_T2.csv")

all_Tecan = rbind(d_1_1, d_1_2, d_1_3, d_1_4, d_1_5, d_1_6, d_1_7,
                  d_2_1, d_2_2, d_2_3, d_2_4, d_2_5, d_2_6, d_2_7,
                  d_2_8, d_2_9, d_2_10, d_5_1, d_5_2, d_5_3, d_5_4, 
                  d_5_5, d_5_6, d_5_7) # 695


vScore = csvReader_T("~/Dropbox/PPiseq_02/Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
count_summary = csvReader_T("~/Dropbox/PPiseq_02/Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")
SD_select = count_summary[which(count_summary[,3] == "1"),]
vScore_select = vScore[which(vScore[,1] %in% SD_select[,1]),]
vScore_select_Tecan = match_both_direction(vScore_select, all_Tecan[,1]) # 412 
vScore_select_fitness = as.numeric(vScore_select_Tecan[,4]) # 4: DMSO column

pdf("~/Desktop/Histogram_training_fitness_PPI.pdf", width =5, height = 5)
hist(vScore_select_fitness, breaks = seq(0.15, 1, by = 0.01))
dev.off()
fitness_bins = c(0.25, seq(0.26, 0.39, by = 0.01), seq(0.42, 0.66, by = 0.04), c(0.8, 1.0))

vScore_Tecan_random = vScore_select_Tecan[sample(1:nrow(vScore_select_Tecan), nrow(vScore_select_Tecan), replace = F),]
vScore_split = list()
vScore_split[[1]] = vScore_Tecan_random[1:103,]
vScore_split[[2]] = vScore_Tecan_random[104:206,]
vScore_split[[3]] = vScore_Tecan_random[207:309,]
vScore_split[[4]] = vScore_Tecan_random[310:412,]

mean_fitness_pos = rep(0, nrow(vScore))
for(i in 1:length(mean_fitness_pos)){
  pos_env_index = which(count_summary[i,3:ncol(count_summary)] == "1")
  pos_env_index_fit = pos_env_index + 3
  mean_fitness_pos[i] = mean(as.numeric(vScore[i, pos_env_index_fit]))
}
fitness_count = data.frame(mean_fitness_pos, count_summary[,2])
group_mean = rep(0, 9)
for(i in 1:9){
  index = which(as.numeric(fitness_count[,2]) == i)
  group_mean[i] = median(fitness_count[index,1])
}
count_normal = 1:9/9
### Calculate validation rate
rep_PPI_matrix = dataFrameReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/Reported_validation_matrix_SD_merge.csv")
unrep_PPI_matrix = dataFrameReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/Unreported_validation_matrix_SD_merge.csv")
merge_validate = rep_PPI_matrix[1, 2:ncol(rep_PPI_matrix)] + unrep_PPI_matrix[1, 2:ncol(rep_PPI_matrix)]
merge_nonvalidate = rep_PPI_matrix[2, 2:ncol(rep_PPI_matrix)] + unrep_PPI_matrix[2, 2:ncol(rep_PPI_matrix)]
merge_sum = merge_validate + merge_nonvalidate
merge_ratio = as.numeric(merge_validate/merge_sum)
external_test = data.frame(val_rate = merge_ratio, bins = group_mean, env_count_normal= count_normal)


split_data_generation = function(training, environment_loc, bins, all_Tecan){
  val_rate = rep(0, length(bins))
  env_count = rep(0, length(bins))
  PPI_count = rep(0, length(bins))
  for(i in 1:length(bins)){
    if(i == 1){
      fit_bin = c(0, bins[i])
    }
    else{
      fit_bin = c(bins[i-1], bins[i])
    }
    training_fitness = as.numeric(training[,environment_loc])
    index_chosen = which(training_fitness >= fit_bin[1] & training_fitness <= fit_bin[2])
    PPI_chosen = training[index_chosen, 1]
    PPI_count[i] = length(which(training_fitness >= fit_bin[1] & training_fitness <= fit_bin[2]))
    env_count[i] = mean(as.numeric(training[index_chosen, 2]))
    Tecan_select = match_both_direction(all_Tecan, PPI_chosen)
    if (length(Tecan_select) > 11){
      validate_PPI_count = length(which(as.numeric(Tecan_select[,11]) <= 0.05))
      val_rate[i] = validate_PPI_count/nrow(Tecan_select)
    }else if(length(Tecan_select) == 11){
      validate_PPI_count = as.numeric(Tecan_select[11]) <= 0.05
      val_rate[i] = validate_PPI_count/1
    }else{
      val_rate[i] = NA
    }
    
  }
  env_count = as.numeric(env_count)
  env_count_normal = env_count/9
  var_matrix = na.omit(data.frame(val_rate, bins, env_count_normal))
  return(var_matrix)
}

MSE= rep(0, 4)
MSE_external = rep(0, 4)
for (i in 1:4){
  test = vScore_split[[i]]
  remaining = which(1:4 != i)
  training = rbind(vScore_split[[remaining[1]]], vScore_split[[remaining[2]]],
                   vScore_split[[remaining[3]]])
  training_data = split_data_generation(training, 4, fitness_bins, all_Tecan)
  test_data = split_data_generation(test, 4, fitness_bins, all_Tecan)
  fitmodel = lm(val_rate ~ bins + env_count_normal, training_data)
  predict_test = predict(fitmodel, test_data, type = "response")
  MSE[i] = sum((test_data$val_rate - predict_test)^2)/length(predict_test)
  predict_external_test = predict(fitmodel, external_test, type = "response")
  MSE_external[i] =  sum((external_test$val_rate - predict_external_test)^2)/length(predict_external_test)
  
}

### Second model
MSE_2= rep(0, 4)
MSE_external_2 = rep(0, 4)
for (i in 1:4){
  test = vScore_split[[i]]
  remaining = which(1:4 != i)
  training = rbind(vScore_split[[remaining[1]]], vScore_split[[remaining[2]]],
                   vScore_split[[remaining[3]]])
  training_data = split_data_generation(training, 4, fitness_bins, all_Tecan)
  test_data = split_data_generation(test, 4, fitness_bins, all_Tecan)
  fitmodel = glm(val_rate ~ bins + env_count_normal, training_data, family = binomial)
  predict_test = predict(fitmodel, test_data, type = "response")
  MSE_2[i] = sum((test_data$val_rate - predict_test)^2)/length(predict_test)
  predict_external_test = predict(fitmodel, external_test, type = "response")
  MSE_external_2[i] =  sum((external_test$val_rate - predict_external_test)^2)/length(predict_external_test)
}

MSE # 0.04727893 0.01251794 0.02689478 0.02007871
MSE_2 # 0.04617380 0.01622315 0.02829342 0.02139248
MSE_external # 0.003877231 0.004597433 0.001143269 0.001644423
MSE_external_2 # 0.005914028 0.021486035 0.008457951 0.002842725

metrics = rbind(MSE,MSE_external, MSE_2, MSE_external_2)
metrics_1 = transform(metrics, Mean=apply(metrics,1, mean, na.rm = TRUE))
metrics_1$SD = apply(metrics,1, sd, na.rm = TRUE)
metrics_1$Group = c("Linear + test", "Linear + external", 
                    "Logistic + test", "Logistic + external")
library(ggplot2)
ggplot(metrics_1, aes(x=Group, y=Mean)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=Mean-SD, ymax=Mean+SD), width=.1,
                position=position_dodge(.9))+
  xlab("") +
  ylab("Mean square error")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, angle= 45, color = "black", hjust = 1),
        axis.text.y.left = element_text(size = 10, color = "black"))
ggsave("~/Desktop/Errors_two_model.pdf", width = 5, height =5)

##### Take linear model and training with 1,2,4
test = vScore_split[[3]]
remaining = c(1,2,4)
training = rbind(vScore_split[[remaining[1]]], vScore_split[[remaining[2]]],
                 vScore_split[[remaining[3]]])
training_data = split_data_generation(training, 4, fitness_bins, all_Tecan)
test_data = split_data_generation(test, 4, fitness_bins, all_Tecan)
fitmodel = lm(val_rate ~ bins + env_count_normal, training_data)
predict_test = predict(fitmodel, test_data, type = "response")
#MSE[i] = sum((test_data$val_rate - predict_test)^2)/length(predict_test)
predict_external_test = predict(fitmodel, external_test, type = "response")
sum((external_test$val_rate - predict_external_test)^2)/length(predict_external_test) # 0.001143269


fitmodel_2 = glm(val_rate ~ bins + env_count_normal, training_data, family = binomial)

predict_external_test_2 = predict(fitmodel_2, external_test, type = "response")
sum((external_test$val_rate - predict_external_test_2)^2)/length(predict_external_test_2) # 0.008457951


pdf("~/Desktop/Predicting_validation_rate.pdf", width =5, height =5)
plot(external_test$val_rate, predict_external_test, pch = 16, col = "blue",
     type = "p", xlim = c(0.5, 1), ylim = c(0.5,1), xlab = "Observed validation rate",
     ylab = "Predicted validation rate")
points(external_test$val_rate, predict_external_test_2, pch =16, col = "red")
lines(seq(0.5, 1, by = 0.1), seq(0.5, 1, by = 0.1), col = "black")
legend(0.5, 1, c("Linear regression", "Logistic regression"), pch = c(16, 16), 
       col = c("blue", "red"), bty = "n")
dev.off()

### Predict each PPI tested by Tecan
all_Tecan_select = all_Tecan[match(vScore_select_Tecan[,1], all_Tecan[,1]),]

Tecan_fit_select = data.frame(vScore_select_Tecan[,1], as.numeric(vScore_select_Tecan[,2])/9,
                              as.numeric(vScore_select_Tecan[,4]), as.numeric(all_Tecan_select[,11]))

colnames(Tecan_fit_select) = c("PPI", "env_count_normal", "bins", "Q_value")

predict_indiv = predict(fitmodel, Tecan_fit_select, type= "response")
Tecan_fit_select$val_rate = predict_indiv
Tecan_fit_select$label = (Tecan_fit_select$Q_value <= 0.05)
Tecan_fit_select = na.omit(Tecan_fit_select) # 398

ggplot(Tecan_fit_select,aes(x = label, y = val_rate, group = label))+
  geom_boxplot(outlier.shape=NA) +
  geom_dotplot(binaxis="y",stackdir="center",binwidth=0.005, alpha=0.2, 
               fill = "blue",col = "blue") +
  xlab("Validated by Tecan") +
  ylab("Predicted validation rate") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y.left = element_text(size = 10, color = "black")) 
#theme(plot.margin = unit(c(1,1,2,1), "cm"))
ggsave("~/Desktop/Predicting_individual_PPI.pdf", width=5, height =5)


  
  
  
  
  
  

pdf("~/Desktop/Validate_rate_versus_fitness.pdf", width = 5, height = 5)
plot(var_matrix$bins, var_matrix$val_rate, type = "p", pch = 16, ylim = c(0.5, 1), xlim = c(0.1, 1),
     xlab = "Fitness", ylab = "Validation rate")
#lines(bins, pred_rate, type = "l", col = "red")
dev.off()

pdf("~/Desktop/Validate_rate_versus_envrionments_number.pdf", width = 5, height = 5)
plot(env_count, var_matrix$val_rate, type = "p", pch = 16, ylim = c(0.5, 1), xlim = c(1,9),
     xlab = "Number of positive environments", ylab = "Validation rate")
#lines(bins, pred_rate, type = "l", col = "red")
dev.off()


# (1) fitmodel = nls(val_rate  ~ SSlogis(bins, Asym, xmid, scal))
# (2) fitmodel = lm(val_rate ~ bins + env_count_normal)

#(3) 
fitmodel = glm(val_rate ~ bins + env_count_normal, var_matrix, family = binomial)
coeff = coef(fitmodel)
# (1) pred_rate = coeff[1]/(1+exp((coeff[2] - bins)/coeff[3]))
# (2) pred_rate = coeff[1] + bins*coeff[2] + env_count_normal*coeff[3]
var_matrix_new = data.frame(bins, env_count_normal)
var_matrix_new$pred_rate = predict(fitmodel, var_matrix_new, type= "response")
#plot(env_count_normal, val_rate, type = "p")
#plot(bins, val_rate, type = "p")
 
 ## (1) First model
pdf("~/Desktop/Validate_rate_versus_fitness.pdf", width = 5, height = 5)
plot(bins, val_rate, type = "p", pch = 16, ylim = c(0.5, 1), xlim = c(0.1, 1),
     xlab = "Fitness", ylab = "Validation rate")
lines(bins, pred_rate, type = "l", col = "red")
dev.off()

## (2) Second model
pdf("~/Desktop/Validate_rate_versus_fitness_env_count.pdf", width = 5, height = 5)
plot(bins, val_rate, type = "p", pch = 16, ylim = c(0.5, 1), xlim = c(0.1, 1),
     xlab = "Fitness", ylab = "Validation rate")
lines(bins, pred_rate, type = "l", col = "red")
dev.off()

## (3) Third model
pdf("~/Desktop/Validate_rate_versus_fitness_env_count_logistic.pdf", width = 5, height = 5)
plot(var_matrix$bins, var_matrix$val_rate, type = "p", pch = 16, ylim = c(0.5, 1), xlim = c(0.1, 1),
     xlab = "Fitness", ylab = "Validation rate")
lines(var_matrix_new$bins, var_matrix_new$pred_rate, type = "l", col = "red")
dev.off()

#### Calculate the mean fitness for each group of environmental specificity and predict their validation rate
#### and then compare the predicted rate with real validation rate


vali_rate_group =coeff[1]/(1+exp((coeff[2] - group_mean)/coeff[3]))

vali_rate_group = coeff[1] + coeff[2]* group_mean + coeff[3] * (1:9)/9

count_normal = 1:9/9
vali_count = data.frame(bins = group_mean, env_count_normal= count_normal)
vali_count$rate = predict(fitmodel, vali_count, type = "response")

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/")


merge_ratio # 0.6428571 0.6666667 0.7297297 0.6896552 0.8076923 0.8412698 0.8983051 0.9756098 0.9767442
vali_rate_group # 0.8408432 0.8471980 0.8294519 0.8344757 0.8448614 0.8580337 0.8790035 0.9156009 0.9287731
vali_rate_group # 0.5875286 0.6419218 0.6989385 0.7535298 0.8075086 0.8610409 0.9131874 0.9605053 1.0108215
vali_count$rate # 0.4503239 0.5548184 0.6566244 0.7443035 0.8155703 0.8702663 0.9101497 0.9377176 0.9576256

all_Tecan_select = all_Tecan[which(all_Tecan[,1] %in% vScore[,1]),]
fitness_Tecan = as.numeric(vScore[match(all_Tecan_select[,1], vScore[,1]), 4])
SD = csvReader_T("~/Dropbox/PPiseq_02/Paper_data/PPI_mean_fitness_calling_files/SD_merge_mean_fitness_positive.csv")
SD[which(SD[,1] == "YCR021C_YHR089C"),]
