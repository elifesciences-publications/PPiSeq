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


# This script is to use the validation predicting model to generate a more accurate PPI estimation
# Here I use the SD replicates to measure if the model improved the accuracy of PPI calling

### (1) Generate the model
### (2) Predicting positive PPIs in SD and SD replicates, 
###  and finding the optimal predicted accuracy to generate the maximum reproducibility
### (3) Predict validation rate for each PPI in each environment and use the thresholds
###  to generate a "true" PPI in each environment. Then recalculate the variation score
### and check its correlation with other features.

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

d_only_1 = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/SD_only/data/Combine_TECAN/Diploid_SD_only_1_T14.csv")
d_only_2 = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/SD_only/data/Combine_TECAN/Diploid_SD_only_2_T7.csv")
d_only_3 = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/SD_only/data/Combine_TECAN/Diploid_SD_only_3_T5_second.csv")
all_Tecan = rbind(d_1_1, d_1_2, d_1_3, d_1_4, d_1_5, d_1_6, d_1_7,
                  d_2_1, d_2_2, d_2_3, d_2_4, d_2_5, d_2_6, d_2_7,
                  d_2_8, d_2_9, d_2_10, d_5_1, d_5_2, d_5_3, d_5_4, 
                  d_5_5, d_5_6, d_5_7, d_only_1, d_only_2, d_only_3) # 785 all single PPIs
csvWriter(all_Tecan, "Tecan_DMSO_positive_combine.csv")

all_Tecan = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/Tecan_DMSO_positive_combine.csv")

vScore = csvReader_T("~/Dropbox/PPiseq_02/Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
count_summary = csvReader_T("~/Dropbox/PPiseq_02/Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")
SD_select = count_summary[which(count_summary[,3] == "1"),]
vScore_select = vScore[which(vScore[,1] %in% SD_select[,1]),]
vScore_select_Tecan = match_both_direction(vScore_select, all_Tecan[,1]) # 412 
#vScore_select_fitness = as.numeric(vScore_select_Tecan[,4]) # 4: DMSO column
fitness_bins = c(0.25, seq(0.26, 0.39, by = 0.01), seq(0.42, 0.66, by = 0.04), c(0.8, 1.0))

### Here I do not need to split the data in to training and test data.
### I split the data into different bins, and obtained a dataframe of fitness, validation rate, and mean evnironment count
### Train the linear model with that data. 
### Will have two test datasets: (1) split data into different bins by stability; 
### (2) directly predict each PPI that was validated by Tecan
### In the end, predict all the PPIs in each environment

## generate test data 1 by splitting data according to their stability
SD_fit_count = data.frame(as.numeric(vScore[,4]), as.numeric(vScore[,2]))
stability_mean = rep(0, 9)
for(i in 1:9){
        index = which(as.numeric(SD_fit_count[,2]) == i)
        stability_mean[i] = mean(SD_fit_count[index,1])
}
count_normal = 1:9/9
# Calculate validation rate
rep_PPI_matrix = dataFrameReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/Reported_validation_matrix_SD_merge.csv")
unrep_PPI_matrix = dataFrameReader_T("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/Unreported_validation_matrix_SD_merge.csv")
merge_validate = rep_PPI_matrix[1, 2:ncol(rep_PPI_matrix)] + unrep_PPI_matrix[1, 2:ncol(rep_PPI_matrix)]
merge_nonvalidate = rep_PPI_matrix[2, 2:ncol(rep_PPI_matrix)] + unrep_PPI_matrix[2, 2:ncol(rep_PPI_matrix)]
merge_sum = merge_validate + merge_nonvalidate
merge_ratio = as.numeric(merge_validate/merge_sum)
external_test = data.frame(val_rate = merge_ratio, bins = stability_mean, env_count_normal= count_normal)

#### Function to prepare the training data
split_data_generation = function(training, environment_loc, bins, all_Tecan){
        val_rate = rep(0, length(bins))
        env_count = rep(0, length(bins))
        PPI_count = rep(0, length(bins))
        fit = rep(0, length(bins))
        for(i in 1:length(bins)){
                if(i == 1){
                        fit_bin = c(0, bins[i])
                }
                else{
                        fit_bin = c(bins[i-1], bins[i])
                }
                training_fitness = as.numeric(training[,environment_loc])
                index_chosen = which(training_fitness > fit_bin[1] & training_fitness <= fit_bin[2])
                PPI_chosen = training[index_chosen, 1]
                PPI_count[i] = length(which(training_fitness > fit_bin[1] & training_fitness <= fit_bin[2]))
                env_count[i] = mean(as.numeric(training[index_chosen, 2])) ## average environment number within each bin
                fit[i] = mean(as.numeric(training[index_chosen, 4])) ## average fitness within each bin
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
        var_matrix = na.omit(data.frame(val_rate, bins= fit, env_count_normal))
        return(var_matrix)
}
training_data = split_data_generation(vScore_select_Tecan, 4, fitness_bins, all_Tecan)
fitmodel = lm(val_rate ~ bins + env_count_normal, training_data)
coeffs = coef(fitmodel)
#(Intercept)             bins env_count_normal 
#0.56568422      -0.09742696       0.49596063 
#predict(fitmodel, data.frame(bins = c(0.2, 0.4), env_count_normal = c(1/9, 5/9)))
#coeffs[1] + coeffs[2]*0.2 + coeffs[3]* 1/9
predict_external_test = predict(fitmodel, external_test, type = "response")

cor(external_test$val_rate, predict_external_test, method = "spearman") # 0.9833333
cor(external_test$val_rate, predict_external_test, method = "pearson") # 0.9684343

pdf("~/Dropbox/PPiSeq_02/Working_figure/SFigures/paper/Method/SFigureSM6_Modeling_validation_rate/Predicting_validation_rate.pdf", width =5, height =5)
plot(external_test$val_rate, predict_external_test, pch = 16, col = "blue",
     type = "p", xlim = c(0.5, 1), ylim = c(0.5,1), xlab = "Observed validation rate",
     ylab = "Predicted validation rate")
#points(external_test$val_rate, predict_external_test_2, pch =16, col = "red")
lines(seq(0.5, 1, by = 0.1), seq(0.5, 1, by = 0.1), col = "black")
text(0.6, 0.95, labels = expression(paste("Spearman's ", italic(r), " = 0.98")))
dev.off()

