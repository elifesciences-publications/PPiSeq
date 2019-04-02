###########################
#### This script is use to generate Figure 2, which show the identified PPI number
#### across multiple environments, PPI network, validation rate, etc.
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

# Figure 2A to show the cumulative PPI number after checking more environments
setwd("~/Dropbox/PPiSeq_02/")
#write a function to subset positive PPIs from mean-fitness file and remove control PPIs
# Output the files into a folder of "Pos_PPI_environment"
extract_pos_PPI = function(DMSO_mean, output_name){
        pos_PPI = DMSO_mean[which(DMSO_mean[,7] == "1"),]
        PPI_RRS = pos_PPI[grep("Neg_PPI", pos_PPI[,1]),1] #3
        PPI_PRS = pos_PPI[grep("Pos_PPI", pos_PPI[,1]),1] #29
        DHFR_pos = pos_PPI[grep("positive_DHFR", pos_PPI[,1]),1] # 1
        DHFR_neg = pos_PPI[grep("negative_non_DHFR", pos_PPI[,1]),1] # 1
        PPI_control = c(PPI_PRS, PPI_RRS, DHFR_pos, DHFR_neg)
        pos_PPI_real = pos_PPI[which(!pos_PPI[,1] %in% PPI_control),] # Remove positive controls strains
        csvWriter(pos_PPI_real, output_name)
        return(pos_PPI_real)
}
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv")
DMSO_only = extract_pos_PPI(DMSO_mean, "Working_data/Positive_PPI_environment/DMSO_Pos_PPI_real.csv")# 5178
DMSO_unique = mark_duplicates_fast(DMSO_only[,1]) # 4786

H2O2_mean = csvReader_T("Paper_data/H2O2_mean_fitness_positive.csv")
H2O2_only = extract_pos_PPI(H2O2_mean, "Working_data/Positive_PPI_environment/H2O2_Pos_PPI_real.csv")# 4999
H2O2_unique = mark_duplicates_fast(H2O2_only[,1]) # 4650

HU_mean = csvReader_T("Paper_data/Hydroxyurea_mean_fitness_positive.csv")
HU_only = extract_pos_PPI(HU_mean, "Working_data/Positive_PPI_environment/Hydroxyurea_Pos_PPI_real.csv")#6299
HU_unique = mark_duplicates_fast(HU_only[,1]) # 5936

Dox_mean = csvReader_T("Paper_data/Doxorubicin_mean_fitness_positive.csv")
Dox_only = extract_pos_PPI(Dox_mean, "Working_data/Positive_PPI_environment/Doxorubicin_Pos_PPI_real.csv")#3492
Dox_unique = mark_duplicates_fast(Dox_only[,1]) # 3315

Forskolin_mean = csvReader_T("Paper_data/Forskolin_mean_fitness_positive.csv")
Forskolin_only = extract_pos_PPI(Forskolin_mean, "Working_data/Positive_PPI_environment/Forskolin_Pos_PPI_real.csv")#4257
Forskolin_unique = mark_duplicates_fast(Forskolin_only[,1]) # 3901

Raffinose_mean = csvReader_T("Paper_data/Raffinose_mean_fitness_positive.csv")
Raffinose_only = extract_pos_PPI(Raffinose_mean, "Working_data/Positive_PPI_environment/Raffinose_Pos_PPI_real.csv")#4378
Raffinose_unique = mark_duplicates_fast(Raffinose_only[,1]) # 4019

NaCl_mean = csvReader_T("Paper_data/NaCl_mean_fitness_positive.csv")
NaCl_only = extract_pos_PPI(NaCl_mean, "Working_data/Positive_PPI_environment/NaCl_Pos_PPI_real.csv")#2725
NaCl_unique = mark_duplicates_fast(NaCl_only[,1]) # 2592

cold_mean = csvReader_T("Paper_data/Cold_16C_mean_fitness_positive.csv")
cold_only = extract_pos_PPI(cold_mean, "Working_data/Positive_PPI_environment/Cold_16C_Pos_PPI_real.csv")#3154
cold_unique = mark_duplicates_fast(cold_only[,1]) # 3102

FK506_mean = csvReader_T("Paper_data/FK506_mean_fitness_positive.csv")
FK506_only = extract_pos_PPI(FK506_mean, "Working_data/Positive_PPI_environment/FK506_Pos_PPI_real.csv")# 3472
FK506_unique = mark_duplicates_fast(FK506_only[,1]) # 3196

#Combine all the PPIs detected in each environment
all = unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
               Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
               cold_only[,1], FK506_only[,1])) # 14564
all_dup = mark_duplicates_fast(all) #13829

# check the order of PPI number detected in different environments
PPI_count_matrix = matrix(0, 9, 2)
PPI_count_matrix[,1] = c("DMSO", "H2O2", "HU", "Dox", "Forskolin",
                         "Raffinose", "NaCl", "Cold_16C", "FK506")
PPI_count_matrix[,2] = c(nrow(DMSO_unique), nrow(H2O2_unique), nrow(HU_unique), nrow(Dox_unique),
                         nrow(Forskolin_unique), nrow(Raffinose_unique), nrow(NaCl_unique), 
                         nrow(cold_unique), nrow(FK506_unique))

# PPI number count order (low to high):NaCl, 16C, FK506, Dox, Forskolin, Raffinose, H2O2, DMSO, HU
# Count_order(low to high) of environment specific PPIs: 

# first generate the unique PPI sets after removing one specific environment
NaCl_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 13094

cold_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                             FK506_only[,1]))) # 11943

FK506_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1]))) # 13580

Dox_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 12697

Forskolin_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 13658

Raffinose_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 13135

H2O2_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 12967

DMSO_remove = mark_duplicates_fast(unique(c(H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 13106

HU_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1],  Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 11919

# PPI number count order (low to high):NaCl, 16C, FK506, Dox, Forskolin, Raffinose, H2O2, DMSO, HU
NaCl_other = match_both_direction(NaCl_unique, NaCl_remove[,1]) # 1857
NaCl_specific = NaCl_unique[which(!NaCl_unique[,1] %in% NaCl_other[,1]),] # 735

cold_other = match_both_direction(cold_unique, cold_remove[,1]) #1216
cold_specific = cold_unique[which(!cold_unique[,1] %in% cold_other[,1]),] # 1886

FK506_other = match_both_direction(FK506_unique, FK506_remove[,1]) #2947
FK506_specific = FK506_unique[which(!FK506_unique[,1] %in% FK506_other[,1]),] #249

Dox_other = match_both_direction(Dox_unique, Dox_remove[,1]) #2183
Dox_specific = Dox_unique[which(!Dox_unique[,1] %in% Dox_other[,1]),] #1132

Forskolin_other = match_both_direction(Forskolin_unique, Forskolin_remove[,1]) #3730
Forskolin_specific = Forskolin_unique[which(!Forskolin_unique[,1] %in% Forskolin_other[,1]),] #171

Raffinose_other = match_both_direction(Raffinose_unique, Raffinose_remove[,1]) #3325
Raffinose_specific = Raffinose_unique[which(!Raffinose_unique[,1] %in% Raffinose_other[,1]),] #694

H2O2_other = match_both_direction(H2O2_unique, H2O2_remove[,1]) #3788
H2O2_specific = H2O2_unique[which(!H2O2_unique[,1] %in% H2O2_other[,1]),] #862

DMSO_other = match_both_direction(DMSO_unique, DMSO_remove[,1]) #4063
nrow(DMSO_other)
DMSO_specific = DMSO_unique[which(!DMSO_unique[,1] %in% DMSO_other[,1]),] #723
nrow(DMSO_specific)

HU_other = match_both_direction(HU_unique, HU_remove[,1]) #4026

HU_specific = HU_unique[which(!HU_unique[,1] %in% HU_other[,1]),] #1910

## Add the environment specific PPI number into the matrix
# check the order of PPI number detected in different environments
environment_count = c(nrow(DMSO_specific), nrow(H2O2_specific), nrow(HU_specific),
                      nrow(Dox_specific), nrow(Forskolin_specific), nrow(Raffinose_specific),
                      nrow(NaCl_specific), nrow(cold_specific), nrow(FK506_specific))
PPI_count_matrix_add = cbind(PPI_count_matrix, environment_count)
colnames(PPI_count_matrix_add) = c("Environment", "PPI_number", "Environment_specific")
csvWriter(PPI_count_matrix_add, "Working_data/Positive_PPI_environment/Count_summary.csv")

PPI_count_matrix_add[order(as.numeric(PPI_count_matrix_add[,3])),1]
# "Forskolin" "FK506"     "Raffinose" "DMSO"      "NaCl"      "H2O2"      "Dox"       "Cold_16C"  "HU"
# Make a barplot to show the cumulative number of PPIs across different environments
Forskolin_FK506_merge = mark_duplicates_fast(unique(c(Forskolin_unique[,1], FK506_unique[,1])))
Forskolin_FK506_Raffinose_merge = mark_duplicates_fast(unique(c(Forskolin_FK506_merge[,1], Raffinose_unique[,1])))
Forskolin_FK506_Raffinose_DMSO_merge = mark_duplicates_fast(unique(c(Forskolin_FK506_Raffinose_merge[,1], DMSO_unique[,1])))
Forskolin_FK506_Raffinose_DMSO_NaCl_merge = mark_duplicates_fast(unique(c(Forskolin_FK506_Raffinose_DMSO_merge[,1], NaCl_unique[,1])))
Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_merge = mark_duplicates_fast(unique(c(Forskolin_FK506_Raffinose_DMSO_NaCl_merge[,1], 
                                                                               H2O2_unique[,1])))
Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_Dox_merge = 
        mark_duplicates_fast(unique(c(Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_merge[,1], Dox_unique[,1])))
Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_Dox_Cold_merge = 
        mark_duplicates_fast(unique(c(Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_Dox_merge[,1], cold_unique[,1])))

Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_Dox_Cold_HU_merge = 
        mark_duplicates_fast(unique(c(Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_Dox_Cold_merge[,1], HU_unique[,1])))

PPI_matrix = matrix(0, 2, 9)
PPI_matrix[1,] = c(nrow(Forskolin_specific), nrow(FK506_specific), nrow(Raffinose_specific),
                   nrow(DMSO_specific), nrow(NaCl_specific), nrow(H2O2_specific),
                   nrow(Dox_specific), nrow(cold_specific), nrow(HU_specific))
cumulative_count = c(nrow(Forskolin_unique), nrow(Forskolin_FK506_merge), nrow(Forskolin_FK506_Raffinose_merge),
                     nrow(Forskolin_FK506_Raffinose_DMSO_merge), nrow(Forskolin_FK506_Raffinose_DMSO_NaCl_merge),
                     nrow(Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_merge), nrow(Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_Dox_merge),
                     nrow(Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_Dox_Cold_merge),
                     nrow(Forskolin_FK506_Raffinose_DMSO_NaCl_H2O2_Dox_Cold_HU_merge))

PPI_matrix[2,] = as.numeric(cumulative_count) - as.numeric(PPI_matrix[2,])
colnames(PPI_matrix) = c("Forskolin", "FK506", "Raffinose", "DMSO", "NaCl", 
                   "H2O2", "Doxorubicin", "16C", "Hydroxyurea")
cumulative_count # 3901  4547  5695  6779  7623  8633  9940 11919 13829
PPI_matrix[1,]  # 171   249    694   723   735   862   1132  1886  1910 (Environment sepcific PPIs)
c(nrow(Forskolin_unique), nrow(FK506_unique), nrow(Raffinose_unique), nrow(DMSO_unique),
  nrow(NaCl_unique), nrow(H2O2_unique), nrow(Dox_unique), nrow(cold_unique), nrow(HU_unique))
# 3901 3196 4019 4786 2592 4650 3315 3102 5936 # number of PPIs in each environment

pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2A_Number_PPIs_across_environments.pdf", height = 5, width = 5)
barCenter = barplot(PPI_matrix, horiz=F, beside=F, ylim=c(0,20000), ylab="Cumulative number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
                    col= apple_colors[c(1,2)], axisnames=F, border=NA, cex.axis=0.8)
legend("topleft", legend=c("Single environment", "Multiple environments"), 
       fill=apple_colors[c(1,2)], cex=0.8, bty="n", border=FALSE)
#text(x= barCenter, y = as.numeric(cumulative_count) + 1200, 
     #labels = as.character(cumulative_count) , cex=0.8, xpd = TRUE) # add cumulative number

#text(x= barCenter-0.3, y = -500, 
     #labels = c("171/3901", "249/3196", "694/4019", "723/4786","735/2529", "862/4650", "1132/3315", "1886/3102", "1910/5936"), 
     #cex=0.6, xpd = TRUE, srt= 45, adj =1) # add environment specific PPI number/ PPIs detected in that environment
text(x= barCenter, y = -500, labels = colnames(PPI_matrix), srt= 45, adj = 1, cex=0.8, xpd = TRUE)
dev.off()

##################################################
# Figure2B make a pie plot to show the number of detected environments for each PPI 
# The idea is attach positive environments to each PPI 
setwd("~/Dropbox/PPiSeq_02/")
DMSO_pos = csvReader_T("Working_data/Positive_PPI_environment/DMSO_Pos_PPI_real.csv")
H2O2_pos = csvReader_T("Working_data/Positive_PPI_environment/H2O2_Pos_PPI_real.csv")
HU_pos = csvReader_T("Working_data/Positive_PPI_environment/Hydroxyurea_Pos_PPI_real.csv")
Dox_pos = csvReader_T("Working_data/Positive_PPI_environment/Doxorubicin_Pos_PPI_real.csv")
Forskolin_pos = csvReader_T("Working_data/Positive_PPI_environment/Forskolin_Pos_PPI_real.csv")
Raffinose_pos = csvReader_T("Working_data/Positive_PPI_environment/Raffinose_Pos_PPI_real.csv")
NaCl_pos = csvReader_T("Working_data/Positive_PPI_environment/NaCl_Pos_PPI_real.csv")
cold_pos = csvReader_T("Working_data/Positive_PPI_environment/Cold_16C_Pos_PPI_real.csv")
FK506_pos = csvReader_T("Working_data/Positive_PPI_environment/FK506_Pos_PPI_real.csv")
PPI_list = list(DMSO_pos[,1], H2O2_pos[,1], HU_pos[,1], Dox_pos[,1], Forskolin_pos[,1],
                Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1]) # store PPIs of one environment as an element in a list
all_PPI = unique(c(DMSO_pos[,1], H2O2_pos[,1], HU_pos[,1], Dox_pos[,1], Forskolin_pos[,1],
                   Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1])) #14564
all_PPI_unique = mark_duplicates_fast(all_PPI) # 13829
all_PPI_matrix = matrix(0, nrow(all_PPI_unique), 10)
all_PPI_matrix[,1] = all_PPI_unique[,1]
colnames(all_PPI_matrix)= c("PPI","DMSO", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
# Check each PPI in each environment, if reported, put a value of 1 into the space.
for(i in 1:nrow(all_PPI_matrix)){
        if(all_PPI_unique[i,2] != "0"){
                for(j in 1:length(PPI_list)){
                        if(all_PPI_unique[i,1] %in% PPI_list[[j]] | all_PPI_unique[i,2] %in% PPI_list[[j]]){
                                all_PPI_matrix[i,j +1] = 1
                        }
                }    
        }else {
                for(j in 1:length(PPI_list)){
                        if(all_PPI_unique[i,1] %in% PPI_list[[j]]){
                                all_PPI_matrix[i,j +1] = 1
                        }
                }    
        }
        
}

environment_number = rep(0, nrow(all_PPI_matrix)) # adding the column of 2:9 in each row will give the number of positive environment
for(i in 1:length(environment_number)){
        environment_number[i] = sum(as.numeric(all_PPI_matrix[i,2:10]))
}
all_PPI_matrix_final = cbind(all_PPI_matrix[,1], environment_number, all_PPI_matrix[,2:10])
csvWriter(all_PPI_matrix_final, "Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")

## make a pie plot for environment_number
environment_count = data.frame(table(environment_number))
## Environment number per positive PPI

library(ggplot2)
colfunc <- colorRampPalette(apple_colors[c(5,3,7)])
ggplot(environment_count, aes(x = factor(1), y = Freq, fill= environment_number)) +
        geom_bar(width =1, size = 0.1, color = "white", stat = "identity") + 
        geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), size = 1.8)+
        coord_polar("y") +
        theme_minimal() +
        theme(axis.text = element_blank(),
              panel.border = element_blank(),
              panel.grid=element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        scale_fill_manual(values=colfunc(9),
                          name = "Number of environments",
                          breaks = 1:9,
                          labels = 1:9)

ggsave("~/Dropbox/PPiSeq_02/Working_figure/Figure2B_PPI_environment_distribution.pdf", width =5 , height = 5)

# Or make a barplot to show how many of them have been reproted
reported_PPI = csvReader_T("Working_data/multiple_validated_PPI.csv")
matrix_PPI_env_rep = matrix(0, 2, 9)
for(i in 1:9){
        all = all_PPI_matrix_final[which(all_PPI_matrix_final[,2] == i),]
        all_reported = match_both_direction(all,reported_PPI[,1])
        all_unreported = all[which(!all[,1] %in% all_reported[,1]),]
        matrix_PPI_env_rep[1,i] = nrow(all_reported)
        matrix_PPI_env_rep[2,i] = nrow(all_unreported)
}
matrix_PPI_env_rep[1,] # 308 144 112 115 148 212 254 338 222
all_PPI_count = matrix_PPI_env_rep[1,] + matrix_PPI_env_rep[2,]
all_PPI_count # 8362 1356  655  535  555  664  570  619  513
ratio = matrix_PPI_env_rep[1,]/all_PPI_count
ratio # 0.03683329 0.10619469 0.17099237 0.21495327 0.26666667 0.31927711 0.44561404 0.54604200 0.43274854
ratio_reported = c("3.7%", "10.6%", "17%", "21.5%", "26.7%", "31.9%", "44.6%", "54.6%", "43.3%")
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2B_02_Number_environments_PPI_reproted.pdf", height = 5, width = 5)
barCenter = barplot(matrix_PPI_env_rep, horiz=F, beside=F, ylim=c(0,10000), ylab="Number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
                    col= apple_colors[c(6,7)], axisnames=F, border=NA, cex.axis=0.8)
legend("topright", legend=c("Reported", "Unreported"), 
       fill=apple_colors[c(6,7)], bty="n", border=FALSE)
text(x= barCenter, y = all_PPI_count + 300, labels = ratio_reported , 
     cex=0.6, xpd = TRUE, col= "darkgreen") # add cumulative number
text(x= barCenter, y = -500, labels = as.character(1:9), xpd = TRUE)
text(median(barCenter), y = -1200, labels = "Number of positive environments", xpd = TRUE)
dev.off()

#######################################
## Figure2C A barplot to show the validation rate
setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/pos_PPI/Combine_TECAN/")
### Combine all the chosen PPIs 
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
                       d_5_5, d_5_6, d_5_7)
reported_PPI = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/multiple_validated_PPI.csv")
all_Tecan_rep = match_both_direction(all_Tecan, reported_PPI[,1]) # 185
all_Tecan_unrep = all_Tecan[which(!all_Tecan[,1] %in% all_Tecan_rep[,1]),] # 510
#split validated PPIs into different groups
PPI_group = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
length(which(all_Tecan[,1] %in% PPI_group[,1])) # 695 all in

rep_PPI_matrix = matrix(0, 4,9)
unrep_PPI_matrix = matrix(0, 4, 9)

for(i in 1:9){
        PPI_select = PPI_group[which(as.numeric(PPI_group[,2]) == i),1]
        reported_PPI_select = all_Tecan_rep[which(all_Tecan_rep[,1] %in% PPI_select),]
        validate_PPI = length(which(as.numeric(reported_PPI_select[,11]) <= 0.05))
        non_val_PPI = length(which(as.numeric(reported_PPI_select[,11]) > 0.05))
        rep_PPI_matrix[1,i] = validate_PPI
        rep_PPI_matrix[2,i] = non_val_PPI
        rep_PPI_matrix[3,i] = validate_PPI + non_val_PPI
        rep_PPI_matrix[4,i] = rep_PPI_matrix[1,i]/rep_PPI_matrix[3,i]
        
        unrep_PPI_select = all_Tecan_unrep[which(all_Tecan_unrep[,1] %in% PPI_select),]
        validate_PPI_unrep = length(which(as.numeric(unrep_PPI_select[,11]) <= 0.05))
        non_val_PPI_unrep = length(which(as.numeric(unrep_PPI_select[,11]) > 0.05))
        unrep_PPI_matrix[1,i] = validate_PPI_unrep
        unrep_PPI_matrix[2,i] = non_val_PPI_unrep
        unrep_PPI_matrix[3,i] = validate_PPI_unrep + non_val_PPI_unrep
        unrep_PPI_matrix[4,i] = unrep_PPI_matrix[1,i]/unrep_PPI_matrix[3,i]
}
colnames(rep_PPI_matrix) = c("Envir_1", "Envir_2", "Envir_3", "Envir_4", "Envir_5",
                             "Envir_6", "Envir_7", "Envir_8", "Envir_9")
rownames(rep_PPI_matrix) = c("Validate", "Nonvalidate", "Total", "Val_ratio")
csvWriter_rownames(rep_PPI_matrix, "Reported_validation_matrix.csv")

colnames(unrep_PPI_matrix) = c("Envir_1", "Envir_2", "Envir_3", "Envir_4", "Envir_5",
                             "Envir_6", "Envir_7", "Envir_8", "Envir_9")
rownames(unrep_PPI_matrix) = c("Validate", "Nonvalidate", "Total", "Val_ratio")
csvWriter_rownames(unrep_PPI_matrix, "Unreported_validation_matrix.csv")

### Make barplot to show the percentage
### put the reported and unreported on to the same figure
ratio_rep = rep_PPI_matrix[4,]
ratio_unrep = unrep_PPI_matrix[4,]
ratio_all = c(ratio_rep[1], ratio_unrep[1], ratio_rep[2], ratio_unrep[2], 
              ratio_rep[3], ratio_unrep[3], ratio_rep[4], ratio_unrep[4],
              ratio_rep[5], ratio_unrep[5], ratio_rep[6], ratio_unrep[6],
              ratio_rep[7], ratio_unrep[7], ratio_rep[8], ratio_unrep[8],
              ratio_rep[9], ratio_unrep[9])
rep_PPI_matrix[1,] #2       5      19      17      10      14      33      34      22 
rep_PPI_matrix[3,] #7      10      21      20      10      21      38      36      22 
unrep_PPI_matrix[1,]#69      44      43      37      37      33      27      42      35 
unrep_PPI_matrix[3,]# 132      89      57      41      41      38      34      42      36 

counts_label = c("2/7", "69/132", "5/10", "44/89", "19/21", "43/57",
                 "17/20", "37/41", "10/10", "37/41", "14/21", "33/38",
                 "33/38", "27/34", "34/36", "42/42", "22/22", "35/36")
pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure2C_Validation_bar_plot.pdf", width= 6, height=5)
barCenter = barplot(ratio_all*100, horiz=F, beside=F, ylim=c(0,100), ylab="Validation rate (%)",
                    space= c(0.4, 0.08, 0.4, 0.08, 0.4, 0.08, 0.4, 0.08, 0.4, 0.08,
                             0.4, 0.08, 0.4, 0.08, 0.4, 0.08, 0.4, 0.08),
                    col= apple_colors[6:7] , axisnames=F, border=NA, cex.axis=0.8)
legend("topleft", legend=c("Reported", "Unreported"),fill=apple_colors[6:7], cex=0.8, bty="n", border=FALSE)
text(x= barCenter, y = ratio_all*100 + 2, labels = counts_label, cex=0.5, xpd = TRUE)
env_num_loc = rep(0, 9)
for(i in 1:9){
        env_num_loc[i] = mean(barCenter[(2*i-1):(2*i)])
}
text(x = env_num_loc, y = -8, labels = as.character(1:9), cex = 0.8, xpd = TRUE)
text(median(barCenter), y = -16, labels = "Number of positive environments", xpd = TRUE)
dev.off()

