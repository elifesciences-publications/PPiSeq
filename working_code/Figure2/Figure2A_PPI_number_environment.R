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
DMSO_only = extract_pos_PPI(DMSO_mean, "Working_data/Positive_PPI_environment/DMSO_Pos_PPI_real.csv")# 5036
DMSO_unique = mark_duplicates_fast(DMSO_only[,1]) # 4645

H2O2_mean = csvReader_T("Paper_data/H2O2_mean_fitness_positive.csv")
H2O2_only = extract_pos_PPI(H2O2_mean, "Working_data/Positive_PPI_environment/H2O2_Pos_PPI_real.csv")# 4992
H2O2_unique = mark_duplicates_fast(H2O2_only[,1]) # 4643

HU_mean = csvReader_T("Paper_data/Hydroxyurea_mean_fitness_positive.csv")
HU_only = extract_pos_PPI(HU_mean, "Working_data/Positive_PPI_environment/Hydroxyurea_Pos_PPI_real.csv")#5150
HU_unique = mark_duplicates_fast(HU_only[,1]) # 4787

Dox_mean = csvReader_T("Paper_data/Doxorubicin_mean_fitness_positive.csv")
Dox_only = extract_pos_PPI(Dox_mean, "Working_data/Positive_PPI_environment/Doxorubicin_Pos_PPI_real.csv")#3382
Dox_unique = mark_duplicates_fast(Dox_only[,1]) #3209

Forskolin_mean = csvReader_T("Paper_data/Forskolin_mean_fitness_positive.csv")
Forskolin_only = extract_pos_PPI(Forskolin_mean, "Working_data/Positive_PPI_environment/Forskolin_Pos_PPI_real.csv")#4232
Forskolin_unique = mark_duplicates_fast(Forskolin_only[,1]) # 3877

Raffinose_mean = csvReader_T("Paper_data/Raffinose_mean_fitness_positive.csv")
Raffinose_only = extract_pos_PPI(Raffinose_mean, "Working_data/Positive_PPI_environment/Raffinose_Pos_PPI_real.csv")#4371
Raffinose_unique = mark_duplicates_fast(Raffinose_only[,1]) # 4013

NaCl_mean = csvReader_T("Paper_data/NaCl_mean_fitness_positive.csv")
NaCl_only = extract_pos_PPI(NaCl_mean, "Working_data/Positive_PPI_environment/NaCl_Pos_PPI_real.csv")#2601
NaCl_unique = mark_duplicates_fast(NaCl_only[,1]) # 2481

cold_mean = csvReader_T("Paper_data/Cold_16C_mean_fitness_positive.csv")
cold_only = extract_pos_PPI(cold_mean, "Working_data/Positive_PPI_environment/Cold_16C_Pos_PPI_real.csv")#3112
cold_unique = mark_duplicates_fast(cold_only[,1]) # 3060

FK506_mean = csvReader_T("Paper_data/FK506_mean_fitness_positive.csv")
FK506_only = extract_pos_PPI(FK506_mean, "Working_data/Positive_PPI_environment/FK506_Pos_PPI_real.csv")# 3458
FK506_unique = mark_duplicates_fast(FK506_only[,1]) # 3183

#Combine all the PPIs detected in each environment
all = unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
               Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
               cold_only[,1], FK506_only[,1])) # 13050
all_dup = mark_duplicates_fast(all) #12333

# check the order of PPI number detected in different environments
PPI_count_matrix = matrix(0, 9, 2)
PPI_count_matrix[,1] = c("DMSO", "H2O2", "HU", "Dox", "Forskolin",
                         "Raffinose", "NaCl", "Cold_16C", "FK506")
PPI_count_matrix[,2] = c(nrow(DMSO_unique), nrow(H2O2_unique), nrow(HU_unique), nrow(Dox_unique),
                         nrow(Forskolin_unique), nrow(Raffinose_unique), nrow(NaCl_unique), 
                         nrow(cold_unique), nrow(FK506_unique))

PPI_count_matrix[order(as.numeric(PPI_count_matrix[,2])),]
# PPI number count order (low to high):NaCl, 16C, FK506, Dox, Forskolin, Raffinose, H2O2, DMSO, HU
# Count_order(low to high) of environment specific PPIs: 

# first generate the unique PPI sets after removing one specific environment
NaCl_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 11665

cold_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                             FK506_only[,1]))) # 10468

FK506_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1]))) # 12090

Dox_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 11265

Forskolin_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 12172

Raffinose_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 11639

H2O2_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 11476

DMSO_remove = mark_duplicates_fast(unique(c(H2O2_only[,1], HU_only[,1], Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 11726

HU_remove = mark_duplicates_fast(unique(c(DMSO_only[,1], H2O2_only[,1],  Dox_only[,1], 
                                            Forskolin_only[,1], Raffinose_only[,1], NaCl_only[,1],
                                            cold_only[,1], FK506_only[,1]))) # 11555

# PPI number count order (low to high):NaCl, 16C, FK506, Dox, Forskolin, Raffinose, H2O2, DMSO, HU
NaCl_other = match_both_direction(NaCl_unique, NaCl_remove[,1]) # 1813
NaCl_specific = NaCl_unique[which(!NaCl_unique[,1] %in% NaCl_other[,1]),] # 668

cold_other = match_both_direction(cold_unique, cold_remove[,1]) #1195
cold_specific = cold_unique[which(!cold_unique[,1] %in% cold_other[,1]),] # 1865

FK506_other = match_both_direction(FK506_unique, FK506_remove[,1]) #2940
FK506_specific = FK506_unique[which(!FK506_unique[,1] %in% FK506_other[,1]),] #243

Dox_other = match_both_direction(Dox_unique, Dox_remove[,1]) #2141
Dox_specific = Dox_unique[which(!Dox_unique[,1] %in% Dox_other[,1]),] #1068

Forskolin_other = match_both_direction(Forskolin_unique, Forskolin_remove[,1]) #3716
Forskolin_specific = Forskolin_unique[which(!Forskolin_unique[,1] %in% Forskolin_other[,1]),] #161

Raffinose_other = match_both_direction(Raffinose_unique, Raffinose_remove[,1]) #3319
Raffinose_specific = Raffinose_unique[which(!Raffinose_unique[,1] %in% Raffinose_other[,1]),] #694

H2O2_other = match_both_direction(H2O2_unique, H2O2_remove[,1]) #3786
H2O2_specific = H2O2_unique[which(!H2O2_unique[,1] %in% H2O2_other[,1]),] #857

DMSO_other = match_both_direction(DMSO_unique, DMSO_remove[,1]) #4038
DMSO_specific = DMSO_unique[which(!DMSO_unique[,1] %in% DMSO_other[,1]),] #607

HU_other = match_both_direction(HU_unique, HU_remove[,1]) #4009
HU_specific = HU_unique[which(!HU_unique[,1] %in% HU_other[,1]),] #778

## Add the environment specific PPI number into the matrix
# check the order of PPI number detected in different environments
environment_count = c(nrow(DMSO_specific), nrow(H2O2_specific), nrow(HU_specific),
                      nrow(Dox_specific), nrow(Forskolin_specific), nrow(Raffinose_specific),
                      nrow(NaCl_specific), nrow(cold_specific), nrow(FK506_specific))
PPI_count_matrix_add = cbind(PPI_count_matrix, environment_count)
colnames(PPI_count_matrix_add) = c("Environment", "PPI_number", "Environment_specific")
csvWriter(PPI_count_matrix_add, "Working_data/Positive_PPI_environment/Count_summary.csv")

PPI_count_matrix_add[order(as.numeric(PPI_count_matrix_add[,3])),1]
# "Forskolin" "FK506"     "DMSO"      "NaCl"      "Raffinose" "HU"        "H2O2"      "Dox"       "Cold_16C" 
# Put DMSO into the first position
# "DMSO", "Forskolin", "FK506", "NaCl", "Raffinose", "HU", "H2O2", "Dox", "Cold_16C"
# Make a barplot to show the cumulative number of PPIs across different environments
DMSO_Forskolin_merge = mark_duplicates_fast(unique(c(DMSO_unique[,1], Forskolin_unique[,1]))) # 5287
DMSO_Forskolin_single_merge = mark_duplicates_fast(unique(DMSO_specific[,1], Forskolin_specific[,1])) # 607

DMSO_Forskolin_FK506_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_merge[,1], FK506_unique[,1]))) #5764
DMSO_Forskolin_FK506_single_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_single_merge[,1], 
                                                                  FK506_specific[,1]))) #850

DMSO_Forskolin_FK506_NaCl_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_merge[,1], NaCl_unique[,1]))) #6520
DMSO_Forskolin_FK506_NaCl_single_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_single_merge[,1], NaCl_specific[,1])))#1518

DMSO_Forskolin_FK506_NaCl_Raffinose_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_merge[,1], Raffinose_unique[,1]))) #7352
DMSO_Forskolin_FK506_NaCl_Raffinose_single_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_single_merge[,1], 
                                                                                 Raffinose_specific[,1]))) #2212

DMSO_Forskolin_FK506_NaCl_Raffinose_HU_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_Raffinose_merge[,1], 
                                                                               HU_unique[,1]))) #8443
DMSO_Forskolin_FK506_NaCl_Raffinose_HU_single_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_Raffinose_single_merge[,1], 
                                                                                 HU_specific[,1]))) #2990

DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_merge[,1], 
                                                                               H2O2_unique[,1]))) #9348
DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_single_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_single_merge[,1], 
                                                                                      H2O2_specific[,1]))) #3847

DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_merge= mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_merge[,1], 
                                                                                   Dox_unique[,1]))) # 10468
DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_single_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_single_merge[,1], 
                                                                                          Dox_specific[,1]))) #4915

DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_cold_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_merge[,1], 
                                                                                        cold_unique[,1]))) # 12333
DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_cold_single_merge = mark_duplicates_fast(unique(c(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_single_merge[,1], 
                                                                                               cold_specific[,1]))) #6780

## "DMSO", "Forskolin", "FK506", "NaCl", "Raffinose", "HU", "H2O2", "Dox", "Cold_16C"
PPI_matrix = matrix(0, 3, 9)
PPI_matrix[1,] = c(nrow(DMSO_specific), nrow(Forskolin_specific), nrow(FK506_specific), nrow(NaCl_specific),
                   nrow(Raffinose_specific), nrow(HU_specific), nrow(H2O2_specific), nrow(Dox_specific), nrow(cold_specific))
cumulative_count = c(nrow(DMSO_unique), nrow(DMSO_Forskolin_merge), nrow(DMSO_Forskolin_FK506_merge), nrow(DMSO_Forskolin_FK506_NaCl_merge),
                     nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_merge), nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_merge),
                     nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_merge), 
                     nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_merge),
                     nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_cold_merge))
cumulative_single_count = c(nrow(DMSO_specific), nrow(DMSO_Forskolin_single_merge), nrow(DMSO_Forskolin_FK506_single_merge), 
                            nrow(DMSO_Forskolin_FK506_NaCl_single_merge),
                            nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_single_merge), 
                            nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_single_merge),
                            nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_single_merge),
                            nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_single_merge),
                            nrow(DMSO_Forskolin_FK506_NaCl_Raffinose_HU_H2O2_Dox_cold_single_merge))

PPI_matrix[2,] = as.numeric(cumulative_single_count) - as.numeric(PPI_matrix[1,])
PPI_matrix[3,] = as.numeric(cumulative_count) - as.numeric(cumulative_single_count)

colnames(PPI_matrix) = c("SD", "Forskolin", "FK506", "NaCl", "Raffinose", "Hydroxyurea",
                         "H2O2", "Doxorubicin", "16 \u00B0C")
csvWriter(PPI_matrix, "Working_data/Positive_PPI_environment/Count_summary_02.csv")

pdf("Working_figure/Figure2/Figure2A_Number_PPIs_across_environments.pdf", height = 5, width = 5)
barCenter = barplot(PPI_matrix, horiz=F, beside=F,  ylab="Cumulative number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6), yaxt="n", ylim = c(0,15000),
                    col= apple_colors[c(7,3,5)], axisnames=F, border=NA)
axis(2,at = c(seq(0, 18000, by = 3000)))
legend("topleft", legend=c("Multiple environments","Unique to another environment","Unique to this environment"), 
       fill=apple_colors[c(5,3,7)],  bty="n", border=FALSE)

text(x= barCenter, y = -500, labels = c("SD", "Forskolin", "FK506", "NaCl", "Raffinose", "Hydroxyurea",  
                                        expression('H'[2]* 'O'[2]), "Doxorubicin", "16 \u00B0C"), 
     srt= 45, adj = 1, xpd = TRUE)
dev.off()





