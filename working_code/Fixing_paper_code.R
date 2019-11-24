
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
setwd("~/Dropbox/PPiSeq_02/")
### This script is to double check numbers mentioned in manuscript

source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

### (1) How many protein-protein pairs we have tested (1.6 million), and how many replicates
setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/")
lineage = csvReader_T("Normalized_multiple_files/All_PPI_environments_normalized_fit_SD_merge.csv")
nrow(na.omit(lineage))
#DMSO_matrix = csvReader_T("Normalized_multiple_files/SD_merge_normalized_multiple.csv")
DMSO_1 = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Normalized_multiple_files/SD_normalized_multiple.csv")
DMSO_2 = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments/Normalized_multiple_files/SD2_normalized_multiple.csv")
Forskolin_norm = csvReader_T("Normalized_multiple_files/Forskolin_normalized_multiple.csv")
FK506_norm = csvReader_T("Normalized_multiple_files/FK506_normalized_multiple.csv")
Raffinose_norm = csvReader_T("Normalized_multiple_files/Raffinose_normalized_multiple.csv")
NaCl_norm = csvReader_T("Normalized_multiple_files/NaCl_normalized_multiple.csv")
H2O2_norm = csvReader_T("Normalized_multiple_files/H2O2_normalized_multiple.csv")
Dox_norm = csvReader_T("Normalized_multiple_files/Doxorubicin_normalized_multiple.csv")
cold_norm = csvReader_T("Normalized_multiple_files/Cold_16C_normalized_multiple.csv")
HU_norm = csvReader_T("Normalized_multiple_files/Hydroxyurea_normalized_multiple.csv")

## Check these one-barcode marked PPIs
DMSO_1_single = csvReader_T("PPI_single_files/SD_PPI_1_p.values.csv")
DMSO_2_single = csvReader_T("PPI_single_files/SD2_PPI_1_p.values.csv")
Forskolin_single = csvReader_T("PPI_single_files/Forskolin_PPI_1_p.values.csv")
FK506_single = csvReader_T("PPI_single_files/FK506_PPI_1_p.values.csv")
Raffinose_single = csvReader_T("PPI_single_files/Raffinose_PPI_1_p.values.csv")
NaCl_single = csvReader_T("PPI_single_files/NaCl_PPI_1_p.values.csv")
H2O2_single = csvReader_T("PPI_single_files/H2O2_PPI_1_p.values.csv")
Dox_single = csvReader_T("PPI_single_files/Doxorubicin_PPI_1_p.values.csv")
cold_single = csvReader_T("PPI_single_files/Cold_16C_PPI_1_p.values.csv")
HU_single = csvReader_T("PPI_single_files/Hydroxyurea_PPI_1_p.values.csv")


shared_PPI_multiple = Reduce(intersect, list(DMSO_1[,1],DMSO_2[,1],Forskolin_norm[,1], FK506_norm[,1],
                                    Raffinose_norm[,1], NaCl_norm[,1], H2O2_norm[,1],
                                    Dox_norm[,1], cold_norm[,1], HU_norm[,1])) 
length(shared_PPI_multiple) # 393470

shared_PPI_single = Reduce(intersect, list(DMSO_1_single[,1],DMSO_2_single[,1],Forskolin_single[,1], 
                                           FK506_single[,1], Raffinose_single[,1], NaCl_single[,1], 
                                           H2O2_single[,1], Dox_single[,1], cold_single[,1], 
                                           HU_single[,1])) # 393470
length(shared_PPI_single) # 2582


all_PPI_multiple = unique(c(DMSO_1[,1],DMSO_2[,1],Forskolin_norm[,1], FK506_norm[,1],
                   Raffinose_norm[,1], NaCl_norm[,1], H2O2_norm[,1],
                   Dox_norm[,1], cold_norm[,1], HU_norm[,1])) 
length(all_PPI_multiple) # 1593779 Use this number

all_PPI_single = unique(c(DMSO_1_single[,1],DMSO_2_single[,1],Forskolin_single[,1], 
                          FK506_single[,1], Raffinose_single[,1], NaCl_single[,1], 
                          H2O2_single[,1], Dox_single[,1], cold_single[,1], 
                          HU_single[,1])) 
length(all_PPI_single) # 1027646

all_PPI= unique(c(all_PPI_multiple, all_PPI_single))
length(all_PPI) # 1750542


measurements_multiple = sum(as.numeric(c(DMSO_1[,2],DMSO_2[,2],Forskolin_norm[,2], FK506_norm[,2],
                                      Raffinose_norm[,2], NaCl_norm[,2], H2O2_norm[,2],
                                      Dox_norm[,2], cold_norm[,2], HU_norm[,2])))
measurements_multiple # 42,197,816 Use this number

measurements_single = sum(as.numeric(c(DMSO_1_single[,2],DMSO_2_single[,2],Forskolin_single[,2], 
                                       FK506_single[,2], Raffinose_single[,2], NaCl_single[,2], 
                                       H2O2_single[,2], Dox_single[,2], cold_single[,2], 
                                       HU_single[,2])))
measurements_single # 2334029

total_measurements = measurements_multiple + measurements_single
total_measurements # 44531845

### (2) How many PPIs we have found and how many of them are environment specific
setwd("~/Dropbox/PPiseq_02/")
PPI_count = csvReader_T("Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge.csv")
nrow(PPI_count) # 13764
# Ratio of newly finding PPIs
reported_PPI = 15026 # (PCA + BioGRID)
overlap_reported = 1915
PPiseq_new = 13764 -1915
PPiseq_new/reported_PPI # 78.8%

# Ratio of environment specific PPIs
PPI_count_filter = csvReader_T("Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")
length(which(as.numeric(PPI_count_filter[,2]) <= 3)) # 9720
nrow(PPI_count_filter) #12981
unstable_ratio = 9720/nrow(PPI_count_filter)
unstable_ratio

nrow(PPI_count_filter) # 12981
length(which(PPI_count_filter[,2] =="1")) # 7724
environment_dependent_ratio = 7724/12981
environment_dependent_ratio # 0.5950235


### (3) Negative control group in t.test in SD media
setwd("/Volumes/zmliu_02/PPiseq_03/DMSO/counts/")
PPI_fitness = csvReader_T("known_PPI_fitness_barcodes_sorted.csv")
MATa_DHFR12 = csvReader_T("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/PPI_barcodes/MATa_genome_combine.csv")
MATalpha_DHFR3 = csvReader_T("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/PPI_barcodes/MATalpha_genome_combine.csv")
PPI_DHFR12 = PPI_fitness[which(PPI_fitness[,3] %in% MATa_DHFR12[,3]),]
PPI_DHFR3 = PPI_fitness[which(PPI_fitness[,3] %in% MATalpha_DHFR3[,3]),]
#PPI_negative_DHFR= PPI_indiv[which(PPI_indiv[,1] == "negative_non_DHFR"),] 
#PPI_negative = rbind(PPI_DHFR12, PPI_DHFR3, PPI_negative_DHFR)
PPI_negative = rbind(PPI_DHFR12, PPI_DHFR3) 
PPI_negative_fitness= as.numeric(PPI_negative[,4])
length(PPI_negative_fitness)#17483
length(unique(PPI_negative[,1])) # 6071

### (4) Number of PPIs covered by marginal PCA and in these PPIs how many are reported by PPiseq
setwd("~/Dropbox/PPiseq_02/")
PCA_lower = as.matrix(read.table("Paper_data/Outside_datasets/Tarassov_PPI_PPV_80.txt", header= T, sep = "\t")) # 10230
PCA_lower = PCA_lower[which(as.numeric(PCA_lower[,9]) >= 80),] # 10230
PCA_lower_PPI = paste(PCA_lower[,1], PCA_lower[,4], sep = "_")
PCA_paper = read.delim(file = "~/Desktop/Big_PPiSeq_paper/PCA_rearray/PPI_set_science.txt", sep=" ")
PCA_called_PPI = paste(PCA_paper[,1], PCA_paper[,3], sep="_") # 2770
PCA_called_PPI_reverse = paste(PCA_paper[,3], PCA_paper[,1], sep="_")
length(intersect(PCA_lower_PPI, PCA_called_PPI)) # 1772
length(intersect(PCA_lower_PPI, PCA_called_PPI_reverse)) # 1660
PCA_called_PPI_matrix = cbind(PCA_called_PPI, PCA_called_PPI_reverse)
nrow(match_both_direction(PCA_called_PPI_matrix, PCA_lower_PPI))# 2770
PCA_called_PPI_both = unique(c(PCA_called_PPI, PCA_called_PPI_reverse))
PCA_lower_PPI_not_calling = PCA_lower_PPI[which(!PCA_lower_PPI %in% PCA_called_PPI_both)] # 7034
PCA_lower_PPI_not_calling_unique = mark_duplicates_fast(PCA_lower_PPI_not_calling) # 6876

fragment_select = csvReader_T("Paper_data/Useful_datasets/Promiscuous_protein_summary_SD_merge.csv")
fragment_protein = fragment_select[which(fragment_select[,2] >= 2),1]

## Remove PPIs that contain promiscuous proteins from BioGRID
promiscuous_PCA_lower_select = rep(0, length(PCA_lower_PPI_not_calling_unique[,1]))
for(i in 1: length(promiscuous_PCA_lower_select)){
        PPI = split_string(PCA_lower_PPI_not_calling_unique[i])
        if(PPI[1] %in% fragment_protein | PPI[2] %in% fragment_protein){
                promiscuous_PCA_lower_select[i] = 1
        }
}
length(which(promiscuous_PCA_lower_select == 1))
PCA_lower_PPI_not_calling_unique_filter = PCA_lower_PPI_not_calling_unique[which(promiscuous_PCA_lower_select != 1),] # 6786

all_PPI = csvReader_T("Paper_data/Useful_datasets/All_PPI_environments_normalized_fit_SD_merge.csv") 
pos_PPI = csvReader_T("Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge.csv")

### Marginal PPIs covered by PPiSeq
nrow(match_both_direction(PCA_lower_PPI_not_calling_unique_filter, all_PPI[,1])) # 5347
nrow(match_both_direction(PCA_lower_PPI_not_calling_unique_filter, pos_PPI[,1])) # 1838


### (5) Checking if one environment PPI outnumbered other PPIs
PPI_count_filter = csvReader_T("Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")
count_freq = as.data.frame(table(PPI_count_filter[,2]))
count_freq[1,] # 7724
sum(as.numeric(count_freq[2:nrow(count_freq),2])) # 5257