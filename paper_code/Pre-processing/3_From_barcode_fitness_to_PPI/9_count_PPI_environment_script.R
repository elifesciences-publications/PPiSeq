###########################
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


#### (1) Generate a count summary for each PPI (single_orientation) across different environments
#### If a PPI has two oritentations and either of them (or both) is positive in a environment, that PPI is positive in that environment
# The basic idea of the code is to attach positive environments to each PPI 

setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Positive_PPI_remove_promiscuous")

DMSO_pos = csvReader_T("SD_merge_Pos_PPI_real.csv")
H2O2_pos = csvReader_T("H2O2_Pos_PPI_real.csv")
HU_pos = csvReader_T("Hydroxyurea_Pos_PPI_real.csv")
Dox_pos = csvReader_T("Doxorubicin_Pos_PPI_real.csv")
Forskolin_pos = csvReader_T("Forskolin_Pos_PPI_real.csv")
Raffinose_pos = csvReader_T("Raffinose_Pos_PPI_real.csv")
NaCl_pos = csvReader_T("NaCl_Pos_PPI_real.csv")
cold_pos = csvReader_T("Cold_16C_Pos_PPI_real.csv")
FK506_pos = csvReader_T("FK506_Pos_PPI_real.csv")

PPI_list = list(DMSO_pos[,1],  H2O2_pos[,1], HU_pos[,1], Dox_pos[,1], Forskolin_pos[,1],
                Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1]) # store PPIs of one environment as an element in a list
all_PPI = unique(c(DMSO_pos[,1],  H2O2_pos[,1], HU_pos[,1], Dox_pos[,1], Forskolin_pos[,1],
                   Raffinose_pos[,1], NaCl_pos[,1], cold_pos[,1], FK506_pos[,1])) #13241
all_PPI_unique = mark_duplicates_fast(all_PPI) # 12552; 689 duplicates
all_PPI_matrix = matrix(0, nrow(all_PPI_unique), 10)
all_PPI_matrix[,1] = all_PPI_unique[,1]
colnames(all_PPI_matrix)= c("PPI","SD", "H2O2", "HU", "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
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
csvWriter(all_PPI_matrix_final, "PPI_environment_count_summary_SD_merge.csv")
### Remove single
PPI_count = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge.csv")
PPI_SD = PPI_count[which(PPI_count[,3] == "1"),]
SD1_single = csvReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO/counts/PPI_1_p.values.csv") # 189787
SD2_single = csvReader_T("/Volumes/zmliu_02/PPiseq_03/DMSO_2/counts/PPI_1_p.values.csv") # 168705
common_single = intersect(SD1_single[,1], SD2_single[,1]) # 100412
SD_pos_single = intersect(common_single, PPI_SD[,1]) # 783
PPI_count_filter = PPI_count[which(!PPI_count[,1] %in% SD_pos_single),]
csvWriter(PPI_count_filter, "/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge_filter.csv")

