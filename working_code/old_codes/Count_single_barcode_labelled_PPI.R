################# Count the PPIs that were marked by one barcode

setwd("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/PPI_single_files")
DMSO_single = csvReader_T("SD_PPI_1_p.values.csv") 
DMSO2_single = csvReader_T("SD2_PPI_1_p.values.csv") 
H2O2_single = csvReader_T("H2O2_PPI_1_p.values.csv")
HU_single = csvReader_T("Hydroxyurea_PPI_1_p.values.csv")
Dox_single = csvReader_T("Doxorubicin_PPI_1_p.values.csv")
Forskolin_single = csvReader_T("Forskolin_PPI_1_p.values.csv")
Raffinose_single = csvReader_T("Raffinose_PPI_1_p.values.csv")
NaCl_single = csvReader_T("NaCl_PPI_1_p.values.csv")
cold_single = csvReader_T("Cold_16C_PPI_1_p.values.csv")
FK506_single = csvReader_T("FK506_PPI_1_p.values.csv")

PPI_list = list(DMSO_single[,1], DMSO2_single[,1], H2O2_single[,1], HU_single[,1], Dox_single[,1], 
                Forskolin_single[,1], Raffinose_single[,1], NaCl_single[,1], cold_single[,1], FK506_single[,1]) # store PPIs of one environment as an element in a list
all_PPI = unique(c(DMSO_single[,1], DMSO2_single[,1], H2O2_single[,1], HU_single[,1], Dox_single[,1], 
                   Forskolin_single[,1], Raffinose_single[,1], NaCl_single[,1], cold_single[,1], FK506_single[,1])) #1027646
all_PPI_matrix = matrix(0, length(all_PPI), 11)
all_PPI_matrix[,1] = all_PPI
colnames(all_PPI_matrix)= c("PPI","SD", "SD2","H2O2", "HU", "Dox", "Forskolin", "Raffinose", "NaCl", "16C", "FK506")
# Check each PPI in each environment, if reported, put a value of 1 into the space.
for(i in 1:nrow(all_PPI_matrix)){
        for(j in 1:length(PPI_list)){
                if(all_PPI[i] %in% PPI_list[[j]]){
                        all_PPI_matrix[i,j +1] = 1
                }
        }
        
}

environment_number = rep(0, nrow(all_PPI_matrix)) # adding the column of 2:9 in each row will give the number of positive environment
for(i in 1:length(environment_number)){
        environment_number[i] = sum(as.numeric(all_PPI_matrix[i,2:11]))
}
all_PPI_matrix_final = cbind(all_PPI_matrix[,1], environment_number, all_PPI_matrix[,2:11])
csvWriter(all_PPI_matrix_final, "PPI_single_environment_count_summary.csv")
