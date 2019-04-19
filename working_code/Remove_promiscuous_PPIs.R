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
setwd("~/Dropbox/PPiSeq_02/")
all_PPI = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv") # 13829
PPI_degree = protein_degree_count(all_PPI[,1]) # 2111
colnames(PPI_degree) = c("protein", "degree")
csvWriter(PPI_degree, "Working_data/Positive_PPI_environment/Protein_all_environment_degree.csv"  )
hist(PPI_degree[,2], breaks = seq(1,1200,by =1), xlim = c(0,100))
#### Check the number of environments in which a promiscuous protein is detected
setwd("~/Dropbox/PPiSeq_02/Working_data/Promiscuous_PPIs/")
DMSO = csvReader_T("DMSO_Promiscuous_proteins.csv")
H2O2 = csvReader_T("H2O2_Promiscuous_proteins.csv")
HU = csvReader_T("HU_Promiscuous_proteins.csv")
Dox = csvReader_T("Dox_Promiscuous_proteins.csv")
Forskolin = csvReader_T("Forskolin_Promiscuous_proteins.csv")
Raffinose = csvReader_T("Raffinose_Promiscuous_proteins.csv")
NaCl = csvReader_T("NaCl_Promiscuous_proteins.csv")
FK506 = csvReader_T("FK506_Promiscuous_proteins.csv")
cold = csvReader_T("16C_Promiscuous_proteins.csv")
extract_promiscuous_protein = function(DMSO){
        temp = unique(c(split_string_vector(DMSO[,1])[,1:2]))
        temp_HO = temp[grep("HO", temp)]
        temp_unique = temp[which(!temp %in% temp_HO)]
        return(temp_unique)
}
DMSO_bad = extract_promiscuous_protein(DMSO)
H2O2_bad = extract_promiscuous_protein(H2O2)
HU_bad = extract_promiscuous_protein(HU)
Dox_bad = extract_promiscuous_protein(Dox)
Forskolin_bad = extract_promiscuous_protein(Forskolin)
Raffinose_bad = extract_promiscuous_protein(Raffinose)
NaCl_bad = extract_promiscuous_protein(NaCl)
FK506_bad = extract_promiscuous_protein(FK506)
cold_bad = extract_promiscuous_protein(cold)

final_bad = unique(c(DMSO_bad, H2O2_bad, HU_bad, Dox_bad, Forskolin_bad,
                   Raffinose_bad, NaCl_bad, FK506_bad, cold_bad)) # 70

PPI_list = list(DMSO_bad, H2O2_bad, HU_bad, Dox_bad, Forskolin_bad,
                Raffinose_bad, NaCl_bad, FK506_bad, cold_bad)
matrix_bad = matrix(0, length(final_bad), 11)
colnames(matrix_bad) = c("PPI","Positive_sum","DMSO", "H2O2", "HU", "Dox", "Forskolin",
                         "Raffinose", "NaCl", "FK506", "Cold")
matrix_bad[,1] = final_bad
for(i in 1:nrow(matrix_bad)){
        for(j in 1:length(PPI_list)){
                if(matrix_bad[i,1] %in% PPI_list[[j]]){
                        matrix_bad[i,j +2] = 1
                }
        }
        matrix_bad[i,2] = sum(as.numeric(matrix_bad[i,3: ncol(matrix_bad)]))
}

matrix_bad_order = matrix_bad[order(matrix_bad[,2], decreasing = T),]
csvWriter(matrix_bad_order, "Promiscuous_protein_summary.csv")

# Remove PPIs that contain a protein that is identified as promiscuous protein in >= 2 environments
setwd("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Pos_PPI_not_removing_promiscuous_protein_different_environments/")
DMSO_pos = csvReader_T("DMSO_Pos_PPI_real.csv") # 5178
H2O2_pos = csvReader_T("H2O2_Pos_PPI_real.csv") # 4999
HU_pos = csvReader_T("Hydroxyurea_Pos_PPI_real.csv") # 6299
Dox_pos = csvReader_T("Doxorubicin_Pos_PPI_real.csv") # 3492
Forskolin_pos = csvReader_T("Forskolin_Pos_PPI_real.csv") # 4257
NaCl_pos = csvReader_T("NaCl_Pos_PPI_real.csv") # 2725
Raffinose_pos = csvReader_T("Raffinose_Pos_PPI_real.csv") # 4378
FK506_pos = csvReader_T("FK506_Pos_PPI_real.csv") #3472
Cold_pos = csvReader_T("Cold_16C_Pos_PPI_real.csv") # 3154
promiscuous_protein = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/Promiscuous_PPIs/Promiscuous_protein_summary.csv")
promiscuous_protein_chosen = promiscuous_protein[which(as.numeric(promiscuous_protein[,2])>=2),1]
remove_PPI_protein = function(DMSO_pos, promiscuous_protein_chosen){
       check = rep(0, nrow(DMSO_pos))
       PPI = split_string_vector(DMSO_pos[,1])
       for(i in 1:nrow(DMSO_pos)){
               if (PPI[i,1] %in% promiscuous_protein_chosen | PPI[i,2] %in% promiscuous_protein_chosen){
                       check[i] = 1
               }
       }
       DMSO_pos_filter = DMSO_pos[which(check != 1),]
       return(DMSO_pos_filter)
}
DMSO_pos_filter = remove_PPI_protein(DMSO_pos, promiscuous_protein_chosen) # 5036
H2O2_pos_filter = remove_PPI_protein(H2O2_pos, promiscuous_protein_chosen) # 4992
HU_pos_filter = remove_PPI_protein(HU_pos, promiscuous_protein_chosen) # 5150
Dox_pos_filter = remove_PPI_protein(Dox_pos, promiscuous_protein_chosen) # 3382
Forskolin_pos_filter = remove_PPI_protein(Forskolin_pos, promiscuous_protein_chosen) # 4232
NaCl_pos_filter = remove_PPI_protein(NaCl_pos, promiscuous_protein_chosen) # 2601
Raffinose_pos_filter = remove_PPI_protein(Raffinose_pos, promiscuous_protein_chosen) # 4371
FK506_pos_filter = remove_PPI_protein(FK506_pos, promiscuous_protein_chosen) # 3458
cold_pos_filter = remove_PPI_protein(Cold_pos, promiscuous_protein_chosen) # 3112

csvWriter(DMSO_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/DMSO_Pos_PPI_real.csv")
csvWriter(H2O2_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/H2O2_Pos_PPI_real.csv")
csvWriter(HU_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Hydroxyurea_Pos_PPI_real.csv")
csvWriter(Dox_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Doxorubicin_Pos_PPI_real.csv")
csvWriter(Forskolin_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Forskolin_Pos_PPI_real.csv")
csvWriter(NaCl_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/NaCl_Pos_PPI_real.csv")
csvWriter(Raffinose_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Raffinose_Pos_PPI_real.csv")
csvWriter(FK506_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/FK506_Pos_PPI_real.csv")
csvWriter(cold_pos_filter, "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Cold_16C_Pos_PPI_real.csv")
