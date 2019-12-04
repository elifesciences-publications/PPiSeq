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

##### DNA templated transcription elongation correlated PPIs
# (1) DNA tempalted transcription, elongation: GO:0006354
# (2) Nucleobase containing small molecule, metabolic process: GO:0055086
# (3) Nucleobase-containing compound transport： GO:0015931
# (4) Peptidyl-amino acid modification: GO:0018193
setwd("~/Dropbox/PPiSeq_02/")
PPI_fit = csvReader_T("Working_data/Positive_PPI_environment/Variation_score_PPI_environment_primary.csv")
GO_slim = as.matrix(read.table("Working_data/GO_term_files/go_slim_mapping_tab_20190405.txt", header = F, sep = "\t"))
Gene_trans_initiation = unique(GO_slim[which(GO_slim[,6] == "GO:0006352"), 1])
Gene_trans_elongation = unique(GO_slim[which(GO_slim[,6] == "GO:0006354"), 1])
Gene_nucleobase_transport = unique(GO_slim[which(GO_slim[,6] == "GO:0015931"), 1])
Gene_nucleobase_metabolic = unique(GO_slim[which(GO_slim[,6] == "GO:0055086"), 1])
Gene_peptidyl_modification = unique(GO_slim[which(GO_slim[,6] == "GO:0018193"), 1])

#Gene_translation = unique(GO_slim[which(GO_slim[,6] == "GO:0006413")])
check_specific_protein = function(PPI, Gene_Carbon){
        PPI_chosen = "0"
        protein_pair = split_string_vector(PPI[,1])
        if(length(Gene_Carbon) > 1){
                for(i in 1:nrow(protein_pair)){
                        if(protein_pair[i,1] %in% Gene_Carbon | protein_pair[i,2] %in% Gene_Carbon){
                                PPI_chosen = c(PPI_chosen, PPI[i,1])
                        }
                }  
        }else {
                for(i in 1:nrow(protein_pair)){
                        if(protein_pair[i,1] == Gene_Carbon | protein_pair[i,2] == Gene_Carbon){
                                PPI_chosen = c(PPI_chosen, PPI[i,1])
                        }
                }  
        }
        
        PPI_chosen = PPI_chosen[2:length(PPI_chosen)]
        return(PPI_chosen)
}

#PPI_carbon = check_specific_protein(PPI_fit, Gene_carbon) # 284

PPI_trans_initiation = check_specific_protein(PPI_fit, Gene_trans_initiation) # 634
PPI_trans_elongation = check_specific_protein(PPI_fit, Gene_trans_elongation ) # 525
PPI_nucleobase_metabolic = check_specific_protein(PPI_fit, Gene_nucleobase_metabolic) #1210
PPI_nucleobase_transport = check_specific_protein(PPI_fit, Gene_nucleobase_transport ) # 832
PPI_peptidyl_modification = check_specific_protein(PPI_fit, Gene_peptidyl_modification) # 557

Protein_trans_initiation = as.data.frame(table(split_string_vector(PPI_trans_initiation )[,1:2]))
Protein_trans_initiation = Protein_trans_initiation[order(Protein_trans_initiation$Freq, decreasing = T),]
# Var1 Freq
# YMR112C  360
# YGL070C  109
# YNL118C   26
# YPR052C   24
# YKL058W   13

Protein_trans_elongation = as.data.frame(table(split_string_vector(PPI_trans_elongation)[,1:2]))
Protein_trans_elongation = Protein_trans_elongation[order(Protein_trans_elongation$Freq, decreasing = T),]
#Var1 Freq
#YGL070C  109
#YKL139W   96
#YKR048C   44
#YNL288W   31
Protein_nucleobase_met = as.data.frame(table(split_string_vector(PPI_nucleobase_metabolic)[,1:2]))
Protein_nucleobase_met = Protein_nucleobase_met[order(Protein_nucleobase_met$Freq, decreasing = T),]
# Var1 Freq
# YKR080W  196
# YLR432W  120
# YER070W  119
Protein_nucleobase_transport = as.data.frame(table(split_string_vector(PPI_nucleobase_transport)[,1:2]))
Protein_nucleobase_transport = Protein_nucleobase_transport[order(Protein_nucleobase_transport$Freq, decreasing = T),]
# Var1 Freq
# YKL002W  215
# YNL007C   80
# YNL064C   69
# YER063W   57
# YER056C   49
Protein_peptidyl = as.data.frame(table(split_string_vector(PPI_peptidyl_modification)[,1:2]))
Protein_peptidyl = Protein_peptidyl[order(Protein_peptidyl$Freq, decreasing = T),]
# Var1 Freq
# YJL002C  166
# YKL139W   96
# YHR068W   29
# YJR070C   20
# YBL071W-A   15

####### Tecan data validation of transcription initiation correlated PPIs
TECAN_data_AUC_validate = function(data_file_name, map_name,time_window,time_window_low, matrix_name, bad_wells){
        ## robot_data_01 as MTX data
        robot_data_01= read.csv(data_file_name, skip = 19)
        robot_data_01 = robot_data_01[1:(nrow(robot_data_01)-2),]
        #Matix of ODs
        od_robot_01 = as.matrix(robot_data_01[, 3:ncol(robot_data_01)])
        map_index = 1:96
        od_robot_01 = od_robot_01[,which(!map_index %in% bad_wells)] # remove the data of bad wells
        #Vector of Times (s)
        y = as.vector(robot_data_01[,1])
        z = strsplit(y, "\t")
        time_robot_01 = 1:length(z)
        for(i in 1:length(time_robot_01)){
                time_robot_01[i] = as.numeric(z[[i]][2])/60
        }
        ## input the map 
        diploid_05_02 = read.csv(map_name)
        PPI_map = diploid_05_02[,1]
        PPI_map = PPI_map[which(!map_index %in% bad_wells)]
        PPI_unique = unique(PPI_map)
        PPI_unique = PPI_unique[which(PPI_unique != "Negative")]
        
        library(pracma)
        AUC_PPI = function(time, od_data, column_PPI){
                auc_PPI = rep(0, length(column_PPI))
                for(i in 1: length(column_PPI)){
                        auc_PPI[i] = trapz(time, od_data[,column_PPI[i]])
                }
                return(auc_PPI)
        }
        
        ## Calculate AUC 
        column_negative_PPI = which(PPI_map == "Negative")
        time_select_01 = time_robot_01[which(time_robot_01/60 <= time_window & time_robot_01/60 >= time_window_low)] ## only consider data before MTX_time_window_threshold
        od_robot_select_01 = od_robot_01[1:length(time_select_01),] # only consider data before MTX_time_window_threshold
        Negative_AUC_01 = AUC_PPI(time_select_01, od_robot_select_01, column_negative_PPI)
        
        
        ### Create a matrix that contain normalized AUCs (subtract negative) 
        ### for each PPI, and calculate mean and SD of these SDs
        auc_PPI_matrix = matrix(0, length(PPI_unique), 6)
        colnames(auc_PPI_matrix) = c("PPI","Mean_Neg_AUC", "SD_Neg_AUC", "Mean_target_AUC", 
                                     "SD_target_AUC", "adjust_pvalue")
        auc_PPI_matrix[,1] = as.character(PPI_unique)
        auc_PPI_matrix[,2] = mean(Negative_AUC_01)
        auc_PPI_matrix[,3] = sd(Negative_AUC_01)
        Neg_AUC = mean(Negative_AUC_01)
        for(i in 1:length(PPI_unique)){
                column_each_PPI = which(PPI_map == PPI_unique[i])
                auc_each_PPI_01 = AUC_PPI(time_select_01, od_robot_select_01, column_each_PPI)
                auc_PPI_matrix[i,4] = mean(auc_each_PPI_01)
                auc_PPI_matrix[i,5] = sd(auc_each_PPI_01)
                auc_PPI_matrix[i,6] = t.test(auc_each_PPI_01, Negative_AUC_01, alternative = "greater")$p.value
        }
        auc_PPI_matrix[,6] = p.adjust(auc_PPI_matrix[,6], "BH")
        write.csv(auc_PPI_matrix, matrix_name, row.names = F)
        #return(length(which(as.numeric(auc_PPI_matrix[,6]) <= 0.05 & as.numeric(auc_PPI_matrix[,5]) >= 0)))
        return(length(which(as.numeric(auc_PPI_matrix[,6]) <= 0.05)))
}

setwd("~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/transcription/data/")
data_file_name = "2019-08-20_16-21_Diploid_transcription_01_SD_MTX_96_T7.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/transcription/Final_map/Diploid_transcription_1.csv" 
time_window = 42
time_window_low = 0
bad_wells = c(94,95) ## Two very different wells were removed
matrix_name = "Validation_diploid_transcription_01_SD_MTX_T7.csv"
TECAN_data_AUC_validate(data_file_name, map_name, time_window, time_window_low, matrix_name, bad_wells) #18 from 30

data_file_name = "2019-08-22_14-17_Diploid_transcription_02_SD_MTX_96_T14.txt"
map_name = "~/Dropbox/PPiSeq_02/Working_data/TECAN_validation/transcription/Final_map/Diploid_transcription_2.csv" 
time_window = 42
time_window_low = 0
bad_wells = 0
matrix_name = "Validation_diploid_transcription_02_SD_MTX_T14.csv"
TECAN_data_AUC_validate(data_file_name, map_name, time_window, time_window_low, matrix_name, bad_wells) #24 from 30

transcription_01 = csvReader_T("Validation_diploid_transcription_01_SD_MTX_T7.csv")
transcription_02 = csvReader_T("Validation_diploid_transcription_02_SD_MTX_T14.csv")
transcription_validate = rbind(transcription_01, transcription_02)
length(which(as.numeric(transcription_validate[,6]) <= 0.05))
csvWriter(transcription_validate, "Validation_diploid_transcription_1_2_SD_MTX.csv")

pdf("~/Dropbox/PPiSeq_02/Working_figure/Figure1_extra/Tecan_validate_transcription.pdf", width = 3, height = 5)
barcenter = barplot(c(60, 42), col = apple_colors[c(1,2)], ylab = c("Number of positive PPIs"))
text(barcenter, -4, labels = c("PPiSeq", "OD595"), xpd = T)
text(barcenter, c(62, 45), labels = c("60", "42"), xpd = T)
dev.off()