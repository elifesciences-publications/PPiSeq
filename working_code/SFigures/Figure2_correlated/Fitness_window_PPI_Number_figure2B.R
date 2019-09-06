source("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/R_code/function.R")
setwd("/Volumes/zmliu_02/PPiseq/DMSO/counts/")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

DMSO = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/SD_mean_fitness_positive.csv")
DMSO_DHFR_pos = as.numeric(DMSO[grep("pos", DMSO[,1]),3]) # 0.970635
DMSO_DHFR_neg = as.numeric(DMSO[grep("negative", DMSO[,1]),3]) # 0.07146862
DMSO_pos = DMSO[which(DMSO[,7] == "1"),]
min(as.numeric(DMSO_pos[,3])) # 0.24941
DMSO_window = DMSO_DHFR_pos - DMSO_DHFR_neg # 0.8991664

Forskolin = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/Forskolin_mean_fitness_positive.csv")
Forskolin_DHFR_pos = as.numeric(Forskolin[grep("pos", Forskolin[,1]),3]) 
Forskolin_DHFR_pos #1.16344
Forskolin_DHFR_neg = as.numeric(Forskolin[grep("negative", Forskolin[,1]),3]) 
Forskolin_DHFR_neg #1.16344
Forskolin_pos = Forskolin[which(Forskolin[,7] == "1"),]
min(as.numeric(Forskolin_pos[,3])) # 0.42258
Forskolin_window = Forskolin_DHFR_pos - Forskolin_DHFR_neg # 1.079712

FK506 = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/FK506_mean_fitness_positive.csv")
FK506_DHFR_pos = as.numeric(FK506[grep("pos", FK506[,1]),3]) 
FK506_DHFR_pos #1.69887
FK506_DHFR_neg = as.numeric(FK506[grep("negative", FK506[,1]),3]) 
FK506_DHFR_neg #0.1773712
FK506_pos = FK506[which(FK506[,7] == "1"),]
min(as.numeric(FK506_pos[,3])) #0.54274
FK506_window = FK506_DHFR_pos - FK506_DHFR_neg #1.521499

NaCl = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/NaCl_mean_fitness_positive.csv")
NaCl_DHFR_pos = as.numeric(NaCl[grep("pos", NaCl[,1]),3]) 
NaCl_DHFR_pos # 0.3372279
NaCl_DHFR_neg = as.numeric(NaCl[grep("negative", NaCl[,1]),3]) 
NaCl_DHFR_neg # 0.0934315
NaCl_pos = NaCl[which(NaCl[,7] == "1"),]
min(as.numeric(NaCl_pos[,3])) # 0.0934315
NaCl_window = NaCl_DHFR_pos - NaCl_DHFR_neg # 0.3293516

Raffinose = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/Raffinose_mean_fitness_positive.csv")
Raffinose_DHFR_pos = as.numeric(Raffinose[grep("pos", Raffinose[,1]),3]) 
Raffinose_DHFR_pos # 1.83227
Raffinose_DHFR_neg = as.numeric(Raffinose[grep("negative", Raffinose[,1]),3]) 
Raffinose_DHFR_neg # 0.1205733
Raffinose_pos = Raffinose[which(Raffinose[,7] == "1"),]
min(as.numeric(Raffinose_pos[,3])) # 0.5568825
Raffinose_window = Raffinose_DHFR_pos - Raffinose_DHFR_neg # 1.711697

HU = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/Hydroxyurea_mean_fitness_positive.csv")
HU_DHFR_pos = as.numeric(HU[grep("pos", HU[,1]),3]) 
HU_DHFR_pos # 0.656245
HU_DHFR_neg = as.numeric(HU[grep("negative", HU[,1]),3]) 
HU_DHFR_neg # 0.01219391
HU_pos = HU[which(HU[,7] == "1"),]
min(as.numeric(HU_pos[,3])) # 0.1920075
HU_window = HU_DHFR_pos - HU_DHFR_neg # 0.6440511

H2O2 = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/H2O2_mean_fitness_positive.csv")
H2O2_DHFR_pos = as.numeric(H2O2[grep("pos", H2O2[,1]),3]) 
H2O2_DHFR_pos # 0.747422
H2O2_DHFR_neg = as.numeric(H2O2[grep("negative", H2O2[,1]),3]) 
H2O2_DHFR_neg # 0.01421095
H2O2_pos = H2O2[which(H2O2[,7] == "1"),]
min(as.numeric(H2O2_pos[,3])) # 0.198405
H2O2_window = H2O2_DHFR_pos - H2O2_DHFR_neg # 0.733211

Dox = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/Doxorubicin_mean_fitness_positive.csv")
Dox_DHFR_pos = as.numeric(Dox[grep("pos", Dox[,1]),3]) 
Dox_DHFR_pos # 0.595295
Dox_DHFR_neg = as.numeric(Dox[grep("negative", Dox[,1]),3]) 
Dox_DHFR_neg #  0.02563448
Dox_pos = Dox[which(Dox[,7] == "1"),]
min(as.numeric(Dox_pos[,3])) # 0.14729
Dox_window = Dox_DHFR_pos - Dox_DHFR_neg # 0.5696605

cold = csvReader_T("~/Dropbox/PPiSeq_02/Paper_data/Cold_16C_mean_fitness_positive.csv")
cold_DHFR_pos = as.numeric(cold[grep("pos", cold[,1]),3]) 
cold_DHFR_pos # 0.281273
cold_DHFR_neg = as.numeric(cold[grep("negative", cold[,1]),3]) 
cold_DHFR_neg # 0.007256943
cold_pos = cold[which(cold[,7] == "1"),]
min(as.numeric(cold_pos[,3])) # 0.0719585
cold_window = cold_DHFR_pos - cold_DHFR_neg # 0.2740161

window_all = c(DMSO_window, Forskolin_window, NaCl_window, FK506_window,
               Dox_window, Raffinose_window, H2O2_window, cold_window, HU_window)
PPI_number_all = c(nrow(DMSO_pos), nrow(Forskolin_pos), nrow(NaCl_pos),
                   nrow(FK506_pos), nrow(Dox_pos), nrow(Raffinose_pos),
                   nrow(H2O2_pos), nrow(cold_pos), nrow(HU_pos))
pdf("~/Dropbox/PPiSeq_02/Working_figure/SFigures/Figure2_related/Fitness_window_PPI_number_related_figure2B.pdf", 
    width =5, height =5)
plot(window_all, PPI_number_all, type = "p", col = apple_colors[1], bty = "n", 
     ylim = c(0, 7000), xlim = c(0, 2), pch = 16, xlab = "Fitness window",
     ylab = "Number of PPIs")
PPI_number_all[7] = PPI_number_all[7] - 600
text(window_all, (PPI_number_all + 300), cex = 0.7,
     c("SD", "Forskolin","NaCl", "FK506", "Doxorubicin", "Raffinose", 
       expression('H'[2]* 'O'[2]),  "16 \u00B0C", "Hydroxyurea"), xpd = TRUE)
dev.off()