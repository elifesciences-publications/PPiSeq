source("/Volumes/zmliu_02/PPiseq/DMSO/all_lintag_files_final/R_code/function.R")
setwd("/Volumes/zmliu_02/PPiseq_02/DMSO/reference_set/")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
SD_metrics = dataFrameReader_T("Threshold_metrics_final.csv")
overlap_PCA = SD_metrics$Overlap_Michnick/SD_metrics$all
overlap_BioGRID = SD_metrics$Overlap_BioGRID/SD_metrics$all
pdf("~/Dropbox/PPiSeq_02/Working_figure/SFigures/Figure2_related/Percentage_reprted_PPI_figure2A.pdf", 
    width =5, height =5)
plot(SD_metrics$Combined_PPV, overlap_PCA, type= "p", pch = 16, col = apple_colors[1],
     xlim = c(0.5, 0.9), ylim = c(0, 0.5), xlab = "Dynamic threshold", 
     ylab = "Percentage of overlapped PPIs", bty = "n")
text(SD_metrics$Combined_PPV, rep(0.03, nrow(SD_metrics)), as.character(SD_metrics$all), 
     srt =90,cex = 0.6)
points(SD_metrics$Combined_PPV, overlap_BioGRID, type = "p", pch = 16, col = apple_colors[2])
legend("topleft", c("PCA", "BioGRID containing PCA"), pch = c(16, 16), col = apple_colors[1:2],
       bty = "n")
dev.off()