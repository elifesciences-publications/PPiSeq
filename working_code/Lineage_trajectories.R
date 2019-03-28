## (3) Make lineage trajectories of representative lineages 
setwd("/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files_final/DBC_counts/max_5/")
library("dplyr")
library("tidyr")
library(ggplot2)
library(plotly)
library(gridExtra)
library(reshape)
library(ggpubr)
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

PPI_lineages = csvReader_T("known_PPI_fitness_barcodes_filtered.csv")

PPI_lineages_select = PPI_lineages[, c(1, 3, 5, 7:11)]
for (i in 4:8){
        PPI_lineages_select[,i] = frequency(as.numeric(PPI_lineages_select[,i]))
}

PPI_lineages_select[PPI_lineages_select == "0"] = 1e-9 # 5081689

csvWriter(PPI_lineages_select, "Lineage_plot/PPI_barcode_trajectories_fitness_plot.csv")

filter_trajectories = function(a){
        matrix_final_select= matrix(0, nrow(a), ncol(a))
        colnames(matrix_final_select) = colnames(a)
        for (i in 1: nrow(a)){
                l= which(as.numeric(a[i,4:8]) == 1e-9)
                if (length(l) == 0){
                        matrix_final_select[i,]= a[i,]
                }
                else if (identical(l, (6-length(l)) : 5)){ # 5 = 8-4+1; 6 = 5 +1
                        matrix_final_select[i,]= a[i,]
                }
                else{
                        matrix_final_select[i,]= NA
                }
        }
        matrix_final_select= na.omit(matrix_final_select) 
        return(matrix_final_select)
}

PPI_lineages_final = filter_trajectories(PPI_lineages_select) # 4993630
csvWriter(PPI_lineages_final, "Lineage_plot/Good_PPI_barcode_trajectories_fitness_plot.csv")

setwd("/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files_final/DBC_counts/max_5/Lineage_plot/")
PPI_lineages = dataFrameReader_T("PPI_barcode_trajectories_fitness_plot.csv") ### Not choose good trajectories

PPI_pos_DHFR = PPI_lineages[which(PPI_lineages[,1] == "positive_DHFR"),]
PPI_pos_DHFR[,1] = "Positive control"
PPI_neg_DHFR = PPI_lineages[which(PPI_lineages[,1] == "negative_non_DHFR"),]
PPI_neg_DHFR[,1] = "Negative control"

# Get PPIs containing IMD3 (YLR432W) as bait
PPI_ERP3 = PPI_lineages[grep("YDL018C_", PPI_lineages[,1]),] # 5258

PPI_IMD3 = PPI_lineages[grep("_YLR432W", PPI_lineages[,1]),] # 3545
PPI_promiscuous_YPL139C = PPI_lineages[grep("YPL139C_", PPI_lineages[,1]),] # 6410
PPI_promiscuous_YIL143C = PPI_lineages[grep("_YIL143C", PPI_lineages[,1]), ] # 3425

csvWriter(PPI_ERP3, "ERP3-DHFR[3] X ORF-DHFR[1,2].csv")
csvWriter(PPI_IMD3, "ORF-DHFR[3] X IMD3-DHFR[1,2].csv")
csvWriter(PPI_promiscuous_YIL143C, "ORF-DHFR[3] X YIL143C-DHFR[1,2].csv")
csvWriter(PPI_promiscuous_YPL139C, "YPL139C-DHFR[3] X ORF-DHFR[1,2].csv ")
#PPI_control = do.call("rbind", list(PPI_pos_DHFR, PPI_neg_DHFR))

PRS = PPI_lineages[grep("Pos", PPI_lineages[,1]),]
RRS = PPI_lineages[grep("Neg", PPI_lineages[,1]),]

PRS[,1] = "PRS"
RRS[,1] = "RRS"

#PPI_reference = do.call("rbind", list(PRS, RRS))

#### Fragment strains
MATa_DHFR12 = csvReader_T("/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files_final/PPI_barcodes/MATa_genome_combine.csv")
MATalpha_DHFR3 = csvReader_T("/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files_final/PPI_barcodes/MATalpha_genome_combine.csv")

PPI_DHFR12 = PPI_lineages[which(PPI_lineages[,2] %in% MATa_DHFR12[,3]),]
PPI_DHFR12[,1] = "ORF X DHFR[1,2]"

PPI_DHFR3 = PPI_lineages[which(PPI_lineages[,2] %in% MATalpha_DHFR3[,3]),]
PPI_DHFR3[,1] = "DHFR[3] X ORF"

PPI_pos_DHFR_transform = transform_lineages(PPI_pos_DHFR, 1,2,3, 4)
PPI_neg_DHFR_transform = transform_lineages(PPI_neg_DHFR, 1,2,3,4)
PRS_transform = transform_lineages(PRS, 1,2,3,4)
RRS_transform = transform_lineages(RRS, 1,2,3,4)
PPI_DHFR12_transform = transform_lineages(PPI_DHFR12, 1,2,3,4)
PPI_DHFR3_transform = transform_lineages(PPI_DHFR3, 1,2,3,4)

PPI_ERP3_transform = transform_lineages(PPI_ERP3, 1,2,3,4)
PPI_IMD3_transform = transform_lineages(PPI_IMD3, 1,2,3,4)

PPI_promiscuous_YPL139C_transform = transform_lineages(PPI_promiscuous_YPL139C, 1,2,3,4)
PPI_promiscuous_YIL143C_transform = transform_lineages(PPI_promiscuous_YIL143C, 1,2,3,4)


library(scales)

Lineage_plot = function(PPI_pos_DHFR_transform, output){
        ggplot(data = PPI_pos_DHFR_transform, aes(x = time_points, y= lineages, group = barcode, color = Fitness)) +
                geom_line(alpha =0.7, size = 0.3) +
                scale_color_gradientn(colors = apple_colors[c(5,10,7)], limits = c(-0.5, 1.3)) +
                scale_y_continuous(name = "Frequency",
                                   limits = c(1e-9, 1e-1),
                                   trans = log10_trans(),
                                   breaks = trans_breaks("log10", function(x) 10^x),
                                   labels = trans_format("log10", math_format(10^.x)))+
                scale_x_continuous(name = "Generation",
                                   breaks = seq(0,12, by =3),
                                   labels = seq(0,12, by= 3)) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                theme(axis.text.x = element_text(size = 10, color = "black"), 
                      axis.text.y.left = element_text(size = 10, color = "black")) + 
                theme(text = element_text(size=12))
        
        ggsave(output, width = 5, height = 4)
        
}

Lineage_plot(PPI_pos_DHFR_transform, "Positive DHFR lineages.pdf")
Lineage_plot(PPI_neg_DHFR_transform, "Negative DHFR lineages.pdf")
Lineage_plot(PRS_transform, "PRS lineages.pdf")
Lineage_plot(RRS_transform, "RRS lineages.pdf")
Lineage_plot(PPI_DHFR12_transform, "ORF X DHFR[1,2] lineages.pdf")
Lineage_plot(PPI_DHFR3_transform, "DHFR[3] X ORF lineages.pdf")

Lineage_plot(PPI_ERP3_transform, "ERP3 X ORF lineages.pdf")
Lineage_plot(PPI_IMD3_transform, "ORF X IMD3 lineages.pdf")
Lineage_plot(PPI_promiscuous_YIL143C_transform, "ORF X YIL143C.pdf")
Lineage_plot(PPI_promiscuous_YPL139C_transform, "YPL139C X ORF.pdf")

#### Choose some lineages and make a plot
all_control = do.call("rbind", list(PPI_pos_DHFR, PPI_neg_DHFR, PRS, RRS, PPI_DHFR12, PPI_DHFR3))

PPI_lineages_remaining = PPI_lineages[which(!PPI_lineages[,2] %in% all_control[,2]),]
PPI_lineages_remaining = PPI_lineages_remaining[order(PPI_lineages_remaining[,3], decreasing = T),]

PPI_lineages_high = PPI_lineages_remaining[which(PPI_lineages_remaining[,3] > 0.4 & PPI_lineages_remaining[,4] > 1e-7),]
PPI_lineages_medium = PPI_lineages_remaining[which(PPI_lineages_remaining[,3] > 0.1 & PPI_lineages_remaining [,3] < 0.4
                                                   &PPI_lineages_remaining[,4] > 1e-8 ),]
PPI_lineages_neutral = PPI_lineages_remaining[which(PPI_lineages_remaining[,3] > -0.2 & PPI_lineages_remaining [,3] < 0.2
                                                    &PPI_lineages_remaining[,4] > 1e-8),]
PPI_lineages_low = PPI_lineages_remaining[which(PPI_lineages_remaining[,3] > -0.5 & PPI_lineages_remaining [,3] < -0.1
                                                &PPI_lineages_remaining[,4] > 1e-8),]

PPI_lineages_final = do.call("rbind", list(PPI_lineages_high[sample(1:nrow(PPI_lineages_high), 1000),],
                                           PPI_lineages_medium[sample(1:nrow(PPI_lineages_medium), 2000),],
                                           PPI_lineages_neutral[sample(1:nrow(PPI_lineages_neutral), 2000),],
                                           PPI_lineages_low[sample(1:nrow(PPI_lineages_low), 1000),]))
PPI_lineages_final_transform = transform_lineages(PPI_lineages_final, 1,2,3,4)

Lineage_plot(PPI_lineages_final_transform, "6000 representative lineages.pdf")