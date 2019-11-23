
source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
####Check the initial DMSO environment sequencing data (protein pair number, barcode number, theoretical number)
setwd("~/Dropbox/PPiseq_02/")
DMSO_lineage = csvReader_T("/Volumes/zmliu_02/PPiseq/DMSO/counts/lineage_trajectories.csv")
nrow(DMSO_lineage) # 6,855,881
barcodes_control = csvReader_T("paper_data/PPiseq_barcodes/PPiseq_control_barcodes.csv")
barcodes_real = csvReader_T("paper_data/PPiseq_barcodes/PPiseq_expected_barcodes_non_control.csv")
barcodes_expected = rbind(barcodes_control, barcodes_real)
nrow(barcodes_expected) # 7,741,081
lineage_T0 = DMSO_lineage[which(as.numeric(DMSO_lineage[,2]) > 0),] # 6185226

counts = as.numeric(lineage_T0[,2])

library(ggplot2)
library(scales)
ggplot() + aes(counts) +
        geom_histogram(breaks = seq(0,46000,by=10), color = apple_colors[5], fill = apple_colors[5]) +
        scale_y_continuous(name = "Frequency",
                           trans = log10_trans(),
                           breaks = trans_breaks("log10", function(x) 10^x),
                           labels = trans_format("log10", math_format(10^.x)))+
        scale_x_continuous(name = "Barcode counts (bin size =10)",
                           limits=c(0, 10000),
                           breaks = seq(0,10000, by =2500),
                           labels = seq(0,10000, by= 2500)) +
        
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 12, color = "black"), 
              axis.text.y.left = element_text(size = 12, color = "black")) + 
        theme(text = element_text(size=14))

ggsave("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/FigureS1A_barcode_counts_distribution_01.pdf", width = 5, height = 5)

#### plot initial barcode counts distribution non-log
ggplot() + aes(counts) +
        geom_histogram(breaks = seq(0,46000,by=1), color = apple_colors[5], fill = apple_colors[5]) +
        #labs(title = "Initial barcode counts of PPiSeq library") + 
        scale_y_continuous(name = "Frequency",
                           limits = c(0,6e5),
                           breaks = seq(0,6e5, by =1e5),
                           labels = fancy_scientific(seq(0,6e5, by =1e5))) +
        scale_x_continuous(name = "Barcode counts (bin size = 1)",
                           limits=c(0, 100),
                           breaks = seq(0,100, by =10),
                           labels = seq(0,100, by= 10)) +
        
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x = element_text(size = 12, color = "black"), 
              axis.text.y.left = element_text(size = 12, color = "black")) + 
        theme(text = element_text(size=14))
ggsave("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/FigureS1A_barcode_counts_distribution_02.pdf", width = 6, height = 4)

### pie plot to show the number of barcoded replicates for each PPI
known_covered = length(which(barcodes_expected[,2] %in% lineage_T0[,1])) # 6048285
ratio_known = known_covered/nrow(barcodes_expected)
ratio_known # 0.781323
barcodes_known_covered = barcodes_expected[which(barcodes_expected[,2] %in% lineage_T0[,1]), ]
PPI_known_covered = unique(barcodes_known_covered[,1])
a = length(PPI_known_covered) # 1793976
b = length(unique(barcodes_expected[,1])) # 1942950
ratio_PPI_known = a/b
ratio_PPI_known # 0.9233259


## Pie plot to show the barcode count distribution for each PPI
barcodes_real_covered = barcodes_real[which(barcodes_real[,2] %in% lineage_T0[,1]), ]
PPI_covered = barcodes_real_covered[,1]
PPI_count = as.data.frame(table(PPI_covered))
PPI_count_hist = as.data.frame(table(PPI_count$Freq))
PPI_count_hist
colnames(PPI_count_hist) = c("barcode_count", "Freq" )
#Var1    Freq
# 1  121857
# 2  240960
# 3  274270
# 4 1150205
# On average 
average = (121857* 1 + 240960*2 + 274270*3 + 1150205*4)/sum(PPI_count_hist$Freq)
average

library(ggplot2)
ggplot(PPI_count_hist, aes(x = factor(1), y = Freq, fill= barcode_count)) +
        geom_bar(width =1, stat = "identity") + 
        geom_text(aes(label = Freq), position = position_stack(vjust = 0.5), cex = 2) +
        coord_polar("y") +
        theme_minimal() +
        theme(axis.text = element_blank(),
              panel.border = element_blank(),
              panel.grid=element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank()) +
        scale_fill_manual(values=apple_colors[c(3,4,5,6)],
                          name = "Barcodes/PPI",
                          breaks = c(1,2,3,4),
                          labels = c("1", "2", "3", "4"))
#guides(fill=FALSE)
ggsave("~/Dropbox/PPiseq_02/Working_figure/SFigures/paper/FigureS1B_barcode_counts_per_PPI_real.pdf", width = 5, height = 5)


