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
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

PPI_lineages = dataFrameReader_T("Paper_data/SD_PPI_barcodes_fitness_counts.csv")
PPI_lineages_select = PPI_lineages[, c(1, 3, 4, 6:10)] 
for (i in 4:8){
        PPI_lineages_select[,i] = frequency(as.numeric(PPI_lineages_select[,i]))
}
PPI_lineages_select[PPI_lineages_select == "0"] = 1e-9 # 5013513

# Extract Positive and Negative control
PPI_pos_DHFR = PPI_lineages_select[which(PPI_lineages_select[,1] == "positive_DHFR"),]
PPI_pos_DHFR[,1] = "DHFR(+)"
PPI_neg_DHFR = PPI_lineages_select[which(PPI_lineages_select[,1] == "negative_non_DHFR"),]
PPI_neg_DHFR[,1] = "DHFR(-)"

# Get PPIs containing IMD3 (YLR432W) as bait, EPR3 as prey, 
# YIL143C (promiscuous protein) as bit, and YPL139C as prey. 
PPI_ERP3 = PPI_lineages_select[grep("YDL018C_", PPI_lineages_select[,1]),] # 5172
PPI_IMD3 = PPI_lineages_select[grep("_YLR432W", PPI_lineages_select[,1]),] # 3410
PPI_promiscuous_YPL139C = PPI_lineages_select[grep("YPL139C_", PPI_lineages_select[,1]),] # 6396
PPI_promiscuous_YIL143C = PPI_lineages_select[grep("_YIL143C", PPI_lineages_select[,1]), ] # 3367

# Extract PRS and RRS strains
PRS = PPI_lineages_select[grep("Pos", PPI_lineages_select[,1]),]
RRS = PPI_lineages_select[grep("Neg", PPI_lineages_select[,1]),]
PRS[,1] = "PRS"
RRS[,1] = "RRS"

# Extract strains that contain ORF X Fragment only strains
MATa_DHFR12 = csvReader_T("~/Dropbox/PPiSeq_02/Working_data/PPI_barcodes/MATa_genome_combine.csv")
MATalpha_DHFR3 = csvReader_T("~/Dropbox/PPiSeq_02/Working_data//PPI_barcodes/MATalpha_genome_combine.csv")
PPI_DHFR12 = PPI_lineages_select[which(PPI_lineages_select[,2] %in% MATa_DHFR12[,3]),]
PPI_DHFR12[,1] = "ORF X DHFR[1,2]"
PPI_DHFR3 = PPI_lineages_select[which(PPI_lineages_select[,2] %in% MATalpha_DHFR3[,3]),]
PPI_DHFR3[,1] = "DHFR[3] X ORF"

##### A function to transform lineage matrix to data frame which can be plotted in ggplot2
transform_lineages <- function(x, PPI_pos, barcode_pos, fitness_pos, G0_pos){
        time_points = rep(c(0,3,6,9,12), nrow(x))
        lineages = as.numeric(c(t(x[, G0_pos:(G0_pos + 4)])))
        PPI = rep("0", length(time_points))
        barcode = rep("0", length(time_points))
        Fitness = rep(0, length(time_points))
        for(i in 1:nrow(x)){
                PPI[(i*5 -4): (i*5)] = as.character(rep(x[i,PPI_pos], 5))
                barcode[(i*5 -4): (i*5)] = as.character(rep(x[i,barcode_pos], 5))
                Fitness[(i*5 -4): (i*5)] = as.numeric(rep(x[i,fitness_pos], 5))
        }
        matrix_lineages = data.frame(PPI, barcode, Fitness, time_points, lineages)
        return (matrix_lineages)
}
# Trasnform the lineages into a format for plot
PPI_pos_DHFR_transform = transform_lineages(PPI_pos_DHFR, 1,2,3,4)
PPI_neg_DHFR_transform = transform_lineages(PPI_neg_DHFR, 1,2,3,4)
PRS_transform = transform_lineages(PRS, 1,2,3,4)
RRS_transform = transform_lineages(RRS, 1,2,3,4)
PPI_DHFR12_transform = transform_lineages(PPI_DHFR12, 1,2,3,4)
PPI_DHFR3_transform = transform_lineages(PPI_DHFR3, 1,2,3,4)
PPI_ERP3_transform = transform_lineages(PPI_ERP3, 1,2,3,4)
PPI_IMD3_transform = transform_lineages(PPI_IMD3, 1,2,3,4)
PPI_promiscuous_YPL139C_transform = transform_lineages(PPI_promiscuous_YPL139C, 1,2,3,4)
PPI_promiscuous_YIL143C_transform = transform_lineages(PPI_promiscuous_YIL143C, 1,2,3,4)

library(ggplot2)
library(scales)

Lineage_plot = function(PPI_pos_DHFR_transform, output){
        PPI_pos_DHFR_transform = PPI_pos_DHFR_transform
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
                theme(text = element_text(size=10, color = "black"))
        
        ggsave(output, width = 5, height = 4)
        
}

Lineage_plot(PPI_pos_DHFR_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/Positive DHFR lineages.pdf")
Lineage_plot(PPI_neg_DHFR_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/Negative DHFR lineages.pdf")
Lineage_plot(PRS_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/PRS lineages.pdf")
Lineage_plot(RRS_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/RRS lineages.pdf")
Lineage_plot(PPI_DHFR12_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/ORF X DHFR[1,2] lineages.pdf")
Lineage_plot(PPI_DHFR3_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/DHFR[3] X ORF lineages.pdf")
Lineage_plot(PPI_ERP3_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/ERP3 X ORF lineages.pdf")
Lineage_plot(PPI_IMD3_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/ORF X IMD3 lineages.pdf")
Lineage_plot(PPI_promiscuous_YIL143C_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/ORF X YIL143C.pdf")
Lineage_plot(PPI_promiscuous_YPL139C_transform, "Working_figure/SFigures/FigureSX_SD_lineage_plot_related_Fig1B/YPL139C X ORF.pdf")

#Randomly choose some lineages after removing these controls strains and make a plot
all_control = do.call("rbind", list(PPI_pos_DHFR, PPI_neg_DHFR, PRS, RRS, PPI_DHFR12, PPI_DHFR3))# Combine all the controls
PPI_lineages_remaining = PPI_lineages_select[which(!PPI_lineages_select[,2] %in% all_control[,2]),]
PPI_lineages_remaining = PPI_lineages_remaining[order(as.numeric(PPI_lineages_remaining[,3]), decreasing = T),] # order the lineage based on their fitness

# Different groups
PPI_lineages_high = PPI_lineages_remaining[which(as.numeric(PPI_lineages_remaining[,3]) > 0.4 & as.numeric(PPI_lineages_remaining[,4]) > 1e-7),] # 5568
PPI_lineages_medium = PPI_lineages_remaining[which(as.numeric(PPI_lineages_remaining[,3]) > 0.1 & as.numeric(PPI_lineages_remaining[,3]) < 0.4
                                                   & as.numeric(PPI_lineages_remaining[,4]) > 1e-8 ),] # 2417400
PPI_lineages_neutral = PPI_lineages_remaining[which(as.numeric(PPI_lineages_remaining[,3]) > -0.2 & as.numeric(PPI_lineages_remaining [,3]) < 0.2
                                                    & as.numeric(PPI_lineages_remaining[,4]) > 1e-8),] # 4460803
PPI_lineages_low = PPI_lineages_remaining[which(as.numeric(PPI_lineages_remaining[,3]) > -0.5 & as.numeric(PPI_lineages_remaining [,3]) < -0.1
                                                & as.numeric(PPI_lineages_remaining[,4]) > 1e-8),] # 97172
# 1000 lineages with fitness > 0.4; 2000 lineages with fitness > 0.1 & < 0.4; 2000 lineages > -0.2 & < 0.2; 1000 lineages < -0.1
PPI_lineages_final = do.call("rbind", list(PPI_lineages_high[sample(1:nrow(PPI_lineages_high), 1000),],
                                           PPI_lineages_medium[sample(1:nrow(PPI_lineages_medium), 2000),],
                                           PPI_lineages_neutral[sample(1:nrow(PPI_lineages_neutral), 2000),],
                                           PPI_lineages_low[sample(1:nrow(PPI_lineages_low), 1000),]))
PPI_lineages_final_transform = transform_lineages(PPI_lineages_final, 1,2,3,4)

Lineage_plot(PPI_lineages_final_transform, "Working_figure/Figure1/Figure1B_6000_representative_lineages.pdf")


