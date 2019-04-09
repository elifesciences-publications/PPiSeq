############# This script contains the supplementary figures around Figure1
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

##### Create initial barcode distribution from G0 data of DMSO environment
setwd("~/Dropbox/PPiSeq_02/")
lineage = dataFrameReader_T("Working_data/DMSO_all_lineage_counts/lineage_trajectories.csv")
known_DBC= dataFrameReader_T("Working_data/PPI_barcodes/PPiseq_all_barcodes.csv")
nrow(known_DBC) # 7741081
lineage_T0 = lineage[which(lineage[,2] != 0),]
nrow(lineage_T0) # 6185226
length(which(known_DBC[,3] %in% lineage_T0[,1])) # 6048285

counts = lineage_T0[,2] # max:45389; min:0; 33.93565; 6185226

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

ggsave("Working_figure/FigureS1/FigureS1A_Distribution_initial_barcode_counts.pdf", width = 5, height = 5)

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
ggsave("Working_figure/FigureS1/FigureS1A_Distribution_initial_barcode_counts_non_log.pdf", width = 5, height = 5)