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

#### Check these PPIs can not be explained by homo-dimer dynamics
setwd("~/Dropbox/PPiSeq_02/")
homo_plot_matrix = dataFrameReader_T("Working_data/homo_dimer/PPI_homo_dimer_correlation_abundance.csv")
self_interacting = dataFrameReader_T("Working_data/homo_dimer/Self_protein_matrix.csv")
self_interacting_high_cor = as.character(self_interacting[which(self_interacting$Mean_cor >= 0.5 & self_interacting$PPI_number >= 5),1])
PPI_special = homo_plot_matrix[which(homo_plot_matrix$Correlation < 0 & homo_plot_matrix$Difference_amount > 0
                                     & as.character(homo_plot_matrix$Homo.dimer) %in% self_interacting_high_cor),] # 62
csvWriter(PPI_special, "Working_data/homo_dimer/PPI_special_not_expression.csv")

