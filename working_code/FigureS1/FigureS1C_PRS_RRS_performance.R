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

################# Figure S1C The performance of PRS and RRS. That is comparable with Y2H an PCA.
setwd("~/Dropbox/PPiSeq_02/")
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv") # 1459163
DMSO_pos = DMSO_mean[which(DMSO_mean[,7] == "1"),] # 5211
DMSO_RRS = DMSO_mean[grep("Neg_PPI", DMSO_mean[,1]),1] #97
DMSO_PRS = DMSO_mean[grep("Pos_PPI", DMSO_mean[,1]),1] #108
DMSO_pos_RRS = DMSO_pos[grep("Neg_PPI", DMSO_pos[,1]),1] # 29
DMSO_pos_PRS = DMSO_pos[grep("Pos_PPI", DMSO_pos[,1]),1] # 3
DMSO_pos_RRS_unique = unique(split_string_vector_line(DMSO_pos_RRS)[,1]) # 2
DMSO_pos_PRS_unique = unique(split_string_vector_line(DMSO_pos_PRS)[,1]) # 20
DMSO_RRS_unique = unique(split_string_vector_line(DMSO_RRS)[,1]) # 67
DMSO_PRS_unique = unique(split_string_vector_line(DMSO_PRS)[,1]) # 71
PRS_percent = length(DMSO_pos_PRS_unique)/length(DMSO_PRS_unique)* 100
RRS_percent = length(DMSO_pos_RRS_unique)/length(DMSO_RRS_unique) * 100
# make a barplot to show the percentage 
PRS_percent #28.17
RRS_percent #2.99

pdf("~/Dropbox/PPiSeq_02/Working_figure/FigureS1/FigureS1C_PRS_RRS_barplot.pdf", width = 3, height =4)
par(mar= c(3,4.5,1,1))
name = c("PRS", "RRS")
barCenter = barplot(c(PRS_percent, RRS_percent), horiz=F, beside=F, ylim = c(0,40),
                    ylab= "Fraction positive by PPiSeq (%)", space = c(0.2),
                    col = c(apple_colors[7], apple_colors[5]), axisnames = F,
                    border = NA)
text(x = barCenter+0.2, y = -3, adj = 1, labels = name, xpd = TRUE)
text(x = barCenter, y = c(PRS_percent + 2, RRS_percent + 2), labels = c("20/71", "2/67"))
dev.off()