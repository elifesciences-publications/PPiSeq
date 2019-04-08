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

################# Figure 2D Comparison among PPiSeq, DHFR_PCA, and BioGRID. 
### Make a venn diagram to show the overlap the current data with PCA and BIOGRID and fitness distributions for different parts
# Generate BIOGRID databse that not contain PCA PPIs
setwd("~/Dropbox/PPiSeq_02/")
DMSO_mean = csvReader_T("Paper_data/DMSO_mean_fitness_positive.csv") # 1459163
yeast_PPI = read.delim("Working_data/BIOGRID-ORGANISM-3.5.165.tab2/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.165.tab2.txt", header = T)
pos_PPI= yeast_PPI[which(yeast_PPI$Experimental.System.Type == "physical"), c(6:7, 12:14)] # 164992
nrow(pos_PPI[which(pos_PPI[,5] == 'Tarassov K (2008)'),]) # 2616
BIOGRID_PPI = pos_PPI[which(pos_PPI[,5] != 'Tarassov K (2008)'),] #remove mDHFR-PCA, 162376
BIOGRID_PPI_unique = unique(paste(BIOGRID_PPI[,1], BIOGRID_PPI[,2], sep = "_")) # 116665

## Obtain promiscuous proteins from PPiSeq data
fragment_select = dataFrameReader_T("Working_data/Promiscuous_PPIs/DMSO_Promiscuous_proteins.csv")
fragment_protein = unique(as.vector(split_string_vector(fragment_select[,1])))

## Remove PPIs that contain promiscuous proteins from BioGRID
promiscuous_BIOGRID_select = rep(0, length(BIOGRID_PPI_unique))
for(i in 1: length(promiscuous_BIOGRID_select)){
  PPI = split_string(BIOGRID_PPI_unique[i])
  if(PPI[1] %in% fragment_protein | PPI[2] %in% fragment_protein){
    promiscuous_BIOGRID_select[i] = 1
  }
}
length(which(promiscuous_BIOGRID_select == 1)) #113
BIOGRID_PPI_unique_filter = BIOGRID_PPI_unique[which(promiscuous_BIOGRID_select != 1)] #116552

# Here we only consider comparison of one direction PPIs among three datasets
BIOGRID_PPI_unique_deduplicate = mark_duplicates_fast(BIOGRID_PPI_unique_filter) # 106837

## Check which PPIs in BIOGRID after removing promisicuous PPIs that were covered by PPiseq
# match_both_direction: output the PPIs in the first input matrix that were matched two orientations of the second input
BIOGRID_PPiseq = match_both_direction(BIOGRID_PPI_unique_deduplicate, as.character(DMSO_mean[,1]))# 11480
BIOGRID_PPI_uncover = BIOGRID_PPI_unique_deduplicate[which(!BIOGRID_PPI_unique_deduplicate[,1] %in% BIOGRID_PPiseq[,1]),] # 95357

# Check which PPI in PCA database after removing promiscuous PPIS that were covered by PPiseq 
PPI_paper= read.delim(file = "Working_data/PPI_set_PCA_science.txt", sep=" ")
PPI_PCA= paste(PPI_paper[,1], PPI_paper[,3], sep="_") # 2770
# Remove PPIs that contain promiscuous proteins
promiscuous_PCA_select = rep(0, length(PPI_PCA))
for(i in 1: length(promiscuous_PCA_select)){
  PPI = split_string(PPI_PCA[i])
  if(PPI[1] %in% fragment_protein | PPI[2] %in% fragment_protein){
    promiscuous_PCA_select[i] = 1
  }
}
length(which(promiscuous_PCA_select == 1)) #35
PCA_PPI_unique_filter = PPI_PCA[which(promiscuous_PCA_select != 1)] #2735
PCA_PPI_unique_deduplicate = mark_duplicates_fast(PCA_PPI_unique_filter) # 2735 no duplicated PPIs (two orientations)
PCA_PPI_matrix = cbind(PCA_PPI_unique_filter, 1:length(PCA_PPI_unique_filter)) # add line number to make it can be input of match_both_direction
# Get PCA data covered by PPiseq
PCA_PPiseq = match_both_direction(PCA_PPI_matrix, as.character(DMSO_mean[,1])) #2310
PCA_PPiseq_uncover = PCA_PPI_unique_filter[which(!PCA_PPI_unique_filter %in% PCA_PPiseq[,1])] # 425

####### Check PPIs in DMSO environment
#Get one direction of positive PPIs called by PPiSeq
pos_PPI = DMSO_mean[which(DMSO_mean[,7] == "1"),]
#Remove all the control strains
pos_RRS =pos_PPI[grep("Neg_PPI", pos_PPI[,1]),1] #29
pos_PRS =pos_PPI[grep("Pos_PPI", pos_PPI[,1]),1] #3
pos_DHFR =pos_PPI[grep("positive_DHFR", pos_PPI[,1]),1] # 1
pos_PPI_real = pos_PPI[which(!pos_PPI[,1] %in% c(pos_RRS, pos_PRS, pos_DHFR)),1] # 5178
pos_PPI_unique_filter = mark_duplicates_fast(pos_PPI_real) # 4786

BIOGRID_PPiseq_overlap = match_both_direction(BIOGRID_PPiseq, pos_PPI_unique_filter[,1])#656
BIOGRID_PPiseq_PCA_overlap = match_both_direction(BIOGRID_PPiseq_overlap, as.character(PCA_PPiseq[,1])) #456
PCA_BIOGRID_overlap = match_both_direction(PCA_PPiseq, BIOGRID_PPiseq[,1]) #652
PCA_PPiseq_overlap = match_both_direction(PCA_PPiseq, pos_PPI_unique_filter[,1]) #1198

area1 = nrow(pos_PPI_unique_filter) #4786
area2 = nrow(PCA_PPiseq) # 2310
area3 = nrow(BIOGRID_PPiseq) # 11480
n12 = nrow(PCA_PPiseq_overlap) # 1198
n23 = nrow(PCA_BIOGRID_overlap) # 652
n13 = nrow(BIOGRID_PPiseq_overlap) # 656
n123 = nrow(BIOGRID_PPiseq_PCA_overlap) #456
library(VennDiagram)
pdf("Working_figure/Figure2/Figure2D_draft_Venn diagram of overlapped PPIs among PPiseq, PCA, BIOGRID.pdf", width = 5, height =5)
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c("PPiSeq", "mDHFR-PCA", "BIOGRID data excluding mDHFR-PCA"), euler.d = TRUE,
                 scaled = TRUE, col = apple_colors[c(5,3,6)], fill = apple_colors[c(5,3,6)], cex = 1.5)
dev.off()

# Also include BIOGRID Uncovered, PCA Uncovered, and the intersection between them
a = nrow(BIOGRID_PPI_uncover) # 95357
b = length(PCA_PPiseq_uncover) # 425
c = length(intersect(BIOGRID_PPI_uncover[,1], PCA_PPiseq_uncover)) # 97
a-c # 95260
b-c # 328
nrow(DMSO_mean) # 1459163
# Manually draw the figure in keynote

####### Check PPIs in all environment
#Get one direction of positive PPIs called by PPiSeq
pos_PPI = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv") # 13829
pos_PPI_unique_filter = pos_PPI[,1] # Already remove control strains, duplicate PPIs

BIOGRID_PPiseq_overlap = match_both_direction(BIOGRID_PPiseq, pos_PPI_unique_filter[,1])#867
BIOGRID_PPiseq_PCA_overlap = match_both_direction(BIOGRID_PPiseq_overlap, as.character(PCA_PPiseq[,1])) #528
PCA_BIOGRID_overlap = match_both_direction(PCA_PPiseq, BIOGRID_PPiseq[,1]) #652
PCA_PPiseq_overlap = match_both_direction(PCA_PPiseq, pos_PPI_unique_filter[,1]) #1507

area1 = nrow(pos_PPI_unique_filter) # 13829
area2 = nrow(PCA_PPiseq) # 2310
area3 = nrow(BIOGRID_PPiseq) # 11480
n12 = nrow(PCA_PPiseq_overlap) # 1507
n23 = nrow(PCA_BIOGRID_overlap) # 652
n13 = nrow(BIOGRID_PPiseq_overlap) # 867
n123 = nrow(BIOGRID_PPiseq_PCA_overlap) #528
library(VennDiagram)
pdf("Working_figure/Figure2/Figure2D_draft_all_Venn diagram of overlapped PPIs among PPiseq, PCA, BIOGRID.pdf", width = 5, height =5)
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c("PPiSeq", "mDHFR-PCA", "BIOGRID data excluding mDHFR-PCA"), euler.d = TRUE,
                 scaled = TRUE, col = apple_colors[c(5,3,6)], fill = apple_colors[c(5,3,6)], cex = 1.5)
dev.off()

# Also include BIOGRID Uncovered, PCA Uncovered, and the intersection between them
a = nrow(BIOGRID_PPI_uncover) # 95357
b = length(PCA_PPiseq_uncover) # 425
c = length(intersect(BIOGRID_PPI_uncover[,1], PCA_PPiseq_uncover)) # 97
a-c # 95260
b-c # 328
nrow(DMSO_mean) # 1459163
