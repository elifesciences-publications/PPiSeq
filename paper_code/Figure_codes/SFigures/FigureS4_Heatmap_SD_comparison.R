#### Make a heatmap to show the overlap coefficient across different datasets
source("~/Desktop/PPiSeq_additional_data/function.R") # Load commonly used functions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

setwd('~/Desktop/PPiSeq_additional_data/')
y2h = as.matrix(read.delim('Outsourced_datasets/CCSB-Y11.txt', header = F))
pca = as.matrix(read.delim('Outsourced_datasets/PPI_set_PCA_science.txt', header = T, sep = ' '))
prs = csvReader_T('SFigure4_related_data/Positive_reference_set_SD.csv') 
prs = mark_duplicates_fast(prs)

ppiseq = csvReader_T('Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge.csv')
ppiseq_sd = ppiseq[which(ppiseq[,3] == '1'),] # 4704
sd_multiple = csvReader_T('PPI_mean_fitness_calling_files/SD_merge_mean_fitness_positive.csv')
ppi_multiple = sd_multiple[,1]

## Generate PPIs for two datasets
y2h_ppi = cbind(paste(y2h[,1], y2h[,2], sep = '_'), y2h)
pca_ppi = cbind(paste(pca[,1], pca[,3], sep = '_'), pca)

## PPIs within PPiSeq search space
y2h_ppi_in = match_both_direction(y2h_ppi, ppi_multiple)
pca_ppi_in = match_both_direction(pca_ppi, ppi_multiple)
prs_ppi_in = match_both_direction(prs, ppi_multiple)

## Only consider one orientation
y2h_ppi_in_unique = mark_duplicates(y2h_ppi_in[,1])[,1]
pca_ppi_in_unique = mark_duplicates(pca_ppi_in[,1])[,1]

# Check the overlap between SD replicates and PCA
SD1_mean = csvReader_T("PPI_mean_fitness_calling_files/SD_mean_fitness_positive.csv")
SD2_mean = csvReader_T("PPI_mean_fitness_calling_files/SD2_mean_fitness_positive.csv")
SD1_pos = SD1_mean[which(SD1_mean[,ncol(SD1_mean)] == "1"),]
SD2_pos = SD2_mean[which(SD2_mean[,ncol(SD2_mean)] == "1"),]
SD1_pos_unique = mark_duplicates(SD1_pos[,1])
SD2_pos_unique = mark_duplicates(SD2_pos[,1])

## Unique PPI sets: ppiseq_sd, SD1_pos_unique, SD2_pos_unique, prs_ppi_in, pca_ppi_in_unique, y2h_ppi_in_unique, 
ppi_set = list('SD_merge' = ppiseq_sd[,1], 'SD1'= SD1_pos_unique, 'SD2'= SD2_pos_unique,
               'PCA' = pca_ppi_in_unique, 'PRS' = prs_ppi_in, 'Y2H' = y2h_ppi_in_unique)

# A function to calculate overlap coefficiency (consider two orientations as the one PPI)
calculate_overlap = function(a, b){
  a_sep = split_string_vector(a)
  a_r = paste(a_sep[,2], a_sep[,1], sep = '_')
  b_sep = split_string_vector(b)
  b_r = paste(b_sep[,2], b_sep[,1], sep = '_')
  a_all = c(a, a_r)
  b_all = c(b, b_r)
  if (length(a_all) <= length(b_all)){
    overlap = length(which(a_all %in% b_all))/2
  }else{
    overlap = length(which(b_all %in% a_all))/2
  }
  min_size = min(c(length(a), length(b)))
  overlap/min_size
  return(overlap/min_size)
}


overlap_coeff = matrix(0,6,6)
for (i in 1:6){
  for (j in 1:6){
    if (i == j){
      overlap_coeff[i,j] = NA
    }else{
      overlap_coeff[i,j] = calculate_overlap(ppi_set[[i]], ppi_set[[j]])
    }
   
  }
}

label = c('SD_merge', 'SD1', 'SD2','PCA', 'PRS', 'Y2H')
colnames(overlap_coeff) = label
#rownames(overlap_coeff) = c('SD_merge', 'SD1', 'SD2', 'PCA', 'PRS', 'Y2H')

write.csv(overlap_coeff, 'SFigure4_related_data/Overlap_coefficient_various_dataset.csv', row.names = F)

PPI_heatmap = dataFrameReader_T("SFigure4_related_data/Overlap_coefficient_various_dataset.csv")

#install.packages("pheatmap")
library(pheatmap)
library(RColorBrewer)
heatmap_matrix = PPI_heatmap
#fitness_all = unique(as.vector(as.matrix(heatmap_matrix)))
bk = seq(0, 1, by = 0.05)

col_chosen = c(apple_colors[c(5,2,7)])

my_palette = c(colorRampPalette(col_chosen)(length(bk)))
fit_heatmap = pheatmap(heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames=FALSE,
                       show_colnames=T, col = my_palette,  display_numbers = TRUE,
                        treeheight_row = 0, angle_col = 45)

save_pheatmap_pdf <- function(x, filename, width=4.5, height=5) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(fit_heatmap, "Figures/SFigures/SFigure4/SFigure4A_Overlap_coefficient_environment_heatmap.pdf")

for (i in 1:length(ppi_set)){
  if(is.vector(ppi_set[[i]])){
    print(length(ppi_set[[i]]))
  }
  else{
    print(nrow(ppi_set[[i]]))
  }
}
# PPI number:
#SD_merge: 4704; SD1: 4913, SD2: 4200, PCA: 2512; PRS: 493; Y2H: 304


# pca_ppi_in
SD1_PCA = match_both_direction(SD1_pos_unique, pca_ppi_in[,1])
SD1_PCA_SD2 = match_both_direction(SD1_PCA, SD2_pos_unique[,1])
SD2_SD1 = match_both_direction(SD2_pos_unique, SD1_pos_unique[,1])
SD2_PCA = match_both_direction(SD2_pos_unique, pca_ppi_in[,1])
area1 = nrow(SD1_pos_unique) #4913
area2 = nrow(SD2_pos_unique) #4200
area3 = nrow(pca_ppi_in) # 2512
n12 = nrow(SD2_SD1) # 3314
n23 = nrow(SD2_PCA) # 1139
n13 = nrow(SD1_PCA) #1237
n123 = nrow(SD1_PCA_SD2) #1073
library(VennDiagram)
pdf("Figures/SFigures/SFigure4/SFigure4B_SD1_SD2_PCA_overlap.pdf", width = 5, height =5)
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c("SD_1", "SD_2", "PCA"), euler.d = TRUE,
                 scaled = TRUE, col = apple_colors[c(5,3,6)], fill = apple_colors[c(5,3,6)], cex = 1.5)
dev.off()

SD1_SD2_overlap_PCA = n123/min(n12,area3) # 0.4271497
SD1_nonoverlap_PCA = 164/min(c(1435 + 164, area3)) # 0.1025641
SD2_nonoverlap_PCA = 66/min(c(820 + 66, area3)) # 0.0744921
SD1_SD2_nonoverlap_PCA = (164+66)/min(c(1435 + 164 + 820+66, area3)) # 0.09255533
SD1_PCA = n13/min(c(area1, area3)) # 0.4924363
SD2_PCA = n23/min(c(area2, area3)) # 0.4534236


## prs_ppi_in
SD1_PCA = match_both_direction(SD1_pos_unique, prs_ppi_in[,1])
SD1_PCA_SD2 = match_both_direction(SD1_PCA, SD2_pos_unique[,1])
SD2_SD1 = match_both_direction(SD2_pos_unique, SD1_pos_unique[,1])
SD2_PCA = match_both_direction(SD2_pos_unique, prs_ppi_in[,1])
area1 = nrow(SD1_pos_unique) #4913
area2 = nrow(SD2_pos_unique) #4200
area3 = nrow(prs_ppi_in) # 493
n12 = nrow(SD2_SD1) # 3314
n23 = nrow(SD2_PCA) # 202
n13 = nrow(SD1_PCA) # 218
n123 = nrow(SD1_PCA_SD2) #193
library(VennDiagram)
pdf("Figures/SFigures/SFigure4/SFigure4B_SD1_SD2_prs_overlap.pdf", width = 5, height =5)
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c("SD_1", "SD_2", "PRS"), euler.d = TRUE,
                 scaled = TRUE, col = apple_colors[c(5,3,6)], fill = apple_colors[c(5,3,6)], cex = 1.5)
dev.off()

SD1_SD2_overlap_PRS = n123/min(n12,area3) 
SD1_SD2_overlap_PRS # 0.3914807
SD1_nonoverlap_PRS = 25/min(c(1574 + 25, area3)) 
SD1_nonoverlap_PRS # 0.05070994
SD2_nonoverlap_PRS = 9/min(c(877 + 9, area3)) 
SD2_nonoverlap_PRS # 0.01825558
SD1_SD2_nonoverlap_PRS = (25+9)/min(c(1574 + 25 + 877+9, area3)) 
SD1_SD2_nonoverlap_PRS # 0.06896552
SD1_PRS = n13/min(c(area1, area3)) 
SD1_PRS # 0.4421907
SD2_PRS = n23/min(c(area2, area3)) 
SD2_PRS # 0.4097363

## y2h_ppi_in

SD1_PCA = match_both_direction(SD1_pos_unique, y2h_ppi_in[,1])
SD1_PCA_SD2 = match_both_direction(SD1_PCA, SD2_pos_unique[,1])
SD2_SD1 = match_both_direction(SD2_pos_unique, SD1_pos_unique[,1])
SD2_PCA = match_both_direction(SD2_pos_unique, y2h_ppi_in[,1])
area1 = nrow(SD1_pos_unique) #4913
area2 = nrow(SD2_pos_unique) #4200
area3 = nrow(y2h_ppi_in) # 304
n12 = nrow(SD2_SD1) # 3314
n23 = nrow(SD2_PCA) # 38
n13 = nrow(SD1_PCA) # 44
n123 = nrow(SD1_PCA_SD2) #36
library(VennDiagram)
pdf("Figures/SFigures/SFigure4/SFigure4B_SD1_SD2_y2h_overlap.pdf", width = 5, height =5)
draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, category = c("SD_1", "SD_2", "Y2H"), euler.d = TRUE,
                 scaled = TRUE, col = apple_colors[c(5,3,6)], fill = apple_colors[c(5,3,6)], cex = 1.5)
dev.off()

SD1_SD2_overlap_Y2H = n123/min(n12,area3) 
SD1_SD2_overlap_Y2H # 0.1184211
SD1_nonoverlap_Y2H = 8/min(c(1591 + 8, area3)) 
SD1_nonoverlap_Y2H # 0.02631579
SD2_nonoverlap_Y2H = 2/min(c(884 + 2, area3)) 
SD2_nonoverlap_Y2H # 0.006578947
SD1_SD2_nonoverlap_Y2H = (8+2)/min(c(1591 + 8 + 884 + 2, area3)) 
SD1_SD2_nonoverlap_Y2H # 0.03289474
SD1_Y2H = n13/min(c(area1, area3)) 
SD1_Y2H # 0.1447368
SD2_Y2H = n23/min(c(area2, area3)) 
SD2_Y2H # 0.125

barplot = c(SD1_SD2_overlap_PCA, SD1_SD2_nonoverlap_PCA, 
            SD1_SD2_overlap_PRS, SD1_SD2_nonoverlap_PRS, 
            SD1_SD2_overlap_Y2H, SD1_SD2_nonoverlap_Y2H)


library(RColorBrewer)

col_chosen = apple_colors[c(5,3)]
pdf("Figures/SFigures/SFigure4/SFigure4C_Bar_plot_overlap_coefficient.pdf", width= 5, height=5)
par(mar = c(5,4,1,1)) 
barCenter = barplot(as.numeric(barplot), horiz=F, beside=F, ylim=c(0,0.5), ylab="Overlap coefficient",
                    space= c(0.3, 0.1, 0.3,  0.1, 0.3, 0.1),axisnames=F, border=NA, 
                    col = col_chosen, cex.axis=1)
legend("topright", legend=c("Overlapped PPIs in two SD replicates", "Non-overlapped PPIs in two SD replicates"), 
       fill=col_chosen[c(1,2)], cex = 0.8, bty="n", border=FALSE)

dev.off()


nrow(y2h_ppi_in)
nrow(prs_ppi_in)

### Check the fitness values of overlapped and non-overlapped PPIs
ppiseq = csvReader_T('Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter.csv')
fit_ppiseq = csvReader_T('Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv')

fit_ppiseq_sd = fit_ppiseq[which(ppiseq[,3] == "1"), ] # 3921
fit_ppiseq_sd_pca = match_both_direction(fit_ppiseq_sd, pca_ppi_in[,1]) # 1114
fit_ppiseq_sd_non_pca = fit_ppiseq_sd[which(!fit_ppiseq_sd[,1] %in% fit_ppiseq_sd_pca[,1]),] # 2807

fit_ppiseq_pca = as.numeric(fit_ppiseq_sd_pca[,4])
fit_ppiseq_nonpca = as.numeric(fit_ppiseq_sd_non_pca[,4])

t.test(fit_ppiseq_pca, fit_ppiseq_nonpca, 'greater')
p1 <- hist(fit_ppiseq_nonpca, breaks = seq(0.1, 1.5, by = 0.02))                     
p2 <- hist(fit_ppiseq_pca,breaks = seq(0.1, 1.5, by = 0.02))
pdf('Figures/SFigures/SFigure4/SFigure4D_Histogram_fitness_overlap_PCA.pdf', width = 5, height = 5)
par(mar = c(5,4,1,1)) 
plot(p1, col=alpha(apple_colors[5], alpha = 0.8), xlim = c(0, 1),
     xlab = 'Fitness', ylab = "Number of PPIs", main = NULL, border = 'NA')  # first histogram
plot(p2, col=alpha(apple_colors[7], alpha = 0.8), border = 'NA', add=T)  # second
dev.off()