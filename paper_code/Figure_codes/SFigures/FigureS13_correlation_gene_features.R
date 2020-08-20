# Code for Figure S12 can be found in the code of Main Figure 4C
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")
paper.colors = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
# Calculate a mutability score for each protein
variation_score = read.csv("Datasets_generated_by_preprocessing/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
PP_pair = split_string_vector(variation_score[,1])
protein_unique = unique(as.vector(PP_pair)) # 2083
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
  index_protein = unique(c(which(PP_pair[,1] == protein_unique[i]),
                           which(PP_pair[,2] == protein_unique[i])))
  vScore_protein[i,2] = mean(as.numeric(variation_score[index_protein,3]))
}

##Get GI count for each ORF from Biogrid database
bg = read.delim2(file = "Outsourced_datasets/BIOGRID/BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.5.165.tab2.txt", header = F)
x = bg[2:nrow(bg), c(6,7, 12, 13)]
colnames(x) = c("A", "B", "Experimental System", "Type")
y = x[which(x[,4] == "genetic"),] #all genetic interactions

u = unique(y[,1])
gi = 1:length(u)
names(gi) = u
for(i in u){
  a = y[which(y[,1] == i), 2] #all partners
  gi[i] = length(unique(a)) #unique partners
}
GIcount = gi

gene_feature = as.matrix(read.table("Outsourced_datasets/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39
gene_feature_matched = gene_feature[match(vScore_protein[,1], gene_feature[,1]),]
gene_feature_matched= cbind(gene_feature_matched, GIcount[vScore_protein[,1]], as.numeric(vScore_protein[,2]))
colnames(gene_feature_matched)[40] = "GIdegree"
colnames(gene_feature_matched)[41] = "vScore"
features = c(colnames(gene_feature), "GIdegree", "vScore")
features_chosen = features[c(18, 19, 12, 7, 25, 15, 11, 20, 14, 16, 3, 2, 8, 21, 23, 4, 29, 40 )]
feature_names = c("dN/dS", "Protein disorder", "Number of complexes" ,  "Protein length", "Phenotypic capacitance", 
                  "PPI degree, Y2H" , "Curated phenotypes", "Deletion fitness defect", 'PPI degree, Tap MS', 
                  "Number of domains", "Multifunctionality", "Yeast conservation", "Broad conservation", 
                  "Codon adaptation index",  "Expression level", "PPI degree, Biogrid", "Protein abundance", "Genetic interaction degree" )

num_complex = gene_feature_matched[,c(1, 4, 12, ncol(gene_feature_matched))]
num_complex_new = num_complex[!is.na(num_complex[,3]),]
PPI_degree = as.numeric(num_complex_new[,2])
num_complex_low = num_complex_new[which(PPI_degree >=0 & PPI_degree <= 4),] # 246
num_complex_medium = num_complex_new[which(PPI_degree >=5 & PPI_degree <= 14),] # 274
num_complex_high = num_complex_new[which(as.numeric(num_complex_new[,2]) >= 15),] # 221
cor(as.numeric(num_complex_high[,3]), as.numeric(num_complex_high[,4]), method = 'spearman', use = 'complete.obs')
# 0.21
cor(as.numeric(num_complex_medium[,3]), as.numeric(num_complex_medium[,4]), method = 'spearman', use = 'complete.obs')
# 0.003479496
cor(as.numeric(num_complex_low[,3]), as.numeric(num_complex_low[,4]), method = 'spearman', use = 'complete.obs')
#  0.1011105
complex_count_low = as.numeric(num_complex_low[,3])
complex_count_medium = as.numeric(num_complex_medium[,3])
complex_count_high = as.numeric(num_complex_high[,3])

p_low = table(complex_count_low)
#complex_count_low
#1   2   3 
#227  16   3 
low_count = c(227, 16, 3, 0, 0, 0)
low_label = c('227', '16', '3', '0', '0', '0')
p_medium = table(complex_count_medium)
#1   2   3   4 
#245  24   3   2 0 0
medium_count = c(245, 24, 3, 2, 0, 0)
medium_label = c('245', '24', '3', '2', '0', '0')
p_high = table(complex_count_high)
#1   2   3   4   5   6 
#139  60  13   2   6   1 
low_count + medium_count

high_count = c(139, 60, 13, 2, 6, 1)
high_label = c('139', '60', '13', '2', '6', '1')
sum(high_count)

matrix_complex = rbind(low_count, medium_count, high_count)
all_count = low_count + medium_count  + high_count
all_count # 611 100  19   4   6   1
library(RColorBrewer)
col_chosen = paper.colors[c(8,6,2)]
pdf('Figures/SFigures/SFigure13/SFigure13A_Histogram_complex_count_protein.pdf', width = 5, height = 5)    
barCenter = barplot(matrix_complex, horiz=F, beside=F, ylim=c(0,800), #ylab="Number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6), 
                    col= col_chosen,  ylab = "Number of proteins", border=NA, cex.axis = 0.6)
legend("topright", legend=c("0-4", "5-14", "> 15"), 
       fill=col_chosen[c(1,2,3)], cex = 0.6, bty="n", border=FALSE)
text(x= barCenter, y = all_count + 15, labels = low_label, 
     cex=0.6, xpd = TRUE, col= col_chosen[1]) 
text(x= barCenter, y = all_count + 40, labels = medium_label, 
     cex=0.6, xpd = TRUE, col= col_chosen[2]) 
text(x= barCenter, y = all_count + 65, labels = high_label, 
     cex=0.6, xpd = TRUE, col= col_chosen[3]) 
text(x= barCenter, y = -40, labels = as.character(1:6), xpd = TRUE, cex = 0.6)
text(median(barCenter), y = -100, labels = "Number of complexes per protein", xpd = TRUE)
dev.off()

####### Check the number of members for each complex
complex_GI = csvReader_T('Outsourced_datasets/Costanzo_2016_Science_Data_File_S12_protein_complex_standard.csv')
gene_ppiseq = as.vector(num_complex_new[,1])
length(gene_ppiseq) # 741
index = 0
intersect_count = 0
protein = 0
#component = unlist(strsplit(as.character(complex_GI[537,6]), "; "))
for (i in 1:nrow(complex_GI)){
  component = unlist(strsplit(as.character(complex_GI[i,6]), "; "))
  protein = c(protein, component)
  overlap = length(intersect(component, gene_ppiseq))
  if(overlap >= 1){
    index = c(index, i)
    intersect_count = c(intersect_count, overlap)
  }
}
index_new = index[2:length(index)]
protein = unique(protein[2:length(protein)]) # 2317
overlap = intersect(protein, gene_ppiseq) # 740

gene_ppiseq[which(!gene_ppiseq %in% overlap)] # 'YDR005C' not in the dataframe by which we calculate the R.

num_complex_new[which(num_complex_new == 'YDR005C'),]
# X                   PPI.degree                Number.of.complexes              vScore 
# "YDR005C"               "  8"               "  1"                  "1.589438755" 

intersect_count_new = intersect_count[2:length(intersect_count)]
complex_GI_ppiseq = complex_GI[index_new,] # 357
complex_GI_ppiseq = cbind(complex_GI_ppiseq, intersect_count_new)


complex_GI_ppiseq_two = complex_GI_ppiseq[]

complex_member_count = as.numeric(complex_GI_ppiseq[,3])
complex_ppiseq_gene_count = as.numeric(complex_GI_ppiseq[,8])


pdf('Figures/SFigures/SFigure13/SFigure13B_Histogram_number_members_Protein_complex.pdf', width = 5, height = 5)
par(mar = c(5,4,1,1)) 
hist(complex_member_count, breaks =seq(1,55, by = 1), col = paper.colors[2],
     xlab = 'Number of genes per complex', ylab = 'Number of complexes', main = NA)
dev.off()

############################################################
## Only consider these proteins that participates in more than two complexes
num_complex_new_two = num_complex_new[which(as.numeric(num_complex_new[,3]) >= 2),]
gene_ppiseq = as.vector(num_complex_new_two[,1])
length(gene_ppiseq) # 130
index = 0
intersect_count = 0
protein = 0
#component = unlist(strsplit(as.character(complex_GI[537,6]), "; "))
for (i in 1:nrow(complex_GI)){
  component = unlist(strsplit(as.character(complex_GI[i,6]), "; "))
  protein = c(protein, component)
  overlap = length(intersect(component, gene_ppiseq))
  if(overlap >= 1){
    index = c(index, i)
    intersect_count = c(intersect_count, overlap)
  }
}
index_new = index[2:length(index)] # 125
protein = unique(protein[2:length(protein)]) # 2317
overlap = intersect(protein, gene_ppiseq) # 130

intersect_count_new = intersect_count[2:length(intersect_count)]
complex_GI_ppiseq = complex_GI[index_new,] # 125
complex_GI_ppiseq = cbind(complex_GI_ppiseq, intersect_count_new)

complex_ppiseq_gene_count = as.numeric(complex_GI_ppiseq[,8]) # 274


pdf('Figures/SFigures/SFigure13/SFigure13C_Histogram_number_PPiSeq_genes_per_protein_complex.pdf', width = 5, height = 5)
par(mar = c(5,4,1,1)) 
hist(complex_ppiseq_gene_count, breaks =seq(1,15, by = 1), col = paper.colors[2], ylim = c(0, 100),
     xlab = 'Number of genes in PPiSeq per complex', ylab = 'Number of complexes', main = NA)
dev.off()


