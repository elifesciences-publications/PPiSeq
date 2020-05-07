setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used functions
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
degree_bins = matrix(c(0, 4, 5, 14, 15, 1000, 0, 1000), 2, 4)
n = 1:ncol(degree_bins)
for(i in 1:ncol(degree_bins)){
  n[i]=paste(degree_bins[1,i], degree_bins[2,i], sep = "-")
}
degree_bin_names = n

#get the real correlation values
vScore = as.numeric(gene_feature_matched[,41])
degree = as.numeric(gene_feature_matched[,4])
f_cor = matrix(0, length(features_chosen), ncol(degree_bins))
colnames(f_cor) =degree_bin_names
rownames(f_cor) = feature_names
for(j in 1:ncol(degree_bins)){
  n = length(which(degree >= degree_bins[1,j] & degree <= degree_bins[2,j]))
  delta = 1.96/sqrt(n - 3)
  print(length(which(degree > degree_bins[1,j] & degree < degree_bins[2,j])))
  for(i in 1:length(features_chosen)){
    x = as.numeric(gene_feature_matched[which(degree > degree_bins[1,j] & degree < degree_bins[2,j]), features_chosen[i]])
    y = vScore[which(degree > degree_bins[1,j] & degree < degree_bins[2,j])]
    f_cor[i,j] = cor(x,y, method = 'spearman', use = 'complete.obs' )
  }
}

#Get bootstrap data (1000x)
f_cor_boot = list()
g = gene_feature_matched
for(k in 1:1000){
  b = sample(1:nrow(g), nrow(g), replace = T) #the new bootstrapped indexes
  g_b = g[b,] #the new bootstrapped gene features
  vScore_b = as.numeric(g_b[,41])
  degree_b = as.numeric(g_b[,4])
  f_cor_b = matrix(0, length(features_chosen), ncol(degree_bins))
  colnames(f_cor_b) = degree_bin_names
  rownames(f_cor_b) = feature_names
  for(j in 1:ncol(degree_bins)){
    #n = length(which(degree_b >= degree_bins[1,j] & degree_b <= degree_bins[2,j]))
    #delta = 1.96/sqrt(n - 3)
    #print(length(which(degree > degree_bins[1,j] & degree < degree_bins[2,j])))
    for(i in 1:length(features_chosen)){
      x = as.numeric(g_b[which(degree_b > degree_bins[1,j] & degree_b < degree_bins[2,j]), features_chosen[i]])
      y = vScore_b[which(degree_b > degree_bins[1,j] & degree_b < degree_bins[2,j])]
      f_cor_b[i,j] = cor(x,y, method = 'spearman', use = 'complete.obs' )
    }
  }
  f_cor_boot[[k]] = f_cor_b
}

#make a mena and SD matrix from bootstrapped data
b.sd = f_cor_boot[[1]]
b.mean = f_cor_boot[[1]]
b.sem = f_cor_boot[[1]]
for(i in 1:nrow(b.sd)){
  for(j in 1:ncol(b.sd)){
    x = sapply(f_cor_boot, function(a) a[i,j])
    b.sd[i,j] = sd(x)
    b.mean[i,j] = mean(x)
    b.sem[i,j] = sd(x)/sqrt(1000)
  }
}


#The supplemenatal plot using all data
i = ncol(b.mean)
o = order(b.mean[,i], decreasing = F)

# Generate Figure S9
pdf(file = "Figures/SFigures/SFigure9/SFigure9_Correlation_variability_other_features.pdf")
par(mar= c(5,5,1,1))
barCenter = barplot(b.mean[o, i], horiz=T, beside=F, xlim=c(-0.5, 0.3), 
                    xlab="Spearman correlation with the mutability score",
                    axisnames=F, border=NA, cex.axis=1, col = paper.colors[2])
arrows( b.mean[o, i] - b.sd[o,i], barCenter,  b.mean[o, i] + b.sd[o,i], barCenter, angle = 90, lwd = 1.5, code = 3, length = 0.05)
text(-0.5, barCenter, labels= rownames(f_cor)[o], cex = 0.8, xpd = TRUE)
dev.off()

#Comparison by degree of some major features(Figure S4C)
fs = c(17, 15, 18, 8,  3,  2, 1)
fc = b.mean[fs, 3:1]
fc_ci_lower = b.mean[fs, 3:1]- b.sd[fs, 3:1] 
fc_ci_upper = b.mean[fs, 3:1]+ b.sd[fs, 3:1] 
pdf(file = "Figures/Figure4/Figure4C_gene_features_PPI_stability.pdf")
par(mar= c(5,5,1,2))
barCenter = barplot(t(fc), horiz=T, beside=T, xlim=c(-0.6, 0.4), 
                    xlab="Correlation with the mutability score",
                    axisnames=F, border=NA, cex.axis=1.2, col = paper.colors[c(2,6,8)])
arrows( t(fc_ci_lower), barCenter, t(fc_ci_upper), barCenter, angle = 90, lwd = 1.5, code = 3, length = 0.05)
text(-0.55, barCenter[2,], labels= rownames(fc), cex = 1.2, xpd = TRUE)
legend(.1, 8, c("0-4", "5-14", "15+"), bty = "n", title = "PPI degree", 
       fill = paper.colors[c(8,6,2)], cex = 1.2)
dev.off()

