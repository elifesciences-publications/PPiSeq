setwd("~/Dropbox/PPiSeq_02/")

apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

#Make a combined data matrix with variation score
load("Working_data/Positive_PPI_environment/variation_score.Rfile")
load("Working_data/Positive_PPI_environment/vScore_protein.Rfile")
load("Working_data/BIOGRID-ORGANISM-3.5.165.tab2/GIcount.Rfile")
gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39
gene_feature_matched = gene_feature[match(vScore_protein[,1], gene_feature[,1]),]
gene_feature_matched= cbind(gene_feature_matched, GIcount[vScore_protein[,1]], as.numeric(vScore_protein[,2]))
colnames(gene_feature_matched)[40] = "GIdegree"
colnames(gene_feature_matched)[41] = "vScore"
features = c(colnames(gene_feature), "GIdegree", "vScore")
features_chosen = features[c(18, 19, 12, 7, 25, 15, 11, 20, 14, 16, 3, 2, 8, 21, 23, 4, 29, 40 )]
feature_names = c("dN/dS", "Protein disorder", "Number of complexes" ,  "Protein length", "Phenotypic capacitance", 
                  "PPI degree, Y2H" , "Curated phenotypes", "Single mutant fitness defect", 'PPI degree, Tap MS', 
                  "Number of domains", "Multifunctionality", "Yeast conservation", "Broad conservation", 
                  "Codon adaptation index",  "Expression level", "PPI degree, Biogrid", "Protein abundance", "GI degree", "vScore" )
degree_bins = matrix(c(0, 4, 5, 14, 15, 1000), 2, 3)
n = 1:ncol(degree_bins)
for(i in 1:ncol(degree_bins)){
  n[i]=paste(degree_bins[1,i], degree_bins[2,i], sep = "-")
}
degree_bin_names = n

#get the real correlation values
vScore = as.numeric(gene_feature_matched[,41])
degree = as.numeric(gene_feature_matched[,4])
f_cor = matrix(0, length(features_chosen), ncol(degree_bins))
colnames(f_cor) = n
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
  rownames(f_cor_b) = features_chosen
  for(j in 1:ncol(degree_bins)){
    #n = length(which(degree_b >= degree_bins[1,j] & degree_b <= degree_bins[2,j]))
    #delta = 1.96/sqrt(n - 3)
    #print(length(which(degree > degree_bins[1,j] & degree < degree_bins[2,j])))
    for(i in 1:length(features_chosen)){
      x = as.numeric(g_b[which(degree_b > degree_bins[1,j] & degree_b < degree_bins[2,j]), features_chosen[i]])
      y = vScore[which(degree_b > degree_bins[1,j] & degree_b < degree_bins[2,j])]
      f_cor_b[i,j] = cor(x,y, method = 'spearman', use = 'complete.obs' )
    }
  }
  f_cor_boot[[k]] = f_cor_b
}

#make a mena and SD matrix from bootstrapped data
sd = f_cor_boot[[1]]
mean = f_cor_boot[[1]]
sem = f_cor_boot[[1]]
for(i in 1:nrow(sd)){
  for(j in 1:ncol(sd)){
    x = sapply(f_cor_boot, function(a) a[i,j])
    sd[i,j] = sd(x)
    mean[i,j] = mean(x)
    sem[i,j] = sd(x)/sqrt(1000)
  }
}

