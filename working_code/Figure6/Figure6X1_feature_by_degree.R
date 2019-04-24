
setwd("~/Dropbox/PPiSeq_02/")


#Make a combined data matrix with variation score
load("Working_data/Positive_PPI_environment/variation_score.Rfile")
load("Working_data/Positive_PPI_environment/vScore_protein.Rfile")
gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39
gene_feature_matched = gene_feature[match(vScore_protein[,1], gene_feature[,1]),]
features = colnames(gene_feature)
features_chosen = features[2:(length(features)-7)]
vScore = as.numeric(vScore_protein[,2])
degree = as.numeric(gene_feature_matched[,4])
degree_bins = matrix(c(0,5,6,10,11,20,21,30,31,40,41,50,51,1000), 2, 7)


f_cor = matrix(0, length(features_chosen), ncol(degree_bins))
n = 1:ncol(degree_bins)
for(i in 1:ncol(degree_bins)){
  n[i]=paste(degree_bins[1,i], degree_bins[2,i], sep = "-")
}
colnames(f_cor) = n
rownames(f_cor) = features_chosen
for(j in 1:ncol(degree_bins)){
  print(length(which(degree > degree_bins[1,j] & degree < degree_bins[2,j])))
  for(i in 1:length(features_chosen)){
    x = as.numeric(gene_feature_matched[which(degree > degree_bins[1,j] & degree < degree_bins[2,j]), features_chosen[i]])
    y = vScore[which(degree > degree_bins[1,j] & degree < degree_bins[2,j])]
    f_cor[i,j] = cor(x,y, method = 'spearman', use = 'complete.obs' )
  }
}

#x = as.data.frame(f_cor)
#x$features = features_chosen

for(i in 1:ncol(f_cor)){
  o = order(f_cor[,i], decreasing = F)
  par(mar= c(5,5,1,1))
  pdf(file = paste("Working_figure/Figure6/Figure6X1_feature_by_degree", colnames(f_cor)[i], ".pdf", sep = ""))
  barCenter = barplot(f_cor[o, i], horiz=T, beside=F, xlim=c(-0.6, 0.6), 
                      space= rep(0.4, 31),
                      xlab="Spearman correlation with protein variation score",
                      axisnames=F, border=NA, cex.axis=0.7, col = apple_colors[5], main = paste("degree = ", colnames(f_cor)[i]))
  text(-0.6, barCenter, labels= rownames(f_cor)[o], cex = 0.5, xpd = TRUE)
  dev.off()
}



