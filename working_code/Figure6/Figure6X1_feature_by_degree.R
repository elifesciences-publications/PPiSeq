
#setwd("~/Dropbox/PPiSeq_02/")

apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

#Make a combined data matrix with variation score
load("Working_data/Positive_PPI_environment/variation_score.Rfile")
load("Working_data/Positive_PPI_environment/vScore_protein.Rfile")
gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39
gene_feature_matched = gene_feature[match(vScore_protein[,1], gene_feature[,1]),]
features = colnames(gene_feature)
features_chosen = features[c(18, 19, 12, 7, 25, 15, 11, 20, 14, 16, 3, 2, 8, 21, 23, 4, 29 )]
feature_names = c("dN/dS", "Protein disorder", "Number of complexes" ,  "Protein length", "Phenotypic capacitance", 
                  "PPI degree, Y2H" , "Curated phenotypes", "Single mutant fitness defect", 'PPI degree, Tap MS', 
                  "Number of domains", "Multifunctionality", "Yeast conservation", "Broad conservation", 
                  "Codon adaptation index",  "Expression level", "PPI degree, Biogrid", "Protein abundance" )
vScore = as.numeric(vScore_protein[,2])
degree = as.numeric(gene_feature_matched[,4])
degree_bins = matrix(c(0,5,6,10,11,20,21,30,31,40,41,50,51,1000, 0, 3, 4, 14, 15, 1000, 0 , 1000), 2, 11)


f_cor = matrix(0, length(features_chosen), ncol(degree_bins))
n = 1:ncol(degree_bins)
for(i in 1:ncol(degree_bins)){
  n[i]=paste(degree_bins[1,i], degree_bins[2,i], sep = "-")
}
colnames(f_cor) = n
rownames(f_cor) = feature_names
f_cor_ci_lower = f_cor
f_cor_ci_upper = f_cor
for(j in 1:ncol(degree_bins)){
  n = length(which(degree > degree_bins[1,j] & degree < degree_bins[2,j]))
  delta = 1.96/sqrt(n - 3)
  print(length(which(degree > degree_bins[1,j] & degree < degree_bins[2,j])))
  for(i in 1:length(features_chosen)){
    x = as.numeric(gene_feature_matched[which(degree > degree_bins[1,j] & degree < degree_bins[2,j]), features_chosen[i]])
    y = vScore[which(degree > degree_bins[1,j] & degree < degree_bins[2,j])]
    f_cor[i,j] = cor(x,y, method = 'spearman', use = 'complete.obs' )
    f_cor_ci_lower[i,j] = tanh(atanh(f_cor[i,j]) - delta) #lower bound of 95% CI, see https://stats.stackexchange.com/questions/18887/how-to-calculate-a-confidence-interval-for-spearmans-rank-correlation
    f_cor_ci_upper[i,j] = tanh(atanh(f_cor[i,j]) + delta) #upper bound of 95% CI
  }
}


# function for error bars
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x-lower,y, x+upper, y, angle=0, code=3, length=length, ...)
}


#A bunch of plots, broken up by degree
for(i in 1:ncol(f_cor)){
  o = order(f_cor[,i], decreasing = F)
  par(mar= c(5,5,1,1))
  pdf(file = paste("Working_figure/Figure6/Figure6X1_feature_by_degree", colnames(f_cor)[i], ".pdf", sep = ""))
  barCenter = barplot(f_cor[o, i], horiz=T, beside=F, xlim=c(-0.6, 0.6), 
                      xlab="Spearman correlation with protein interaction dynamicity",
                      axisnames=F, border=NA, cex.axis=0.7, col = apple_colors[i], main = paste("degree = ", colnames(f_cor)[i]))
  arrows( f_cor_ci_lower[o, i], barCenter, f_cor_ci_upper[o, i], barCenter, angle = 90, lwd = 1.5, code = 3, length = 0.05)
  text(-0.6, barCenter, labels= rownames(f_cor)[o], cex = 0.8, xpd = TRUE)
  dev.off()
}

#The main plot using all data
i = ncol(f_cor)
o = order(f_cor[,i], decreasing = F)

pdf(file = "Working_figure/Figure6/Figure6X1_feature_by_degree_main.pdf")
par(mar= c(5,8,1,1))
barCenter = barplot(f_cor[o, i], horiz=T, beside=F, xlim=c(-0.6, 0.3), 
                    xlab="Spearman correlation with protein interaction dynamicity",
                    axisnames=F, border=NA, cex.axis=1, col = apple_colors[1])
arrows( f_cor_ci_lower[o, i], barCenter, f_cor_ci_upper[o, i], barCenter, angle = 90, lwd = 1.5, code = 3, length = 0.05)
text(-0.65, barCenter, labels= rownames(f_cor)[o], cex = 0.8, xpd = TRUE)
dev.off()

#Comparison by degree of some major features
fs = c(17, 15, 3, 2, 1)
f_cor = f_cor[fs, 10:8]
f_cor_ci_lower = f_cor_ci_lower[fs, 10:8]
f_cor_ci_upper = f_cor_ci_upper[fs, 10:8]
pdf(file = "Working_figure/Figure6/Figure6X1_feature_by_degree_main_comparison.pdf")
par(mar= c(5,8,1,1))
barCenter = barplot(t(f_cor), horiz=T, beside=T, xlim=c(-0.6, 0.4), 
                    xlab="Spearman correlation with protein interaction dynamicity",
                    axisnames=F, border=NA, cex.axis=1, col = apple_colors[3:5])
arrows( t(f_cor_ci_lower), barCenter, t(f_cor_ci_upper), barCenter, angle = 90, lwd = 1.5, code = 3, length = 0.05)
text(-0.65, barCenter[2,], labels= rownames(f_cor), cex = 1, xpd = TRUE)
legend(.2, 5, c("0-3", "4-14", "15+"), bty = "n", title = "PPI degree", 
       fill = apple_colors[5:3])
dev.off()

#x = as.data.frame(f_cor)
#x$features = features_chosen




