### Remake Figure S6 after removing data of 16C
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R") # Load commonly used fuctions
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

### barplot to show how many of them have been reproted
all_PPI_matrix_final = csvReader_T("Datasets_generated_by_preprocessing/PPI_environment_count_summary_SD_merge_filter.csv")
all_PPI_matrix_final[,2] = as.numeric(all_PPI_matrix_final[,2]) - as.numeric(all_PPI_matrix_final[,10])
nrow(all_PPI_matrix_final) # 12981
all_PPI_matrix_final = all_PPI_matrix_final[which(all_PPI_matrix_final[,2] > 0),]
nrow(all_PPI_matrix_final) # 9887
all_PPI_matrix_final_remove = all_PPI_matrix_final[,c(1:9, 11)]
csvWriter(all_PPI_matrix_final_remove, 'SFigure9_related_data/PPI_environment_count_summary_SD_merge_filter_removing_cold.csv')

all_PPI_matrix_final = csvReader_T('SFigure9_related_data/PPI_environment_count_summary_SD_merge_filter_removing_cold.csv')
reported_PPI = csvReader_T("Outsourced_datasets/BIOGRID//multiple_validated_PPI.csv")
PCA_lower = as.matrix(read.table("Outsourced_datasets/Tarassov_PPI_PPV_80.txt", header= T, sep = "\t"))
PCA_lower_PPI = paste(PCA_lower[,1], PCA_lower[,4], sep = "_")
PCA_lower_PPI_matrix = cbind(PCA_lower_PPI, rep(1, length(PCA_lower_PPI)))
PCA_lower_PPI_reported = match_both_direction(PCA_lower_PPI_matrix, reported_PPI[,1])
PCA_lower_PPI_unreported = PCA_lower_PPI_matrix[which(!PCA_lower_PPI_matrix[,1] %in% PCA_lower_PPI_reported[,1]),]

matrix_PPI_env_rep = matrix(0, 3, 8)
for(i in 1:8){
  all = all_PPI_matrix_final[which(as.numeric(all_PPI_matrix_final[,2]) == i),]
  all_reported_PCA_low = match_both_direction(all, PCA_lower_PPI_unreported[,1])
  all_reported_BioGrid = match_both_direction(all,reported_PPI[,1])
  all_reported = c(all_reported_PCA_low[,1], all_reported_BioGrid[,1])
  all_unreported = all[which(!all[,1] %in% all_reported),]
  matrix_PPI_env_rep[1,i] = nrow(all_unreported)
  matrix_PPI_env_rep[2,i] = nrow(all_reported_PCA_low)
  matrix_PPI_env_rep[3,i] = nrow(all_reported_BioGrid)
}
matrix_PPI_env_rep[3,] # 248 139 131 159 155 189 252 474
all_PPI_count = matrix_PPI_env_rep[1,] + matrix_PPI_env_rep[2,] + matrix_PPI_env_rep[3,]
all_PPI_count # 4882 1111  696  578  615  584  605  816
ratio_BioGrid = matrix_PPI_env_rep[3,]/all_PPI_count
ratio_BioGrid # 0.05079885 0.12511251 0.18821839 0.27508651 0.25203252 0.32363014 0.41652893 0.58088235
ratio_BioGrid_reported = c("5.1%", "12.5%", "18.8%", "27.5%", "25.2%", "32.3%", "41.7%", "58%")
ratio_PCA_low = matrix_PPI_env_rep[2,]/all_PPI_count
ratio_PCA_low # 0.05428103 0.13681368 0.22270115 0.25951557 0.33170732 0.35787671 0.43966942 0.32720588
ratio_PCA_low_overlapped = c("5.4%", "13.7%", "22.3%", "26.0%", "33.2%", "35.8%", "44.0%", "32.7%")

library(RColorBrewer)
col_chosen = apple_colors[c(1,3,4)]
pdf("Figures/SFigures/SFigure9/FigureS9A_Number_environments_PPI_reproted.pdf", height = 5, width = 5)
par(mar = c(3,5,2,1))
barCenter = barplot(matrix_PPI_env_rep, horiz=F, beside=F, ylim=c(0,6000), ylab="Number of PPIs",
                    space= c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6),
                    col= col_chosen, axisnames=F, border=NA)
legend("topright", legend=c("BioGRID", "mDHFR-PCA (80% < PPV < 98.2%)", "Previously unreported"), 
       fill=col_chosen[c(3,2,1)], cex = 0.8, bty="n", border=FALSE)
text(x= barCenter, y = all_PPI_count + 100, labels = ratio_PCA_low_overlapped, 
     cex=0.7, xpd = TRUE, col= col_chosen[2]) 
text(x= barCenter, y = all_PPI_count + 250, labels = ratio_BioGrid_reported, 
     cex=0.7, xpd = TRUE, col= col_chosen[3]) 
text(x= barCenter, y = -150, labels = as.character(1:8), xpd = TRUE)
text(median(barCenter), y = -350, labels = "Number of environments in which a PPI is identified", xpd = TRUE)
dev.off()

########### Recalculate the PPI mutability score after removing 16 C data (the output file was in the PPiSeq_additional_data)
## Can be neglected if only making the figure
setwd('/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/')
PPI_norm = csvReader_T("Normalized_multiple_files/Normalized_fitness_PPI_all_primary_SD_merge.csv")
PPI_norm = PPI_norm[,c(1:9,11:19, 21)]
PPI_count = csvReader_T("Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge.csv")
PPI_count[,2] = as.numeric(PPI_count[,2]) - as.numeric(PPI_count[,10])
environment_number = as.numeric(PPI_count[match(PPI_norm_matrix[,1], PPI_count[,1]),2])

PPI_norm_matrix = matrix(0, nrow(PPI_norm), 9)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
  for(j in 3:10){
    a = as.numeric(PPI_norm[i, j])
    b = as.numeric(PPI_norm[i, j + 9])
    if ((!is.na(a)) & a < 0){ ### Consider negative values to be zero
      a = 0
    }
    if ((!is.na(b)) & b < 0){
      b = 0
    }
    mean_fit = mean(as.numeric(na.omit(c(a,b))))
    if(is.na(mean_fit)){
      PPI_norm_matrix[i, j-1] = 0
    }
    else{
      PPI_norm_matrix[i,j-1] = mean_fit
    }
    
  }  
}
variation_score = rep(0, nrow(PPI_norm_matrix))
for(i in 1:length(variation_score)){
  fitness = as.numeric(PPI_norm_matrix[i,2:9])
  variation_score[i] = sd(fitness)/mean(fitness)
}



# Include the variation score into the matrix, and output the matrix
PPI_norm_matrix_final = cbind(PPI_norm_matrix[,1],environment_number, variation_score,
                              PPI_norm_matrix[,2:9])
colnames(PPI_norm_matrix_final) = c("PPI", "Environment_number", "Variation_score",
                                    "SD", "H2O2", "HU", "Dox", "Forskolin", "Raffinose",
                                    "NaCl", "FK506") # 13764
PPI_count_filter = csvReader_T("/Volumes/zmliu_02/PPiseq_03/Combine_environments_SD_merge/Positive_PPI_remove_promiscuous/PPI_environment_count_summary_SD_merge_filter.csv") 
PPI_norm_matrix_final_filter = PPI_norm_matrix_final[which(PPI_norm_matrix_final[,1] %in% PPI_count_filter[,1]),] # 12981
PPI_norm_matrix_final_filter = PPI_norm_matrix_final_filter[which(PPI_norm_matrix_final_filter[,2] != "0"),] # 9887
cor(as.numeric(PPI_norm_matrix_final_filter[,2]), as.numeric(PPI_norm_matrix_final_filter[,3]), method="spearman") # -0.7469293

csvWriter(PPI_norm_matrix_final_filter, "~/Desktop/PPiSeq_additional_data/SFigure9_related_data/Variation_score_PPI_environment_neg_zero_SD_merge_filter_removing_cold.csv")

### Make second figure
setwd("~/Desktop/PPiSeq_additional_data/")
#A plot of co-expression mutual rank by environment number
load("Outsourced_datasets/Coexpression/CoExpressDB.Rfile") # load coexpresion data
#load data
ppi_coex = coex
v = csvReader_T("SFigure9_related_data/Variation_score_PPI_environment_neg_zero_SD_merge_filter_removing_cold.csv")
a = sapply(as.character(v[,1]), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
pairs = x
dynamicity = as.numeric(v[,3])
env.number = as.numeric(v[,2])

#Get CoExpression for all PPIs
x = 1:length(dynamicity)
x[] = NA
for(i in 1:length(x)){
  if(pairs[i,1] %in% colnames(coex) & pairs[i,2] %in% colnames(coex) )x[i] = coex[pairs[i,1], pairs[i,2]]
}
ppi_coex = x

#make plot of co-expression by environment number
df = data.frame(ppi_coex, env.number)
df$env.number = as.factor(df$env.number)
library(ggplot2)
pdf(file = "Figures/SFigures/Sfigure9/FigureS9B_coexpression_by_env_number.pdf", height = 6, width = 4)
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Environments in which a PPI is observed", y = "Co-expression mutual rank")
p + theme_classic() + theme(axis.text.y = element_text(size=10, angle=45), 
                            axis.text.x = element_text(size=10), 
                            panel.background = element_rect(fill = "white"), legend.position="none") + scale_fill_manual(values = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")) 
dev.off()


#### Colocalization 
setwd("~/Desktop/PPiSeq_additional_data/")
load("Outsourced_datasets/Colocalization/yolonda_benign_loc_matrix.Rfile") # Load GFP colocalizaiton data
load("Outsourced_datasets/Colocalization/GO_loc_matrix.Rfile") # Load GO colocalization data
GO_loc_matrix = GO_loc_matrix[,c(1, 5:24)] #remove general terms  "other", "membrane",  "cellular_component" 

#Import PPIs by environemnt
PPI_count_filter = read.csv("SFigure9_related_data/PPI_environment_count_summary_SD_merge_filter_removing_cold.csv")
PPI_variation = read.csv("SFigure9_related_data//Variation_score_PPI_environment_neg_zero_SD_merge_filter_removing_cold.csv")
ppi_by_env = PPI_count_filter[,3:ncol(PPI_count_filter)]
rownames(ppi_by_env) <- PPI_count_filter[,1]

a = sapply(as.character(rownames(ppi_by_env)), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
ppi = x #detected PPIs across all environments


#Find Colacalization b/w PPI functions
find_colacalized = function(ppi, yolonda_benign_loc_matrix){
  wt = yolonda_benign_loc_matrix
  x = 1:nrow(ppi)
  x[] = NA
  for(i in 1:length(x)){
    a = ppi[i,1]
    b = ppi[i,2]
    if(a %in% rownames(wt) & b %in% rownames(wt)){
      t = wt[c(a,b),]
      t = t > 0
      u = apply(t, 2, sum)
      if(sum(u == 2) > 0){
        x[i] = 1
      }else{
        x[i] = 0
      }
    }
  }
  return(x)
}

percent_colacalized = function(x){
  #print(sum(x == 1, na.rm = T))
  #print(sum(x == 0, na.rm = T))
  sum(x == 1, na.rm = T)/ (sum(x == 1, na.rm = T) + sum(x == 0, na.rm = T))
}

env.number = apply(ppi_by_env, 1, sum, na.rm = T)
number.of.env.colcalized.in.DMSO = 1:ncol(ppi_by_env)
names(number.of.env.colcalized.in.DMSO) = 1:8
for(i in 1:ncol(ppi_by_env)){
  x = ppi[env.number == i ,] #PPIs seen in i env
  y = find_colacalized(x, GO_loc_matrix)
  number.of.env.colcalized.in.DMSO[i] = percent_colacalized(y)
}
number.of.env.colcalized.bootstrap = matrix(NA, 1000, ncol(ppi_by_env))
colnames(number.of.env.colcalized.bootstrap) = 1:8

for(j in 1:1000){
  b = sample(1:nrow(ppi_by_env), nrow(ppi_by_env), replace = T) #the new bootstrapped indexes
  ppi_by_env_b = ppi_by_env[b,] #the new bootstrapped ppis
  ppi_b = ppi[b,]
  env.number_b = apply(ppi_by_env_b, 1, sum, na.rm = T)
  for(i in 1:ncol(ppi_by_env_b)){
    x = ppi_b[env.number_b == i ,] #PPIs seen in i env
    y = find_colacalized(x, yolonda_benign_loc_matrix)
    number.of.env.colcalized.bootstrap[j,i] = percent_colacalized(y)
  }
}
m = apply(number.of.env.colcalized.bootstrap, 2, mean, na.rm = T)*100
sd = apply(number.of.env.colcalized.bootstrap, 2, sd, na.rm = T)*100
sem = sd/sqrt(1000)
u = m+sd
l = m-sd

#Do a bootstrap with GO to get error bars -- sample PPIs with replacement
number.of.env.colcalized.bootstrap = matrix(NA, 1000, ncol(ppi_by_env))
colnames(number.of.env.colcalized.bootstrap) = 1:8
for(j in 1:1000){
  b = sample(1:nrow(ppi_by_env), nrow(ppi_by_env), replace = T) #the new bootstrapped indexes
  ppi_by_env_b = ppi_by_env[b,] #the new bootstrapped ppis
  ppi_b = ppi[b,]
  env.number_b = apply(ppi_by_env_b, 1, sum, na.rm = T)
  for(i in 1:ncol(ppi_by_env_b)){
    x = ppi_b[env.number_b == i ,] #PPIs seen in i env
    y = find_colacalized(x, GO_loc_matrix)
    number.of.env.colcalized.bootstrap[j,i] = percent_colacalized(y)
  }
}
mgo = apply(number.of.env.colcalized.bootstrap, 2, mean, na.rm = T)*100
sdgo = apply(number.of.env.colcalized.bootstrap, 2, sd, na.rm = T)*100
semgo = sd/sqrt(1000)
ugo = mgo+sdgo
lgo = mgo-sdgo

pdf(file = "Figures/SFigures/SFigure9/FigureS9C_colocalization_rate.pdf", height = 7, width = 4)
plot(1:8, m, xlab = "Environments in which a PPI is observed", 
     ylab = "Percent colocalized", type = 'l', ylim = c(30,100), pch = 16, cex = 1.7, lwd =2, col = "black" ) 
points(1:8, m,  ylim = c(30,90), pch = 16, cex = 1.7, col = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027"))
points(1:8, mgo,type = 'l', lwd = 2, lty = "dashed")
points(1:8, mgo,   pch = 16, cex = 1.7, 
       col = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027"))
arrows( 1:8, u, 1:8,l,  lwd = 1.5, code = 3, length = 0.05, angle = 90)
arrows( 1:8, ugo, 1:8,lgo,  lwd = 1.5, code = 3, length = 0.05, angle = 90)
legend(3, 40, legend = c("Gene Ontology", "Fluoresence"),lty=c("dashed", "solid"), bty= "n" )
dev.off()

########################################## Correlation with other gene features
setwd("~/Desktop/PPiSeq_additional_data/")
source("function.R")
paper.colors = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027")
#Make a combined data matrix with variation score
variation_score = read.csv("SFigure9_related_data//Variation_score_PPI_environment_neg_zero_SD_merge_filter_removing_cold.csv") # 5928
PP_pair = split_string_vector(variation_score[,1])
protein_unique = unique(as.vector(PP_pair)) 
vScore_protein = cbind(protein_unique, rep(0, length(protein_unique)))
for(i in 1:length(protein_unique)){
  index_protein = unique(c(which(PP_pair[,1] == protein_unique[i]),
                           which(PP_pair[,2] == protein_unique[i])))
  vScore_protein[i,2] = mean(as.numeric(variation_score[index_protein,3]))
}

load("Outsourced_datasets/BIOGRID/GIcount.Rfile")
gene_feature = as.matrix(read.table("Outsourced_datasets/geneFeatures_022415_EK.txt", header = T, sep = "\t")) 
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

pdf(file = "Figures/SFigures/SFigure9/FigureS9D_all_feature_correlation.pdf")
par(mar= c(5,5,1,1))
barCenter = barplot(b.mean[o, i], horiz=T, beside=F, xlim=c(-0.5, 0.5), 
                    xlab="Spearman correlation with PPI mutability",
                    axisnames=F, border=NA, cex.axis=1, col = paper.colors[2])
arrows( b.mean[o, i] - b.sd[o,i], barCenter,  b.mean[o, i] + b.sd[o,i], barCenter, angle = 90, lwd = 1.5, code = 3, length = 0.05)
text(0.3, barCenter, labels= rownames(f_cor)[o], cex = 0.8, xpd = TRUE)
dev.off()

#Comparison by degree of some major features
fs = c(17, 15, 18, 8,  3,  2, 1)
fc = b.mean[fs, 3:1]
fc_ci_lower = b.mean[fs, 3:1]- b.sd[fs, 3:1] 
fc_ci_upper = b.mean[fs, 3:1]+ b.sd[fs, 3:1] 
pdf(file = "Figures/SFigures/SFigure9/FigureS9E_gene_features_PPI_stability.pdf")
par(mar= c(5,5,1,2))
barCenter = barplot(t(fc), horiz=T, beside=T, xlim=c(-0.8, 0.6), 
                    xlab="Correlation with the mutability score",
                    axisnames=F, border=NA, cex.axis=1.2, col = paper.colors[c(2,6,8)])
arrows( t(fc_ci_lower), barCenter, t(fc_ci_upper), barCenter, angle = 90, lwd = 1.5, code = 3, length = 0.05)
text(-0.7, barCenter[2,], labels= rownames(fc), cex = 1.2, xpd = TRUE)
legend(.15, 8, c("0-4", "5-14", "15+"), bty = "n", title = "PPI degree", 
       fill = paper.colors[c(8,6,2)], cex = 1.2)
dev.off()



