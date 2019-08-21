require(ggplot2)

#get PPI pairs and dynamicity
setwd("~/Dropbox/PPiSeq_02/")
#load("Working_data/Positive_PPI_environment/variation_score.Rfile")
v = as.matrix(read.csv("Working_data/Positive_PPI_environmentVariation_score_PPI_environment_neg_zero.csv"))
a = sapply(as.character(v[,1]), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
pairs = x
dynamicity = as.numeric(v[,3])
env.number = as.numeric(v[,2])

#Make a yeast systematic name to Entrez table
x = read.csv("Working_data/SC_gene_info.txt", sep = "\t", header = T)
entrez = x[,3]
y = x[,7]
z = sapply(as.character(y), strsplit, ",")
yid = 1:length(z)
for(i in 1:length(z)){
  yid[i] = z[[i]][1]
}
entrez_to_yid = yid
names(entrez_to_yid) = entrez  
yid_to_entrez = entrez
names(yid_to_entrez) = yid

#Make a coexpression matrix
setwd("~/Dropbox/PPiSeq_02/Working_data/CoexpressDB_R")
f = list.files()
m = matrix(NA, length(f), length(f))
colnames(m) = f
rownames(m) = f
for(i in f){
  x = read.delim(i, header = T)
  m[i,as.character(x[,1])] = x[,2]
}
#Rename rows and columns
colnames(m) = entrez_to_yid[colnames(m)]
rownames(m) = entrez_to_yid[rownames(m)]
coex = m
#Save this file
save(coex, file = "~/Dropbox/PPiSeq_02/Working_data/Coexpression/CoExpressDB.Rfile")


#Get CoExpression for all PPIs
x = 1:length(dynamicity)
x[] = NA
for(i in 1:length(x)){
  if(pairs[i,1] %in% colnames(coex) & pairs[i,2] %in% colnames(coex) )x[i] = coex[pairs[i,1], pairs[i,2]]
}
ppi_coex = x
cor(dynamicity, ppi_coex, use = "complete.obs")
plot(dynamicity, ppi_coex)

#make barplot by environment number
df = data.frame(ppi_coex, env.number)
df$env.number = as.factor(df$env.number) 
p <- ggplot(df, aes(x=env.number, y=ppi_coex, fill=env.number)) +   
  geom_boxplot(notch = T) +  labs(x="Number of environments in which PPI is detected", y = "Co-expression Mutual Rank")
p +theme_classic() + scale_fill_brewer(palette = "Blues") + theme(legend.position="none") 


#Get the geometric mean co-expression score for each protein
geomean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  #gets the geometric mean
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}
g = unique(c(pairs[,1], pairs[,2]))
g_coex = 1:length(g); g_coex[] = NA
for(i in 1:length(g)){
  x = pairs[pairs[,1] == g[i] | pairs[,2] == g[i], ]
  if(length(x) == 2){
    x= matrix(x, 1, 2)
  }
  for(j in 1:nrow(x)){
    if(x[j,1] %in% colnames(coex) & x[j,2] %in% colnames(coex)) y[j] = coex[x[j,1], x[j,2]]
  }
  g_coex[i] = geomean(y) #the geomtric mean co-expression value for each protein
}
names(g_coex) = g

#get number of PPIs for each gene
gene_feature = as.matrix(read.table("Working_data/geneFeatures_022415_EK.txt", header = T, sep = "\t")) # 6438, 39
ppi.degree = as.numeric(gene_feature[,4])
names(ppi.degree) = gene_feature[,1]
g_ppi = ppi.degree[g] 

