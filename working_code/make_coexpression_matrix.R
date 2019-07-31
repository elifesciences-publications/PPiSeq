#Make coexpression matrix from mutual rank scores 


#get PPI pairs and dynamicity
setwd("~/Dropbox/PPiSeq_02/")
load("Working_data/Positive_PPI_environment/variation_score.Rfile")
v = variation_score
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
