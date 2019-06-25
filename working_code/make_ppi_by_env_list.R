load(file = "~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/ppi_by_env.Rfile")
a = sapply(rownames(ppi_by_env), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
ppi = x #detected PPIs across all environments

#make 0/1 PPI matrix in each environment
x = unique(ppi[,1])
y = unique(ppi[,2])
#reorder by compartment
load(file="~/Dropbox/PPiSeq_02/Working_data/GOSlim_CC.Rfile" ) #cc, compartment for each ORF
xc = cc[x]
names(xc) = x
x = names(sort(xc))
yc = cc[y]
names(yc) = y
y = names(sort(yc))
#make matrix
z = matrix(0, length(x), length(y))
rownames(z) = x
colnames(z) = y
l = list(z, z, z, z, z, z, z, z, z)
names(l) = colnames(ppi_by_env)
for(i in 1:nrow(ppi)){
  if(ppi[i,1] %in% x & ppi[i,2] %in% y){ 
    l[[1]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,1]
    l[[2]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,2]
    l[[3]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,3]
    l[[4]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,4]
    l[[5]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,5]
    l[[6]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,6]
    l[[7]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,7]
    l[[8]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,8]
    l[[9]][ppi[i,1], ppi[i,2]] = ppi_by_env[i,9]
  }
}
ppi_by_env_list = l #list of matrices of whether or not each PPI is detected
save(ppi_by_env_list, file = "~/Dropbox/PPiSeq_02/Working_data/ppi_by_env_list.Rfile")