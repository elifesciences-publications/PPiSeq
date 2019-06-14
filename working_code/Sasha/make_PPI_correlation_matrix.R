a = read.csv("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/All_PPI_environments_normalized_fit.csv")

#convert to a matrix
b = data.matrix(a[, 2:ncol(a)]) 
class(b) <- "numeric"
colnames(b) = colnames(a)[2:ncol(a)]
rownames(b) = as.matrix(a[,1])
fit = b

#make PPI matrix
a = sapply(rownames(fit), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
ppi = x #detected PPIs across all environments

#make fitness matrix in each environment
x = unique(ppi[,1])
y = unique(ppi[,2])
z = matrix(NA, length(x), length(y))
rownames(z) = x
colnames(z) = y
l = list(z, z, z, z, z, z, z, z, z)
names(l) = colnames(fit)
for(i in 1:nrow(ppi)){
    l[[1]][ppi[i,1], ppi[i,2]] = fit[i,1]
    l[[2]][ppi[i,1], ppi[i,2]] = fit[i,2]
    l[[3]][ppi[i,1], ppi[i,2]] = fit[i,3]
    l[[4]][ppi[i,1], ppi[i,2]] = fit[i,4]
    l[[5]][ppi[i,1], ppi[i,2]] = fit[i,5]
    l[[6]][ppi[i,1], ppi[i,2]] = fit[i,6]
    l[[7]][ppi[i,1], ppi[i,2]] = fit[i,7]
    l[[8]][ppi[i,1], ppi[i,2]] = fit[i,8]
    l[[9]][ppi[i,1], ppi[i,2]] = fit[i,9]
}
fitness_by_env_list = l
save(fitness_by_env_list, file = "~/Dropbox/PPiSeq_02/Working_data/fitness_by_env_list.Rfile")

