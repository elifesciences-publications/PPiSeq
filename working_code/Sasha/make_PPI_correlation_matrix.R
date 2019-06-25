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
#reorder by compartment
load(file="~/Dropbox/PPiSeq_02/Working_data/GOSlim_CC.Rfile" ) #cc, compartment for each ORF
xc = cc[x]
names(xc) = x
x = names(sort(xc))
yc = cc[y]
names(yc) = y
y = names(sort(yc))
#make matrix
z = matrix(NA, length(x), length(y))
rownames(z) = x
colnames(z) = y
l = list(z, z, z, z, z, z, z, z, z)
names(l) = colnames(fit)
for(i in 1:nrow(ppi)){
  if(ppi[i,1] %in% x & ppi[i,2] %in% y){ 
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
}
fitness_by_env_list = l
save(fitness_by_env_list, file = "~/Dropbox/PPiSeq_02/Working_data/fitness_by_env_list.Rfile")
#calculate breaks for plotting
u = unique(sort(xc))
bb = 1:(length(u) + 1) #bait bins
bb[1] = 0
for(i in 1:length(u)){
  bb[i + 1] = max(which(sort(xc) == u[i]))
}
bcc = u
u = unique(sort(yc))
pb = 1:(length(u) + 1) #prey bins
pb[1] = 0
for(i in 1:length(u)){
  pb[i + 1] = max(which(sort(yc) == u[i]))
}
pcc = u
cc_bins = list(bait_cellular_component = bcc, bait_bins = bb, prey_cellular_component = pcc, prey_bins = pb)
save(cc_bins, file = "~/Dropbox/PPiSeq_02/Working_data/cc_bins.Rfile")

#make a big matrix for bait
x = NULL
for(i in 1:length(fitness_by_env_list)){
  x = cbind(x, fitness_by_env_list[[i]])
}
n = NULL
for(i in 1:length(fitness_by_env_list)){
  n = c(n, paste(names(fitness_by_env_list)[i], colnames(fitness_by_env_list[[i]]), sep = '_'))
}
colnames(x) = n
rownames(x) = rownames(fitness_by_env_list[[1]])
bait_fitness_all_env = x
save(bait_fitness_all_env , file = "~/Dropbox/PPiSeq_02/Working_data/bait_fitness_all_env.Rfile")

#make a big matrix for prey
x = NULL
for(i in 1:length(fitness_by_env_list)){
  x = rbind(x, fitness_by_env_list[[i]])
}
n = NULL
for(i in 1:length(fitness_by_env_list)){
  n = c(n, paste(names(fitness_by_env_list)[i], rownames(fitness_by_env_list[[i]]), sep = '_'))
}
rownames(x) = n
colnames(x) = colnames(fitness_by_env_list[[1]])
prey_fitness_all_env = x
save(prey_fitness_all_env , file = "~/Dropbox/PPiSeq_02/Working_data/prey_fitness_all_env.Rfile")

#make bait correlation matrix
x = matrix(NA, nrow(bait_fitness_all_env), nrow(bait_fitness_all_env))
y = rownames(bait_fitness_all_env)
rownames(x) = y
colnames(x) = y
for(i in 1:nrow(bait_fitness_all_env)){
  for(j in 1:nrow(bait_fitness_all_env)){
    x[i,j] = cor(bait_fitness_all_env[y[i],], bait_fitness_all_env[y[j],], use = "na.or.complete")
  }
}
#remove empty rows, columns
r = names(which(apply(x, 1, sum, na.rm = T) == 0))
x = x[which(!rownames(x) %in% r), which(!colnames(x) %in% r)]
#save
bait_cor_matrix = x
save(bait_cor_matrix, file = "~/Dropbox/PPiSeq_02/Working_data/bait_cor_matrix.Rfile")

#make prey correlation matrix
load("~/Dropbox/PPiSeq_02/Working_data/prey_fitness_all_env.Rfile")
x = matrix(NA, ncol(prey_fitness_all_env), ncol(prey_fitness_all_env))
y = colnames(prey_fitness_all_env)
rownames(x) = y
colnames(x) = y
for(i in 1:ncol(prey_fitness_all_env)){
  for(j in 1:ncol(prey_fitness_all_env)){
    x[i,j] = cor(prey_fitness_all_env[,y[i]], prey_fitness_all_env[,y[j]], use = "na.or.complete")
  }
}
r = names(which(apply(x, 1, sum, na.rm = T) == 0))
x = x[which(!rownames(x) %in% r), which(!colnames(x) %in% r)]
prey_cor_matrix = x
save(prey_cor_matrix, file = "~/Dropbox/PPiSeq_02/Working_data/prey_cor_matrix.Rfile")


#visualize
load("~/Dropbox/PPiSeq_02/Working_data/bait_cor_matrix.Rfile")
require(ggplot2)
require(gplots)
x = bait_cor_matrix
y = as.data.frame(x)
p = 

pdf(file = "~/Desktop/test.pdf" )
heatmap.2(x, na.color = "black", dendrogram = "col", labRow = F, labCol = F, trace = "none")
dev.off()


