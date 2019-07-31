# This script compares PPIs discoverd by PPiSeq to Localization data from Yolanda Chong et al.
# It finds that PPIs that are found in  more environments are more likely to colocalize 


#import localization data
require(readxl)
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.matrix)
  names(x) <- sheets
  x
}
x = read_excel_allsheets("~/Dropbox/PPiSeq_02/Working_data/Yolanda_localization/TableS2_LOC_scores.xlsx")

#convert to a list of matrices
y = x
for(i in 1:length(x)){
  a = x[[i]]
  b = data.matrix(a[4:nrow(a), 3:ncol(a)]) 
  class(b) <- "numeric"
  colnames(b) = as.matrix(a[3,3:ncol(a)])
  rownames(b) = as.matrix(a[4:nrow(a),1])
  y[[i]] = b
}

#Save the list
yolanda_loc = y
save(yolanda_loc, file = "~/Dropbox/PPiSeq_02/Working_data/Yolanda_localization/yolonda_loc.Rfile")


#Get mean WT localization
a = unique(c(rownames(y[[1]]), rownames(y[[2]]), rownames(y[[3]])))
wt = matrix(NA, length(a), ncol(y[[1]]))
colnames(wt) = colnames(y[[1]])
rownames(wt) = a
for(i in 1:length(a)){
  x = matrix(NA, 3, ncol(y[[1]]))
  if(a[i] %in% rownames(y[[1]]))  x[1,] = y[[1]][a[i],]
  if(a[i] %in% rownames(y[[2]]))  x[2,] = y[[2]][a[i],]
  if(a[i] %in% rownames(y[[3]]))  x[3,] = y[[3]][a[i],]
  wt[a[i], ] = apply(x, 2, mean, na.rm = T)
}
yolonda_benign_loc_matrix = wt
save(yolonda_benign_loc_matrix, file = "~/Dropbox/PPiSeq_02/Working_data/Yolanda_localization/yolonda_benign_loc_matrix.Rfile")

#Import PPIs by environemnt
load("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/ppi_by_env.Rfile")
a = sapply(as.character(rownames(ppi_by_env)), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
ppi = x #detected PPIs across all environments


#Find Colacalization b/w PPI functions
find_colacalized = function(ppi){
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
  print(sum(x == 1, na.rm = T))
  print(sum(x == 0, na.rm = T))
  sum(x == 1, na.rm = T)/ (sum(x == 1, na.rm = T) + sum(x == 0, na.rm = T))
}

#Find co-localization rate in DMSO for PPIs unique to each env
env.number = apply(ppi_by_env, 1, sum, na.rm = T)
env.specific.PPIs.colcalized.in.DMSO = 1:ncol(ppi_by_env)
names(env.specific.PPIs.colcalized.in.DMSO) = colnames(ppi_by_env)
for(i in 1:ncol(ppi_by_env)){
  x = ppi[env.number == 1 & ppi_by_env[,i] == 1,] #PPIs unique to env i
  y = find_colacalized(x)
  env.specific.PPIs.colcalized.in.DMSO[i] = percent_colacalized(y)
}
barplot(env.specific.PPIs.colcalized.in.DMSO, xlab = "Environemnt in which a PPI is uniquely observed", 
        ylab = "Percent colocalized in Synthethic Media")

#Find co-localization rate in DMSO for PPIs found in different numbers of environments
env.number = apply(ppi_by_env, 1, sum, na.rm = T)
number.of.env.colcalized.in.DMSO = 1:ncol(ppi_by_env)
names(number.of.env.colcalized.in.DMSO) = 1:9
for(i in 1:ncol(ppi_by_env)){
  x = ppi[env.number == i ,] #PPIs seen in i env
  y = find_colacalized(x)
  number.of.env.colcalized.in.DMSO[i] = percent_colacalized(y)
}
barplot(number.of.env.colcalized.in.DMSO, xlab = "Number of environments in which a PPI is observed", 
        ylab = "Percent colocalized in Synthethic Media")

#Find co-localization rate in DMSO for PPIs found in different numbers of environments + found in DMSO
env.number = apply(ppi_by_env, 1, sum, na.rm = T)
number.of.env.colcalized.in.DMSO = 1:ncol(ppi_by_env)
names(number.of.env.colcalized.in.DMSO) = 1:9
for(i in 1:ncol(ppi_by_env)){
  x = ppi[env.number == i & ppi_by_env[,"DMSO"] == 1,] #PPIs seen in i env
  y = find_colacalized(x)
  number.of.env.colcalized.in.DMSO[i] = percent_colacalized(y)
}
barplot(number.of.env.colcalized.in.DMSO, xlab = "Number of environemnts in which a PPI is observed, found in DMSO", 
        ylab = "Percent colocalized in Synthethic Media")


#Do a bootstrap to get error bars -- sample PPIs with replacement
number.of.env.colcalized.bootstrap = matrix(NA, 1000, ncol(ppi_by_env))
colnames(number.of.env.colcalized.bootstrap) = 1:9
for(j in 1:1000){
  b = sample(1:nrow(ppi_by_env), nrow(ppi_by_env), replace = T) #the new bootstrapped indexes
  ppi_by_env_b = ppi_by_env[b,] #the new bootstrapped ppis
  ppi_b = ppi[b,]
  env.number_b = apply(ppi_by_env_b, 1, sum, na.rm = T)
  for(i in 1:ncol(ppi_by_env_b)){
    x = ppi_b[env.number_b == i ,] #PPIs seen in i env
    y = find_colacalized(x)
    number.of.env.colcalized.bootstrap[j,i] = percent_colacalized(y)
  }
}
m = apply(number.of.env.colcalized.bootstrap, 2, mean, na.rm = T)*100
sd = apply(number.of.env.colcalized.bootstrap, 2, sd, na.rm = T)*100
sem = sd/sqrt(1000)
u = m+sd
l = m-sd

pdf(file = "~/Dropbox/PPiSeq_02/Working_figure/Figure3/Figure3B_colocalization_rate.pdf", height = 6, width = 4)
plot(1:9, m, xlab = "Environments in which a PPI is observed", 
        ylab = "Percent colocalized in synthethic media", type = 'b', ylim = c(30,90), pch = 16)
arrows( 1:9, u, 1:9,l,  lwd = 1.5, code = 3, length = 0.05, angle = 90)
dev.off()
