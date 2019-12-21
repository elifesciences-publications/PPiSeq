# This script compares PPIs discoverd by PPiSeq to Localization data from Yolanda Chong et al. and Gene Ontology
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
x = read_excel_allsheets("~/Dropbox/PPiSeq_02/Paper_data/Outside_datasets/Yolanda_localization/TableS2_LOC_scores.xlsx")

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
save(yolanda_loc, file = "~/Dropbox/PPiSeq_02/Paper_data/Outside_datasets/Yolanda_localization/yolonda_loc.Rfile")


#Get mean WT localization for Yolanda
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
save(yolonda_benign_loc_matrix, file = "~/Dropbox/PPiSeq_02/Paper_data/Outside_datasets/Yolanda_localization/yolonda_benign_loc_matrix.Rfile")

#Get GO localization
x = as.matrix(read.csv("~/Dropbox/PPiSeq_02/Paper_data/Outside_datasets/GO_term_files/go_slim_mapping_tab_20190405.txt", sep = '\t'))
x = x[which(x[,4] == "C"),]
comp = unique(x[,5])
genes = unique(x[,1])
y = matrix(0, length(genes), length(comp))
colnames(y) = comp
rownames(y) = genes
for(i in 1:length(comp)){
  z = unique(x[which(x[,5] == comp[i])])
  y[z, i] = 1
}
GO_loc_matrix = y
save(GO_loc_matrix, file = "~/Dropbox/PPiSeq_02/Paper_data/Outside_datasets/Yolanda_localization/GO_loc_matrix.Rfile")
GO_loc_matrix = GO_loc_matrix[,c(1, 5:24)] #remove general terms  "other", "membrane",  "cellular_component" 

#Import PPIs by environemnt
#load("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/ppi_by_env.Rfile")
PPI_count_filter = read.csv("~/Dropbox/PPiSeq_02/Paper_data/Useful_datasets/PPI_environment_count_summary_SD_merge_filter.csv")
PPI_variation = read.csv("~/Dropbox/PPiSeq_02/Paper_data/Useful_datasets/Variation_score_PPI_environment_neg_zero_SD_merge_filter.csv")
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

# #Find co-localization rate in DMSO for PPIs unique to each env
# env.number = apply(ppi_by_env, 1, sum, na.rm = T)
# env.specific.PPIs.colcalized.in.DMSO = 1:ncol(ppi_by_env)
# names(env.specific.PPIs.colcalized.in.DMSO) = colnames(ppi_by_env)
# for(i in 1:ncol(ppi_by_env)){
#   x = ppi[env.number == 1 & ppi_by_env[,i] == 1,] #PPIs unique to env i
#   y = find_colacalized(x, yolonda_benign_loc_matrix)
#   env.specific.PPIs.colcalized.in.DMSO[i] = percent_colacalized(y)
# }
# barplot(env.specific.PPIs.colcalized.in.DMSO, xlab = "Environemnt in which a PPI is uniquely observed", 
#         ylab = "Percent colocalized in Synthethic Media")

#Find co-localization rate in DMSO for PPIs found in different numbers of environments
env.number = apply(ppi_by_env, 1, sum, na.rm = T)
number.of.env.colcalized.in.DMSO = 1:ncol(ppi_by_env)
names(number.of.env.colcalized.in.DMSO) = 1:9
for(i in 1:ncol(ppi_by_env)){
  x = ppi[env.number == i ,] #PPIs seen in i env
  y = find_colacalized(x, GO_loc_matrix)
  number.of.env.colcalized.in.DMSO[i] = percent_colacalized(y)
}
barplot(number.of.env.colcalized.in.DMSO, xlab = "Number of environments in which a PPI is observed", 
        ylab = "Percent colocalized in Synthethic Media")

# #Find co-localization rate in DMSO for PPIs found in different numbers of environments + found in DMSO
# env.number = apply(ppi_by_env, 1, sum, na.rm = T)
# number.of.env.colcalized.in.DMSO = 1:ncol(ppi_by_env)
# names(number.of.env.colcalized.in.DMSO) = 1:9
# for(i in 1:ncol(ppi_by_env)){
#   x = ppi[env.number == i & ppi_by_env[,"DMSO"] == 1,] #PPIs seen in i env
#   y = find_colacalized(x, yolonda_benign_loc_matrix)
#   number.of.env.colcalized.in.DMSO[i] = percent_colacalized(y)
# }
# barplot(number.of.env.colcalized.in.DMSO, xlab = "Number of environemnts in which a PPI is observed, found in DMSO", 
#         ylab = "Percent colocalized in Synthethic Media")


#Do a bootstrap on Yolanda to get error bars -- sample PPIs with replacement
number.of.env.colcalized.bootstrap = matrix(NA, 1000, ncol(ppi_by_env))
colnames(number.of.env.colcalized.bootstrap) = 1:9
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
colnames(number.of.env.colcalized.bootstrap) = 1:9
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


pdf(file = "~/Dropbox/PPiSeq_02/Working_figure/Figure4_date_party_hubs/Figure4B_colocalization_rate.pdf", height = 7, width = 4 )
plot(1:9, m, xlab = "Environments in which a PPI is observed", 
     ylab = "Percent colocalized", type = 'l', ylim = c(30,100), pch = 16, cex = 1.7, lwd =2, col = "black" ) 
points(1:9, m,  ylim = c(30,90), pch = 16, cex = 1.7, col = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027"))
points(1:9, mgo,type = 'l', lwd = 2, lty = "dashed")
points(1:9, mgo,   pch = 16, cex = 1.7, 
      col = c("#4575b4","#74add1","#abd9e9","#e0f3f8","#ffffbf","#fee090", "#fdae61","#f46d43","#d73027"))
arrows( 1:9, u, 1:9,l,  lwd = 1.5, code = 3, length = 0.05, angle = 90)
arrows( 1:9, ugo, 1:9,lgo,  lwd = 1.5, code = 3, length = 0.05, angle = 90)
legend(3, 40, legend = c("Gene Ontology", "Fluoresence"),lty=c("dashed", "solid") )
dev.off()



