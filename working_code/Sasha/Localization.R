# This script compares PPIs discoverd by PPiSeq to Localization data from Yolanda et al.
# It finds that PPIs that are found in DMSO are more likely to colocalize when compared to PPIs discovered across all environments
# 65% colocalization vs 29%
# p< 2.2e-16 t-test after 1000 bootstraps
# Need to test: Is this becasue most PPIs not in DMSO are "dynamic" PPIs and these have different properties from "static" PPIs that are overrepresented in DMSO



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
barplot(number.of.env.colcalized.in.DMSO, xlab = "Number of environemnts in which a PPI is observed", 
        ylab = "Percent colocalized in Synthethic Media")

#Find co-localization rate in DMSO for PPIs found in different numbers of environments + found in DMSO
env.number = apply(ppi_by_env, 1, sum, na.rm = T)
number.of.env.colcalized.in.DMSO = 1:ncol(ppi_by_env)
names(number.of.env.colcalized.in.DMSO) = 1:9
for(i in 1:ncol(ppi_by_env)){
  x = ppi[env.number == i ,] #PPIs seen in i env
  y = find_colacalized(x)
  number.of.env.colcalized.in.DMSO[i] = percent_colacalized(y)
}
barplot(number.of.env.colcalized.in.DMSO, xlab = "Number of environemnts in which a PPI is observed", 
        ylab = "Percent colocalized in Synthethic Media")

x = find_colacalized(ppi)


ppi_pos = sum(x == 1, na.rm = T) #PPIs with detected colocalization
ppi_neg = sum(x == 0, na.rm = T) #PPIs with no detected colocalization
ppi_na = sum(is.na(x)) #PPIs without localization data


#Load full set of PPIs across all environments
all = read.csv("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/All_PPI_for_coannotation.csv")
a = sapply(as.character(all[,1]), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
control = x
#Find Overlaps
x = 1:nrow(control)
x[] = NA
for(i in 1:length(x)){
  a = control[i,1]
  b = control[i,2]
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
control_pos = sum(x == 1, na.rm = T) #PPIs with detected colocalization
control_neg = sum(x == 0, na.rm = T) #PPIs with no detected colocalization
control_na = sum(is.na(x))

#Testing for enrichment in co-localization in the benign env for PPIs found in the benign env
control_pos/(control_neg + control_pos) #Fraction of all PPIs that are colocalized in benign env = 0.29
ppi_pos/(ppi_neg + ppi_pos) #Fraction of benign enc PPIs that are colocalized in benign env = 0.65
#conclusion: PPIs found in DMSO are more likely to co-localize in DMSO
#Need to test: This may be becasue most PPIs not in DMSO are "dynamic" PPIs and these have different properties from static PPIs


#Do a bootstrap to see if this difference is significant 
y = 1:1000
for(j in 1:1000){
  x = sample(1:nrow(control), nrow(ppi))
  z = x
  x[] = NA
  for(i in z){
    a = control[i,1]
    b = control[i,2]
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
  y[j] = sum(x == 1, na.rm = T)/(sum(x == 1, na.rm = T) + sum(x == 0, na.rm = T))
}
mean(y) #.2913
sd(y) #.0087
t.test( y, mu = ppi_pos/(ppi_neg + ppi_pos)) # p< 2.2e-16