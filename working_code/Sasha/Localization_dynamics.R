#load localization matrices of HU and benign from yolanda et al.
load("~/Dropbox/PPiSeq_02/Working_data/Yolanda_localization/yolonda_benign_loc_matrix.Rfile")
load("~/Dropbox/PPiSeq_02/Working_data/Yolanda_localization/yolonda_hu_loc_matrix.Rfile")

#Load PPIs from HU and Benign environments
bp = read.csv("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/DMSO_Pos_PPI_real.csv")
hup = read.csv("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Hydroxyurea_Pos_PPI_real.csv")

#Find those unique to hu benign and shared
b = as.character(bp[,1])
hu = as.character(hup[,1])
hu_specific  = hu[!hu %in% b]
b_specific = b[!b %in% hu]
shared = hu[hu %in% b]
a = sapply(hu_specific, strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
hus_ppi = x
a = sapply(shared, strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
shared_ppi = x
a = sapply(b_specific, strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
bs_ppi = x

#Find colacalization in both environments
coloc = function(ppi, loc){
  x = 1:nrow(ppi)
  x[] = NA
  for(i in 1:length(x)){
    a = ppi[i,1]
    b = ppi[i,2]
    if(a %in% rownames(loc) & b %in% rownames(loc)){
      t = loc[c(a,b),]
      t = t > 0
      u = apply(t, 2, sum)
      if(sum(u == 2) > 0){
        x[i] = 1
      }else{
        x[i] = 0
      }
    }
  }
  y = c(sum(x == 1, na.rm = T), sum(x == 0, na.rm = T), sum(is.na(x)), 
        sum(x == 1, na.rm = T)/(sum(x == 0, na.rm = T) + sum(x == 1, na.rm = T)))
  names(y) = c("Colocalized", "Not Colocalized", "No Information", "Percent Colocalized")
  return(y)
}
coloc(shared_ppi, yolonda_benign_loc_matrix)
coloc(shared_ppi, yolonda_hu_loc_matrix)
coloc(hus_ppi, yolonda_benign_loc_matrix)
coloc(hus_ppi, yolonda_hu_loc_matrix)
coloc(bs_ppi, yolonda_benign_loc_matrix)
coloc(bs_ppi, yolonda_hu_loc_matrix)
#conclusion: no enrichment for environment-specific colocalization amoung environment-specific PPIs

#Find identities of environment-specific colocalization and PPI
coloc_id = function(ppi, loc){
  x = matrix(NA, nrow(ppi), 2)
  for(i in 1:nrow(x)){
    a = ppi[i,1]
    b = ppi[i,2]
    if(a %in% rownames(loc) & b %in% rownames(loc)){
      t = loc[c(a,b),]
      t = t > 0
      u = apply(t, 2, sum)
      if(sum(u == 2) > 0){
        x[i,] = c(a,b)
      }
    }
  }
  return(x)
}
x = coloc_id(hus_ppi, yolonda_benign_loc_matrix)
y = coloc_id(hus_ppi, yolonda_hu_loc_matrix)
z = y[which(is.na(x[,1]) & !is.na(y[,1])),] #hu sprecific PPI and localization

## Look at dynamics -- work in progress
a = unique(as.vector(z))
b = 1:length(a)
names(b) = a
b[] = NA
for(i in 1:length(a)){
  b[i] = sum(z == a[i])
}
y = yolonda_hu_loc_matrix[names(b),]
y = y>0
apply(y, 2, sum) # a lot in the cell perephery
z = yolonda_benign_loc_matrix[which(rownames(yolonda_benign_loc_matrix) %in% names(b)),]
z = z>0
apply(z, 2, sum)

periph = rownames(y[y[,4],])
a = z[rownames(z)[rownames(z) %in% periph], ]
b = y[rownames(z)[rownames(z) %in% periph], ]

b-a
b+a
