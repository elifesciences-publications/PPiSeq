#load localization data from yolanda et al.
load("~/Dropbox/PPiSeq_02/Working_data/Yolanda_localization/yolonda_loc.Rfile")
y = yolanda_loc

#Get mean HU localization
a = unique(c(rownames(y[[13]]), rownames(y[[14]]), rownames(y[[15]])))
hu = matrix(NA, length(a), ncol(y[[13]]))
colnames(hu) = colnames(y[[13]])
rownames(hu) = a
for(i in 1:length(a)){
  x = matrix(NA, 3, ncol(y[[13]]))
  if(a[i] %in% rownames(y[[13]]))  x[1,] = y[[13]][a[i],]
  if(a[i] %in% rownames(y[[14]]))  x[2,] = y[[14]][a[i],]
  if(a[i] %in% rownames(y[[15]]))  x[3,] = y[[15]][a[i],]
  hu[a[i], ] = apply(x, 2, mean, na.rm = T)
}
yolonda_hu_loc_matrix = hu
save(yolonda_hu_loc_matrix, file = "~/Dropbox/PPiSeq_02/Working_data/Yolanda_localization/yolonda_hu_loc_matrix.Rfile")


#Import PPIs in HU
p = read.csv("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/Hydroxyurea_Pos_PPI_real.csv")
a = sapply(as.character(p[,1]), strsplit, "_")
x = matrix(NA, length(a), 2)
for(i in 1:length(a)){
  x[i,] = a[[i]]
}
huppi = x

#Find Overlaps
x = 1:nrow(huppi)
x[] = NA
for(i in 1:length(x)){
  a = huppi[i,1]
  b = huppi[i,2]
  if(a %in% rownames(hu) & b %in% rownames(hu)){
    t = hu[c(a,b),]
    t = t > 0
    u = apply(t, 2, sum)
    if(sum(u == 2) > 0){
      x[i] = 1
    }else{
      x[i] = 0
    }
  }
}
huppi_pos = sum(x == 1, na.rm = T)
huppi_neg = sum(x == 0, na.rm = T)
huppi_na = sum(is.na(x))


#Load full set of PPIs
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
  if(a %in% rownames(hu) & b %in% rownames(hu)){
    t = hu[c(a,b),]
    t = t > 0
    u = apply(t, 2, sum)
    if(sum(u == 2) > 0){
      x[i] = 1
    }else{
      x[i] = 0
    }
  }
}
control_pos = sum(x == 1, na.rm = T)
control_neg = sum(x == 0, na.rm = T)
control_na = sum(is.na(x))

#Enriched colacaliztion
control_pos/(control_neg + control_pos)
huppi_pos/(huppi_neg + huppi_pos)

#Do a bootstrap
y = 1:1000
for(j in 1:1000){
  x = sample(1:nrow(control), nrow(huppi))
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
mean(y) #.317
sd(y) #.0085
t.test( y, mu = huppi_pos/(huppi_neg + huppi_pos)) #p<10e-16
