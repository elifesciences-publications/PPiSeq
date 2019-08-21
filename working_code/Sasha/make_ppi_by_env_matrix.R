#Make a matrix of PPI vs. env

setwd("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/")
source("~/PPiseq/HU/R_code/function.R")

#Import data
x = list()
x[[1]] = read.csv("DMSO_Pos_PPI_real.csv")
x[[2]] = read.csv("Forskolin_Pos_PPI_real.csv")
x[[3]] = read.csv("FK506_Pos_PPI_real.csv")
x[[4]] = read.csv("Raffinose_Pos_PPI_real.csv")
x[[5]] = read.csv("NaCl_Pos_PPI_real.csv")
x[[6]] = read.csv("H2O2_Pos_PPI_real.csv")
x[[7]] = read.csv("Doxorubicin_Pos_PPI_real.csv")
x[[8]] = read.csv("Cold_16C_Pos_PPI_real.csv")
x[[9]] = read.csv("Hydroxyurea_Pos_PPI_real.csv")

#Find all PPIs
y = NULL
for(i in 1:length(x)){
  y = c(y, as.character(x[[i]][,1]))
}
y = unique(y) 

#Make a matrix 
z = matrix(0, length(y), 9)
colnames(z) = c("DMSO", "Forskolin", "FK506", "Raffinose", "NaCl", "H2O2", "Doxorubicin", "16C", "Hydroxyurea")
rownames(z) = y
for(i in 1:length(x)){
  z[y %in% x[[i]][,1], i] = 1
}

ppi_by_env = z
save(ppi_by_env, file = "ppi_by_env.Rfile")