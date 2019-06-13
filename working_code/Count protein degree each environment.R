# Source some basic functions froma function.R in Github repository
source_https <- function(u, unlink.tmp.certs = FALSE) {
        # load package
        require(RCurl)
        # read script lines from website using a security certificate
        if(!file.exists("cacert.pem")) download.file(url="http://curl.haxx.se/ca/cacert.pem", destfile = "cacert.pem")
        script <- getURL(u, followlocation = TRUE, cainfo = "cacert.pem")
        if(unlink.tmp.certs) unlink("cacert.pem")
        
        # parase lines and evealuate in the global environement
        eval(parse(text = script), envir= .GlobalEnv)
}
source_https("https://raw.githubusercontent.com/sashaflevy/PPiSeq/master/working_code/function.R", unlink.tmp.certs = TRUE)

## This function will generate a degree matrix for PPI network in each environment
setwd("~/Dropbox/PPiSeq_02/Working_data/Positive_PPI_environment/")
PPI_count = as.matrix(read.csv("PPI_environment_count_summary.csv"))
protein = unique(as.vector(split_string_vector(PPI_count[,1]))) # 2051
protein_matrix = matrix(0, length(protein), 10)
protein_matrix[,1] = protein
for(i in 1:9){
        PPI_indiv = PPI_count[which(PPI_count[, 2 + i] == "1"), 1]
        protein_degree = protein_degree_count(PPI_indiv)
        protein_matrix[,1+i] = protein_degree[match(protein, protein_degree[,1]),2]
        
}
protein_matrix[is.na(protein_matrix)] = 0
degree_var = rep(0, length(protein))
for(i in 1:length(degree_var)){
        degree= as.numeric(protein_matrix[i,2:10])
        degree_var[i] = sd(degree)/mean(degree)
}
protein_matrix_final = cbind(protein_matrix, degree_var)
colnames(protein_matrix_final) = c("protein", colnames(PPI_count)[3:11], "CV_degree")
protein_matrix_final_order = protein_matrix_final[order(as.numeric(protein_matrix_final[,ncol(protein_matrix_final)]), decreasing = T),]
csvWriter(protein_matrix_final_order, "Protein_degreee_in_each_environment.csv")
