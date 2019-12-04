#### Caculate an variation score for each PPI based on their fitness values across
#### different environments. Use the normalized value and Keep the primary normalized fitess even a PPI is called negative
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

setwd("~/Dropbox/PPiSeq_02/")
PPI_norm = csvReader_T("Working_data/Positive_PPI_environment/Normalized_fitness_PPI_all_primary.csv")
PPI_count = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv")
PPI_norm_matrix = matrix(0, nrow(PPI_norm), 10)
PPI_norm_matrix[,1] = PPI_norm[,1]
for(i in 1:nrow(PPI_norm)){
        for(j in 3:11){
                mean_fit = mean(as.numeric(na.omit(PPI_norm[i,c(j, j + 10)])))
                if(is.na(mean_fit)| mean_fit < 0){
                        PPI_norm_matrix[i, j-1] = 0
                }
                else{
                        PPI_norm_matrix[i,j-1] = mean_fit
                }
                
        }  
}
variation_score = rep(0, nrow(PPI_norm_matrix))
for(i in 1:length(variation_score)){
        fitness = as.numeric(PPI_norm_matrix[i,2:10])
        variation_score[i] = sd(fitness)/mean(fitness)
}
environment_number = as.numeric(PPI_count[match(PPI_norm_matrix[,1], PPI_count[,1]),2])
# Include the variation score into the matrix, and output the matrix
PPI_norm_matrix_final = cbind(PPI_norm_matrix[,1],environment_number, variation_score,
                              PPI_norm_matrix[,2:10])
colnames(PPI_norm_matrix_final) = c("PPI", "Environment_number", "Variation_score",
                                    "SD", "H2O2", "HU", "Dox", "Forskolin", "Raffinose",
                                    "NaCl", "16C", "FK506")
save(PPI_norm_matrix_final, file = "Working_data/Positive_PPI_environment/variation_score.Rfile") #for quick retrieval