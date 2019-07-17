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
#Commonly used colors
apple_colors = c("#5AC8FA", "#FFCC00", "#FF9500", "#FF2D55", "#007AFF", "#4CD964", "#FF3B30",
                 "#8E8E93", "#EFEFF4", "#CECED2", "#000000", "007AFF")

#(2) Get all the possible orders of 9 environments, and check the PPI number changes
setwd("~/Dropbox/PPiSeq_02/")
count_summary = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv") # 12333
#count_summary= count_summary[which(count_summary[,2] != "1"),] # 5392
### Followsing this order DMSO, H2O2, HU, Dox, Forskolin, Raffinose, NaCl, 16C, FK506
row_num = c(factorial(9), factorial(8), factorial(7), factorial(6), factorial(5), factorial(4), factorial(3), factorial(2),factorial(1))

matrix_count_1 = matrix(0, factorial(9), 9)
matrix_count_2 = matrix(0, factorial(9), 9)
matrix_count_3 = matrix(0, factorial(9), 9)
matrix_count_4 = matrix(0, factorial(9), 9)
matrix_count_5 = matrix(0, factorial(9), 9)
matrix_count_6 = matrix(0, factorial(9), 9)
matrix_count_7 = matrix(0, factorial(9), 9)
matrix_count_8 = matrix(0, factorial(9), 9)
matrix_count_9 = matrix(0, factorial(9), 9)

matrix_constructer = rep(0, nrow(count_summary))
all_env = 3:11
# First column
bulk_1 = row_num[2]
index_1 = bulk_1 -1 
for(i in 1:9){
        column_chosen = as.numeric(count_summary[,i +2])
        matrix_chosen = data.frame(column_chosen, matrix_constructer)
        sum_row = rowSums(matrix_chosen)
        environ_1_1 = length(which(sum_row >= 1))
        environ_1_2 = length(which(sum_row >= 2))
        environ_1_3 = length(which(sum_row >= 3))
        environ_1_4 = length(which(sum_row >= 4))
        environ_1_5 = length(which(sum_row >= 5))
        environ_1_6 = length(which(sum_row >= 6))
        environ_1_7 = length(which(sum_row >= 7))
        environ_1_8 = length(which(sum_row >= 8))
        environ_1_9 = length(which(sum_row >= 9))
        index= (bulk_1 * i - index_1):(bulk_1*i)
        matrix_count_1[index,1] = environ_1_1
        matrix_count_2[index,1] = environ_1_2
        matrix_count_3[index,1] = environ_1_3
        matrix_count_4[index,1] = environ_1_4
        matrix_count_5[index,1] = environ_1_5
        matrix_count_6[index,1] = environ_1_6
        matrix_count_7[index,1] = environ_1_7
        matrix_count_8[index,1] = environ_1_8
        matrix_count_9[index,1] = environ_1_9
}

# Second column
bulk_2 = row_num[3]
index_2 = bulk_2-1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h + 2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                column_chosen_2 = as.numeric(count_summary[,env_remaining_1[i]])
                matrix_chosen = data.frame(column_chosen_1, column_chosen_2)
                sum_row = rowSums(matrix_chosen)
                environ_2_1 = length(which(sum_row >= 1))
                environ_2_2 = length(which(sum_row >= 2))
                environ_2_3 = length(which(sum_row >= 3))
                environ_2_4 = length(which(sum_row >= 4))
                environ_2_5 = length(which(sum_row >= 5))
                environ_2_6 = length(which(sum_row >= 6))
                environ_2_7 = length(which(sum_row >= 7))
                environ_2_8 = length(which(sum_row >= 8))
                environ_2_9 = length(which(sum_row >= 9))
                index = (bulk_2 * i - index_2):(bulk_2*i) + env1_index
                matrix_count_1[index,2] = environ_2_1
                matrix_count_2[index,2] = environ_2_2
                matrix_count_3[index,2] = environ_2_3
                matrix_count_4[index,2] = environ_2_4
                matrix_count_5[index,2] = environ_2_5
                matrix_count_6[index,2] = environ_2_6
                matrix_count_7[index,2] = environ_2_7
                matrix_count_8[index,2] = environ_2_8
                matrix_count_9[index,2] = environ_2_9
        }
}


# Third column
bulk_3 = row_num[4] #720
index_3 = bulk_3 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                        column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                        column_chosen_3 = as.numeric(count_summary[,env_remaining_2[j]])
                        matrix_chosen = data.frame(column_chosen_1, column_chosen_2, column_chosen_3)
                        sum_row = rowSums(matrix_chosen)
                        environ_3_1 = length(which(sum_row >= 1))
                        environ_3_2 = length(which(sum_row >= 2))
                        environ_3_3 = length(which(sum_row >= 3))
                        environ_3_4 = length(which(sum_row >= 4))
                        environ_3_5 = length(which(sum_row >= 5))
                        environ_3_6 = length(which(sum_row >= 6))
                        environ_3_7 = length(which(sum_row >= 7))
                        environ_3_8 = length(which(sum_row >= 8))
                        environ_3_9 = length(which(sum_row >= 9))
                        index = (bulk_3 * j - index_3):(bulk_3*j) + env2_index
                        matrix_count_1[index, 3] = environ_3_1 
                        matrix_count_2[index, 3] = environ_3_2
                        matrix_count_3[index, 3] = environ_3_3
                        matrix_count_4[index, 3] = environ_3_4
                        matrix_count_5[index, 3] = environ_3_5
                        matrix_count_6[index, 3] = environ_3_6
                        matrix_count_7[index, 3] = environ_3_7
                        matrix_count_8[index, 3] = environ_3_8
                        matrix_count_9[index, 3] = environ_3_9
                }
                
        } 
}

sum(matrix_count_1[,3] > matrix_count_1[,2])

matrix_count_1[720,]
matrix_count_1[721,]
matrix_count_1[40320,]
matrix_count_1[40320 -720,]
matrix_count_1[40320 -720 + 1,]
length(unique(matrix_count_1[,3])) # 80

# Fourth column
bulk_4 = row_num[5] # 120
index_4 = bulk_4 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                column_chosen_4 = as.numeric(count_summary[,env_remaining_3[k]])
                                matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                           column_chosen_3, column_chosen_4)
                                sum_row = rowSums(matrix_chosen)
                                environ_4_1 = length(which(sum_row >= 1))
                                environ_4_2 = length(which(sum_row >= 2))
                                environ_4_3 = length(which(sum_row >= 3))
                                environ_4_4 = length(which(sum_row >= 4))
                                environ_4_5 = length(which(sum_row >= 5))
                                environ_4_6 = length(which(sum_row >= 6))
                                environ_4_7 = length(which(sum_row >= 7))
                                environ_4_8 = length(which(sum_row >= 8))
                                environ_4_9 = length(which(sum_row >= 9))
                                index = (bulk_4 * k - index_4):(bulk_4*k) + env3_index
                                matrix_count_1[index, 4] = environ_4_1 
                                matrix_count_2[index, 4] = environ_4_2 
                                matrix_count_3[index, 4] = environ_4_3 
                                matrix_count_4[index, 4] = environ_4_4 
                                matrix_count_5[index, 4] = environ_4_5 
                                matrix_count_6[index, 4] = environ_4_6 
                                matrix_count_7[index, 4] = environ_4_7 
                                matrix_count_8[index, 4] = environ_4_8 
                                matrix_count_9[index, 4] = environ_4_9 
                        }
                          
                }
                
        } 
}


matrix_count_1[362880 -120,]
matrix_count_1[362880-120 + 1,]
length(unique(matrix_count_1[,4])) # 121
sum(matrix_count_1[,4] > matrix_count_1[,3])

# Fifth column
bulk_5 = row_num[6] # 24
index_5 = bulk_5 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                env4_index = (k -1) * bulk_4 + env3_index
                                chosen_env_4 = env_remaining_3[k]
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                for(l in 1:5){
                                        column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                        column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                        column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                        column_chosen_4 = as.numeric(count_summary[,chosen_env_4])
                                        column_chosen_5 = as.numeric(count_summary[,env_remaining_4[l]])
                                        matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                                   column_chosen_3, column_chosen_4,
                                                                   column_chosen_5)
                                        sum_row = rowSums(matrix_chosen)
                                        environ_5_1 = length(which(sum_row >= 1))
                                        environ_5_2 = length(which(sum_row >= 2))
                                        environ_5_3 = length(which(sum_row >= 3))
                                        environ_5_4 = length(which(sum_row >= 4))
                                        environ_5_5 = length(which(sum_row >= 5))
                                        environ_5_6 = length(which(sum_row >= 6))
                                        environ_5_7 = length(which(sum_row >= 7))
                                        environ_5_8 = length(which(sum_row >= 8))
                                        environ_5_9 = length(which(sum_row >= 9))
                                        index = (bulk_5 * l - index_5):(bulk_5*l) + env4_index
                                        matrix_count_1[index, 5] = environ_5_1
                                        matrix_count_2[index, 5] = environ_5_2
                                        matrix_count_3[index, 5] = environ_5_3
                                        matrix_count_4[index, 5] = environ_5_4
                                        matrix_count_5[index, 5] = environ_5_5
                                        matrix_count_6[index, 5] = environ_5_6
                                        matrix_count_7[index, 5] = environ_5_7
                                        matrix_count_8[index, 5] = environ_5_8
                                        matrix_count_9[index, 5] = environ_5_9
                                }
                                    
                        }
                        
                }
                
        } 
}
matrix_count_1[362880-24,]
matrix_count_1[362880 - 24 +1,]
matrix_count_1[362880,]
sum(matrix_count_1[,5] > matrix_count_1[,4])


# Sixth column
bulk_6 = row_num[7] # 6
index_6 = bulk_6 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                env4_index = (k -1) * bulk_4 + env3_index
                                chosen_env_4 = env_remaining_3[k]
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                for(l in 1:5){
                                        env5_index = (l -1) * bulk_5 + env4_index
                                        chosen_env_5 = env_remaining_4[l]
                                        env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                        for(m in 1:4){
                                                column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                                column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                                column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                                column_chosen_4 = as.numeric(count_summary[,chosen_env_4])
                                                column_chosen_5 = as.numeric(count_summary[,chosen_env_5])
                                                column_chosen_6 = as.numeric(count_summary[,env_remaining_5[m]])
                                                matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                                           column_chosen_3, column_chosen_4,
                                                                           column_chosen_5, column_chosen_6)
                                                sum_row = rowSums(matrix_chosen)
                                                environ_6_1 = length(which(sum_row >= 1))
                                                environ_6_2 = length(which(sum_row >= 2))
                                                environ_6_3 = length(which(sum_row >= 3))
                                                environ_6_4 = length(which(sum_row >= 4))
                                                environ_6_5 = length(which(sum_row >= 5))
                                                environ_6_6 = length(which(sum_row >= 6))
                                                environ_6_7 = length(which(sum_row >= 7))
                                                environ_6_8 = length(which(sum_row >= 8))
                                                environ_6_9 = length(which(sum_row >= 9))
                                                index = (bulk_6 * m - index_6):(bulk_6*m) + env5_index
                                                matrix_count_1[index, 6] = environ_6_1 
                                                matrix_count_2[index, 6] = environ_6_2 
                                                matrix_count_3[index, 6] = environ_6_3 
                                                matrix_count_4[index, 6] = environ_6_4 
                                                matrix_count_5[index, 6] = environ_6_5 
                                                matrix_count_6[index, 6] = environ_6_6 
                                                matrix_count_7[index, 6] = environ_6_7 
                                                matrix_count_8[index, 6] = environ_6_8 
                                                matrix_count_9[index, 6] = environ_6_9 
                                        }
                                            
                                }
                                
                        }
                        
                }
                
        } 
}
matrix_count_1[1:7,]
sum(matrix_count_1[,6] > matrix_count_1[,5])


# Seventh column
bulk_7 = row_num[8] # 2
index_7 = bulk_7 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                env4_index = (k -1) * bulk_4 + env3_index
                                chosen_env_4 = env_remaining_3[k]
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                for(l in 1:5){
                                        env5_index = (l -1) * bulk_5 + env4_index
                                        chosen_env_5 = env_remaining_4[l]
                                        env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                        for(m in 1:4){
                                                env6_index = (m -1) * bulk_6 + env5_index
                                                chosen_env_6 = env_remaining_5[m]
                                                env_remaining_6 = env_remaining_5[which(env_remaining_5 != chosen_env_6)]
                                                for(n in 1:3){
                                                        column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                                        column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                                        column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                                        column_chosen_4 = as.numeric(count_summary[,chosen_env_4])
                                                        column_chosen_5 = as.numeric(count_summary[,chosen_env_5])
                                                        column_chosen_6 = as.numeric(count_summary[,chosen_env_6])
                                                        column_chosen_7 = as.numeric(count_summary[,env_remaining_6[n]])
                                                        matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                                                   column_chosen_3, column_chosen_4,
                                                                                   column_chosen_5, column_chosen_6,
                                                                                   column_chosen_7)
                                                        sum_row = rowSums(matrix_chosen)
                                                        environ_7_1 = length(which(sum_row >= 1))
                                                        environ_7_2 = length(which(sum_row >= 2))
                                                        environ_7_3 = length(which(sum_row >= 3))
                                                        environ_7_4 = length(which(sum_row >= 4))
                                                        environ_7_5 = length(which(sum_row >= 5))
                                                        environ_7_6 = length(which(sum_row >= 6))
                                                        environ_7_7 = length(which(sum_row >= 7))
                                                        environ_7_8 = length(which(sum_row >= 8))
                                                        environ_7_9 = length(which(sum_row >= 9))
                                                        index = (bulk_7 * n - index_7):(bulk_7*n) + env6_index
                                                        matrix_count_1[index, 7] = environ_7_1 
                                                        matrix_count_2[index, 7] = environ_7_2  
                                                        matrix_count_3[index, 7] = environ_7_3  
                                                        matrix_count_4[index, 7] = environ_7_4  
                                                        matrix_count_5[index, 7] = environ_7_5  
                                                        matrix_count_6[index, 7] = environ_7_6 
                                                        matrix_count_7[index, 7] = environ_7_7  
                                                        matrix_count_8[index, 7] = environ_7_8  
                                                        matrix_count_9[index, 7] = environ_7_9  
                                                }
                                                 
                                        }
                                        
                                }
                                
                        }
                        
                }
                
        } 
}
matrix_count_1[1:7,]
matrix_count_1[40300:40320,]
sum(matrix_count_1[,7] > matrix_count_1[,6])

# Eighth column
bulk_8 = row_num[9] # 1
index_8 = bulk_8 -1
for(h in 1:9){
        env1_index = (h-1) * bulk_1
        chosen_env_1 = h +2
        env_remaining_1 = all_env[all_env != chosen_env_1]
        for(i in 1:8){
                env2_index = (i-1) * bulk_2 + env1_index
                chosen_env_2 = env_remaining_1[i]
                env_remaining_2 = env_remaining_1[which(env_remaining_1 != chosen_env_2)]
                for(j in 1:7){
                        env3_index = (j -1) * bulk_3 + env2_index
                        chosen_env_3 = env_remaining_2[j]
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        for(k in 1:6){
                                env4_index = (k -1) * bulk_4 + env3_index
                                chosen_env_4 = env_remaining_3[k]
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                for(l in 1:5){
                                        env5_index = (l -1) * bulk_5 + env4_index
                                        chosen_env_5 = env_remaining_4[l]
                                        env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                        for(m in 1:4){
                                                env6_index = (m -1) * bulk_6 + env5_index
                                                chosen_env_6 = env_remaining_5[m]
                                                env_remaining_6 = env_remaining_5[which(env_remaining_5 != chosen_env_6)]
                                                for(n in 1:3){
                                                        env7_index = (n -1) * bulk_7 + env6_index
                                                        chosen_env_7 = env_remaining_6[n]
                                                        env_remaining_7 = env_remaining_6[which(env_remaining_6 != chosen_env_7)]
                                                        for(o in 1:2){
                                                                column_chosen_1 = as.numeric(count_summary[,chosen_env_1])
                                                                column_chosen_2 = as.numeric(count_summary[,chosen_env_2])
                                                                column_chosen_3 = as.numeric(count_summary[,chosen_env_3])
                                                                column_chosen_4 = as.numeric(count_summary[,chosen_env_4])
                                                                column_chosen_5 = as.numeric(count_summary[,chosen_env_5])
                                                                column_chosen_6 = as.numeric(count_summary[,chosen_env_6])
                                                                column_chosen_7 = as.numeric(count_summary[,chosen_env_7])
                                                                column_chosen_8 = as.numeric(count_summary[,env_remaining_7[o]])
                                                                matrix_chosen = data.frame(column_chosen_1, column_chosen_2, 
                                                                                           column_chosen_3, column_chosen_4,
                                                                                           column_chosen_5, column_chosen_6,
                                                                                           column_chosen_7, column_chosen_8)
                                                                sum_row = rowSums(matrix_chosen)
                                                                environ_8_1 = length(which(sum_row >= 1))
                                                                environ_8_2 = length(which(sum_row >= 2))
                                                                environ_8_3 = length(which(sum_row >= 3))
                                                                environ_8_4 = length(which(sum_row >= 4))
                                                                environ_8_5 = length(which(sum_row >= 5))
                                                                environ_8_6 = length(which(sum_row >= 6))
                                                                environ_8_7 = length(which(sum_row >= 7))
                                                                environ_8_8 = length(which(sum_row >= 8))
                                                                environ_8_9 = length(which(sum_row >= 9))
                                                                index = (bulk_8 * o - index_8):(bulk_8*o) + env7_index
                                                                matrix_count_1[index, 8] = environ_8_1
                                                                matrix_count_2[index, 8] = environ_8_2 
                                                                matrix_count_3[index, 8] = environ_8_3 
                                                                matrix_count_4[index, 8] = environ_8_4 
                                                                matrix_count_5[index, 8] = environ_8_5 
                                                                matrix_count_6[index, 8] = environ_8_6 
                                                                matrix_count_7[index, 8] = environ_8_7 
                                                                matrix_count_8[index, 8] = environ_8_8 
                                                                matrix_count_9[index, 8] = environ_8_9 
                                                        }
                                                        
                                                }
                                                
                                        }
                                        
                                }
                                
                        }
                        
                }
                
        } 
}
matrix_count[1:10,]
matrix_count[40300:40320,]
sum(matrix_count[,8] > matrix_count[,7])
length(unique(matrix_count[,8]))
# Ninth column
matrix_count_1[,9] = nrow(count_summary)
matrix_count_2[,9] = length(which(as.numeric(count_summary[,2]) >= 2))
matrix_count_3[,9] = length(which(as.numeric(count_summary[,2]) >= 3))
matrix_count_4[,9] = length(which(as.numeric(count_summary[,2]) >= 4))
matrix_count_5[,9] = length(which(as.numeric(count_summary[,2]) >= 5))
matrix_count_6[,9] = length(which(as.numeric(count_summary[,2]) >= 6))
matrix_count_7[,9] = length(which(as.numeric(count_summary[,2]) >= 7))
matrix_count_8[,9] = length(which(as.numeric(count_summary[,2]) >= 8))
matrix_count_9[,9] = length(which(as.numeric(count_summary[,2]) >= 9))

colnames(matrix_count_1) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
colnames(matrix_count_2) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
colnames(matrix_count_3) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
colnames(matrix_count_4) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
colnames(matrix_count_5) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
colnames(matrix_count_6) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
colnames(matrix_count_7) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
colnames(matrix_count_8) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
colnames(matrix_count_9) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")

csvWriter(matrix_count_1, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_1.csv")
csvWriter(matrix_count_2, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_2.csv")
csvWriter(matrix_count_3, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_3.csv")
csvWriter(matrix_count_4, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_4.csv")
csvWriter(matrix_count_5, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_5.csv")
csvWriter(matrix_count_6, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_6.csv")
csvWriter(matrix_count_7, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_7.csv")
csvWriter(matrix_count_8, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_8.csv")
csvWriter(matrix_count_9, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_9.csv")


setwd("~/Dropbox/PPiSeq_02/")
matrix_count = csvReader_T("Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order.csv")
matrix_count_high = csvReader_T("Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order_larger_two.csv")

matrix_count_vec = as.numeric(as.vector(t(matrix_count)))
num_env = rep(1:9, nrow(matrix_count))

matrix_count_high_vec = as.numeric(as.vector(t(matrix_count_high[,2:9])))
num_env_high = rep(2:9, nrow(matrix_count_high))
# linear regression for all data
relation = lm(matrix_count_vec ~ num_env)
linear_count = relation$coefficients[2]* (1:9) + relation$coefficients[1]

"""
# Bionomial regression
polynomial_relation = lm(matrix_count_vec ~ poly(num_env, 2, raw=TRUE))
polynomial_count = polynomial_relation$coefficients[1] + polynomial_relation$coefficients[2] * (1:9) +
        polynomial_relation$coefficients[3] * ((1:9)^2)
"""
# Trinomial regression for high-quality data
polynomial_relation_3 = lm(matrix_count_high_vec ~ poly(num_env_high, 3, raw=TRUE))
polynomial_count_3 = polynomial_relation_3$coefficients[1] + polynomial_relation_3$coefficients[2] * (2:9) +
        polynomial_relation_3$coefficients[3] * ((2:9)^2) +  polynomial_relation_3$coefficients[4] * ((2:9)^3)
# Mean count
mean_count = colMeans(matrix_count)
mean_count_high = colMeans(matrix_count_high)

pdf("Working_figure/Figure2/Figure2B_PPI_number_accumulation_environments_combine.pdf", width= 5, height =5)
par(mar = c(4,4,2,4))
plot(1:9, as.numeric(matrix_count[1,]), xlim = c(1,9), ylim = c(4000,14000), type = "l",
     col = rgb(0.5,0.5,0.5,alpha = 0.05),lwd = 0.05,
     xlab = "Number of environments assayed",
     ylab = "Number of PPIs", main = NA, bty = "n")
axis(2, at= seq(4000, 14000, by = 2000), labels = seq(4000, 14000, by = 2000))
axis(1, at= 1:9, labels = 1:9)
for(i in 2:nrow(matrix_count)){
        lines(1:9, as.numeric(matrix_count[i,]), col = rgb(0.5,0.5,0.5, alpha = 0.05), lwd = 0.05)
}

par(new = T)

par(mar = c(4,4,2,4))
plot(2:9, as.numeric(matrix_count_high[1,2:9]), axes =F, xlim = c(1,9), ylim = c(4000, 6000),
     xlab = NA, ylab = NA, type = "l",
     col = rgb(0.5,0.5,0.5,alpha = 1),lwd = 0.05)
axis(4, at= seq(4000, 6000, by = 250), labels = seq(4000, 6000, by = 250), 
     col = "darkgreen", col.axis = "darkgreen")

mtext(side = 4,  ">= 2 envrionments", line = 2.5, col = "darkgreen")
for(i in 1:nrow(matrix_count_high)){
        lines(2:9, as.numeric(matrix_count_high[i,2:9]), col = rgb(0,0.5,0.5, alpha = 0.05), lwd = 0.05)
}

dev.off()

#### Separate two figures, and make a small one for high-quality data
pdf("Working_figure/Figure2/Figure2B_PPI_number_accumulation_environments_large_two_small.pdf", width= 3, height =3)
par(mar = c(2,1,1,4))
plot(2:9, as.numeric(matrix_count_high[1,2:9]), axes =F, xlim = c(1,9), ylim = c(4000, 6000),
     xlab = NA, ylab = NA, type = "l",
     col = rgb(0.5,0.5,0.5,alpha = 0.05),lwd = 0.05)
axis(4, at= seq(4000, 6000, by = 500), labels = seq(4000, 6000, by = 500), 
     cex.axis = 0.6, mgp = c(3,0.5,0))
     #col = "darkgreen", col.axis = "darkgreen")
axis(1, at= 1:9, labels = 1:9, cex.axis = 0.6, mgp = c(3, 0.5, 0))
mtext(side = 4,  "Number of PPIs", line = 2, cex = 0.6)
for(i in 1:nrow(matrix_count_high)){
        lines(2:9, as.numeric(matrix_count_high[i,2:9]), col = rgb(0.5,0.5,0.5, alpha = 0.05), lwd = 0.05)
}
lines(2:9, mean_count_high[2:9], col = apple_colors[7])
dev.off()





pdf("Working_figure/Figure2/Figure2B_PPI_number_accumulation_environments_large_two.pdf", width= 5, height =5)
par(mar = c(4,4,2,1))
plot(2:9, as.numeric(matrix_count_high[1,2:9]), xlim = c(1,9), ylim = c(4000,6000), type = "l",
     col = rgb(0.5,0.5,0.5,alpha = 0.05),lwd = 0.05,
     xlab = "Number of environments assayed",
     ylab = "Number of PPIs", main = NA, bty = "n")
axis(2, at= seq(4000, 6000, by = 250), labels = seq(4000, 6000, by = 250))
axis(1, at= 1:9, labels = 1:9)
for(i in 2:nrow(matrix_count_high)){
        lines(2:9, as.numeric(matrix_count_high[i,2:9]), col = rgb(0.5,0.5,0.5, alpha = 0.05), lwd = 0.05)
}

lines(2:9, mean_count_high[2:9], col = apple_colors[7])
#lines(c(1,9), as.numeric(matrix_count[1,c(1, 9)]), col = "orange", lty = 1)
#lines(1:9, linear_count, col = apple_colors[5])
#lines(1:9, polynomial_count_3, col = apple_colors[5])
lines(2:9, polynomial_count_3, col = apple_colors[5])
#legend("topleft", c("Mean count", "Linear regression", "Bionomial regression", "Trinomial regression", "Straight line"),
#lty = c(1,1,2,2,1), col = c("red", "green", "cyan", "blue", "orange"), bty= "n")
legend("topleft", c("Mean count", "Trinomial regression"),lty = c(1,1), col = apple_colors[c(7,5)], bty= "n")

dev.off()

######################################################
# Remove PPIs that have been detected only in one environment and then make the same figure
setwd("~/Dropbox/PPiSeq_02/")
count_summary = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv") # 1233
count_summary[1,]
count_summary = count_summary[which(count_summary[,2] != "1"),] # 5392
# Then run the above code and generate a new matrix

