#(2) Get all the possible orders of 9 environments, and check the PPI number changes
count_summary = csvReader_T("Working_data/Positive_PPI_environment/PPI_environment_count_summary.csv") # 1233
### Followsing this order DMSO, H2O2, HU, Dox, Forskolin, Raffinose, NaCl, 16C, FK506
row_num = c(factorial(8), factorial(7), factorial(6), factorial(5), factorial(4), factorial(3), factorial(2),factorial(1))
matrix_count = matrix(0, factorial(8), 9)
matrix_count[,1] = length(which(count_summary[,3] == "1"))

all_env = 3:11
# Second column
bulk_2 = row_num[2]
index_2 = bulk_2-1
for(i in 1:8){
        environ_2 = length(which(count_summary[,3] == "1" | count_summary[,i +3] == "1"))
        matrix_count[(bulk_2 * i - index_2):(bulk_2*i),2] = environ_2
}

# Third column
bulk_3 = row_num[3] #720
index_3 = bulk_3 -1
for(i in 1:8){
        env2_index = (i-1) * bulk_2
        for(j in 1:7){
                chosen_env = c(3, i + 3)
                env_remaining = all_env[which(!all_env %in% chosen_env)]
                environ_3 = length(which(count_summary[,3] == "1" | count_summary[,i +3] == "1" 
                                         | count_summary[,env_remaining[j]] == "1"))
                matrix_count[((bulk_3 * j - index_3):(bulk_3*j) + env2_index), 3] = environ_3    
        }
        
}
sum(matrix_count[,3] > matrix_count[,2])

matrix_count[720,]
matrix_count[721,]
matrix_count[40320,]
matrix_count[40320 -720,]
matrix_count[40320 -720 + 1,]
length(unique(matrix_count[,3]))

# Fourth column
bulk_4 = row_num[4] # 120
index_4 = bulk_4 -1
for(i in 1:8){
        env2_index = (i-1) * bulk_2
        chosen_env_2 = c(3, i + 3)
        env_remaining_2 = all_env[which(!all_env %in% chosen_env_2)]
        for(j in 1:7){
                env3_index = env2_index + (j-1) * bulk_3
                chosen_env_3 = env_remaining_2[j]
                for(k in 1:6){
                        env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                        environ_4 = length(which(count_summary[,3] == "1" | 
                                                         count_summary[,i +3] == "1" | 
                                                         count_summary[, chosen_env_3] == "1"|
                                                         count_summary[, env_remaining_3[k]] == "1"))
                        matrix_count[((bulk_4 * k - index_4):(bulk_4*k) + env3_index), 4] = environ_4   
                }
                
        }
        
}
matrix_count[120,]
matrix_count[240,]
matrix_count[40320,]
matrix_count[40320 -120,]
matrix_count[40320 -120 + 1,]
length(unique(matrix_count[,4]))
sum(matrix_count[,4] > matrix_count[,3])

# Fifth column
bulk_5 = row_num[5] # 24
index_5 = bulk_5 -1
for(i in 1:8){
        env2_index = (i-1) * bulk_2
        chosen_env_2 = c(3, i + 3)
        env_remaining_2 = all_env[which(!all_env %in% chosen_env_2)]
        for(j in 1:7){
                env3_index = env2_index + (j-1) * bulk_3
                chosen_env_3 = env_remaining_2[j]
                env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                for(k in 1:6){
                        env4_index = env3_index + (k -1) * bulk_4
                        chosen_env_4 = env_remaining_3[k]
                        for(l in 1:5){
                                env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                                environ_5 = length(which(count_summary[,3] == "1" | 
                                                                 count_summary[,i +3] == "1" | 
                                                                 count_summary[, chosen_env_3] == "1"|
                                                                 count_summary[, chosen_env_4] == "1"|
                                                                 count_summary[, env_remaining_4[l]] == "1"))
                                matrix_count[((bulk_5 * l - index_5):(bulk_5 * l) + env4_index), 5] = environ_5    
                        }
                        
                }
                
        }
        
}
matrix_count[40320-24,]
matrix_count[40320 - 24 +1,]
matrix_count[40320,]
sum(matrix_count[,5] > matrix_count[,4])


# Sixth column
bulk_6 = row_num[6] # 6
index_6 = bulk_6 -1
for(i in 1:8){
        env2_index = (i-1) * bulk_2
        chosen_env_2 = c(3, i + 3)
        env_remaining_2 = all_env[which(!all_env %in% chosen_env_2)]
        for(j in 1:7){
                env3_index = env2_index + (j-1) * bulk_3
                chosen_env_3 = env_remaining_2[j]
                env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                for(k in 1:6){
                        env4_index = env3_index + (k -1) * bulk_4
                        chosen_env_4 = env_remaining_3[k]
                        env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                        for(l in 1:5){
                                env5_index = env4_index + (l -1) * bulk_5
                                chosen_env_5 = env_remaining_4[l]
                                for(m in 1:4){
                                        env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                        environ_6 = length(which(count_summary[,3] == "1" | 
                                                                         count_summary[,i +3] == "1" | 
                                                                         count_summary[, chosen_env_3] == "1"|
                                                                         count_summary[, chosen_env_4] == "1"|
                                                                         count_summary[, chosen_env_5] == "1"|
                                                                         count_summary[, env_remaining_5[m]] == "1"))
                                        matrix_count[((bulk_6 * m - index_6):(bulk_6 * m) + env5_index), 6] = environ_6
                                }
                                
                        }
                        
                }
                
        }
        
}
matrix_count[1:7,]
matrix_count[(40320 - 7):40320,]
sum(matrix_count[,6] > matrix_count[,5])


# Seventh column
bulk_7 = row_num[7] # 2
index_7 = bulk_7 -1
for(i in 1:8){
        env2_index = (i-1) * bulk_2
        chosen_env_2 = c(3, i + 3)
        env_remaining_2 = all_env[which(!all_env %in% chosen_env_2)]
        for(j in 1:7){
                env3_index = env2_index + (j-1) * bulk_3
                chosen_env_3 = env_remaining_2[j]
                env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                for(k in 1:6){
                        env4_index = env3_index + (k -1) * bulk_4
                        chosen_env_4 = env_remaining_3[k]
                        env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                        for(l in 1:5){
                                env5_index = env4_index + (l -1) * bulk_5
                                chosen_env_5 = env_remaining_4[l]
                                env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                for(m in 1:4){
                                        env6_index = env5_index + (m -1) * bulk_6
                                        chosen_env_6 = env_remaining_5[m]
                                        for(n in 1:3){
                                                env_remaining_6 = env_remaining_5[which(env_remaining_5 != chosen_env_6)]
                                                environ_7 = length(which(count_summary[,3] == "1" | 
                                                                                 count_summary[,i +3] == "1" | 
                                                                                 count_summary[, chosen_env_3] == "1"|
                                                                                 count_summary[, chosen_env_4] == "1"|
                                                                                 count_summary[, chosen_env_5] == "1"|
                                                                                 count_summary[, chosen_env_6] == "1"|
                                                                                 count_summary[, env_remaining_6[n]] == "1"))
                                                matrix_count[((bulk_7 * n - index_7):(bulk_7 * n) + env6_index), 7] = environ_7
                                        }
                                        
                                }
                                
                        }
                        
                }
                
        }
        
}
matrix_count[1:7,]
matrix_count[40300:40320,]
sum(matrix_count[,7] > matrix_count[,6])

# Eighth column
bulk_8 = row_num[8] # 1
index_8 = bulk_8 -1
for(i in 1:8){
        env2_index = (i-1) * bulk_2
        chosen_env_2 = c(3, i + 3)
        env_remaining_2 = all_env[which(!all_env %in% chosen_env_2)]
        for(j in 1:7){
                env3_index = env2_index + (j-1) * bulk_3
                chosen_env_3 = env_remaining_2[j]
                env_remaining_3 = env_remaining_2[which(env_remaining_2 != chosen_env_3)]
                for(k in 1:6){
                        env4_index = env3_index + (k -1) * bulk_4
                        chosen_env_4 = env_remaining_3[k]
                        env_remaining_4 = env_remaining_3[which(env_remaining_3 != chosen_env_4)]
                        for(l in 1:5){
                                env5_index = env4_index + (l -1) * bulk_5
                                chosen_env_5 = env_remaining_4[l]
                                env_remaining_5 = env_remaining_4[which(env_remaining_4 != chosen_env_5)]
                                for(m in 1:4){
                                        env6_index = env5_index + (m -1) * bulk_6
                                        chosen_env_6 = env_remaining_5[m]
                                        env_remaining_6 = env_remaining_5[which(env_remaining_5 != chosen_env_6)]
                                        for(n in 1:3){
                                                env7_index = env6_index + (n -1) * bulk_7
                                                chosen_env_7 = env_remaining_6[n]
                                                for(o in 1:2){
                                                        env_remaining_7 = env_remaining_6[which(env_remaining_6 != chosen_env_7)]
                                                        environ_8 = length(which(count_summary[,3] == "1" | 
                                                                                         count_summary[,i +3] == "1" | 
                                                                                         count_summary[, chosen_env_3] == "1"|
                                                                                         count_summary[, chosen_env_4] == "1"|
                                                                                         count_summary[, chosen_env_5] == "1"|
                                                                                         count_summary[, chosen_env_6] == "1"|
                                                                                         count_summary[, chosen_env_7] == "1"|
                                                                                         count_summary[, env_remaining_7[o]] == "1"))
                                                        matrix_count[((bulk_8 * o - index_8):(bulk_8 * o) + env7_index), 8] = environ_8 
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
matrix_count[,9] = nrow(count_summary)
matrix_count[40311: 40320,]
colnames(matrix_count) = c("One", "Two", "Three", "Four", "Five", "Six", "Seven", "Eight", "Nine")
csvWriter(matrix_count, "Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order.csv")

setwd("~/Dropbox/PPiSeq_02/")
matrix_count = csvReader_T("Working_data/Positive_PPI_environment/PPI_number_environment_accumulation_all_order.csv")
## fit the data with regression
matrix_count_vec = as.numeric(as.vector(t(matrix_count)))
num_env = rep(1:9, nrow(matrix_count))
# linear regression
relation = lm(matrix_count_vec ~ num_env)
linear_count = relation$coefficients[2]* (1:9) + relation$coefficients[1]
# Bionomial regression
polynomial_relation = lm(matrix_count_vec ~ poly(num_env, 2, raw=TRUE))
polynomial_count = polynomial_relation$coefficients[1] + polynomial_relation$coefficients[2] * (1:9) +
        polynomial_relation$coefficients[3] * ((1:9)^2)
# Trinomial regression
polynomial_relation_3 = lm(matrix_count_vec ~ poly(num_env, 3, raw=TRUE))
polynomial_count_3 = polynomial_relation_3$coefficients[1] + polynomial_relation_3$coefficients[2] * (1:9) +
        polynomial_relation_3$coefficients[3] * ((1:9)^2) +  polynomial_relation_3$coefficients[4] * ((1:9)^3)
# Mean count
mean_count = colMeans(matrix_count)

pdf("Working_figure/Figure2/Figure2B_PPI_number_accumulation_environments.pdf", width= 5, height =5)
plot(1:9, as.numeric(matrix_count[1,]), xlim = c(1,9), ylim = c(4000,16000), type = "l",
     col = rgb(0.5,0.5,0.5,alpha = 0.05),lwd = 0.05,
     xlab = "Number of environments assayed",
     ylab = "Number of PPIs", main = NA)
axis(2, at= seq(4000, 16000, by = 2000), labels = seq(4000, 16000, by = 2000))
axis(1, at= 1:9, labels = 1:9)
for(i in 2:nrow(matrix_count)){
        lines(1:9, as.numeric(matrix_count[i,]), col = rgb(0.5,0.5,0.5, alpha = 0.05), lwd = 0.05)
}
lines(1:9, mean_count, col = "red", lwd = 2)
lines(c(1,9), as.numeric(matrix_count[1,c(1, 9)]), col = "orange", lty = 1)
lines(1:9, linear_count, col = "green")
lines(1:9, polynomial_count, col = "cyan", lty = 2)
lines(1:9, polynomial_count_3, col = "blue", lty = 2)
legend("topleft", c("Mean count", "Linear regression", "Bionomial regression", "Trinomial regression", "Straight line"),
       lty = c(1,1,2,2,1), col = c("red", "green", "cyan", "blue", "orange"), bty= "n")
dev.off()

