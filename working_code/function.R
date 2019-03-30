csvReader_T <- function(csv_file){
        return(as.matrix(read.csv(csv_file, header = T)))
}

csvReader_F <- function(csv_file){
        return(as.matrix(read.csv(csv_file, header = F)))
}

tableReader_T <- function(csv_file){
        return(read.table(csv_file, header= T, sep= ","))
}

dataFrameReader_T <- function(csv_file){
        return(read.csv(csv_file, header= T, sep= ","))
}

dataFrameReader_F <- function(csv_file){
        return(read.csv(csv_file, header= F, sep= ","))
}

csvWriter <- function(x, file_name){
        write.csv(x, file_name, row.names=F, quote=F)
}

csvWriter_rownames <- function(x, file_name){
        write.csv(x, file_name,  quote=F)
}

txtWriter <- function (x, file_name){
        write.table(x, file_name, row.names = F, col.names = F, quote = F)
}

frequency <- function(x){
        rs= sum(x)
        return(x/rs)
}

split_string <- function(x){
        a=unlist(strsplit(as.character(x), "_")) 
        return (c(a[1], a[2]))
}

split_string_line <- function(x){
        a=unlist(strsplit(as.character(x), "-")) 
        return (c(a[1], a[2]))
}


split_string_comma <- function(x){
        a= unlist(strsplit(x, ","))
        return(a[1:length(a)])
}

split_string_vector <- function(x){
        PPI= matrix(0, length(x), 2)
        for(i in 1: length(x)){
                a= unlist(strsplit(as.character(x[i]), "_"))
                PPI[i,1]= a[1]
                PPI[i,2]= a[2]
        }
        return (PPI)
}

split_string_vector_line <- function(x){
        PPI= matrix(0, length(x), 2)
        for(i in 1: length(x)){
                a= unlist(strsplit(as.character(x[i]), "-"))
                PPI[i,1]= a[1]
                PPI[i,2]= a[2]
        }
        return (PPI)
}

split_string_vector_comma <- function(x){
        PPI= matrix(0, length(x), 2)
        for(i in 1: length(x)){
                a= unlist(strsplit(as.character(x[i]), ","))
                PPI[i,1]= a[1]
                PPI[i,2]= a[2]
        }
        return (PPI)
}

split_string_vector_ata <- function(x){
    PPI= matrix(0, length(x), 2)
    for(i in 1: length(x)){
        a= unlist(strsplit(as.character(x[i]), "@"))
        PPI[i,1]= a[1]
        PPI[i,2]= a[2]
    }
    return (PPI)
}

split_string_vector_dot <- function(x){
    PPI= matrix(0, length(x), 2)
    for(i in 1: length(x)){
        a= unlist(strsplit(as.character(x[i]), ".", fixed = TRUE))
        PPI[i,1]= a[1]
        PPI[i,2]= a[2]
    }
    return (PPI)
}


protein_count_PPI = function(PPI_vector, threshold){
        all_protein = split_string_vector(PPI_vector)
        all_protein_fr = as.data.frame(table(as.character(c(all_protein[,1], all_protein[,2]))))
        all_protein_fr_order = all_protein_fr[order(all_protein_fr[,2], decreasing = T),]
        all_protein_fr_order_large = all_protein_fr_order[which(all_protein_fr_order[,2] >= threshold),]
        return(all_protein_fr_order_large)
}

# Convert 1e6 into 1 X 10^+6
fancy_scientific <- function(l) {
        # turn in to character string in scientific notation
        l <- format(l, scientific = TRUE)
        l <- gsub("0e\\+00","0",l)
        # quote the part before the exponent to keep all the digits
        l <- gsub("^(.*)e", "'\\1'e", l)
        # turn the 'e+' into plotmath format
        l <- gsub("e+", "%*%10^", l)
        l <- gsub("e-", "%*%10^-", l)
        # return this as an expression
        parse(text=l)
}

### Convert a vector of standard genes to their systematic names; database[,1] = systematic names, database[,2]= standard names
standard_convert_systematic <- function(standard, database){
        database= csvReader_T(database)
        systematic = database[match(standard, database[,2]),1]
        return (systematic)
}

### Take the overlapped systematic and standard genes
Find_standard <- function(gene_list, output){
        gene_list= as.matrix(read.delim(gene_list, header=F))
        systematic_standard = csvReader_T("/Volumes/Zhimin/PPiseq/DMSO/new_pipeline/reference_set/Systematic_standard_protein.csv")
        gene_standard= systematic_standard[which(systematic_standard[,1] %in% gene_list[,1]),]
        write.csv(gene_standard, output, row.names = F,  quote=F)
}

convert_PPI_unique_gene <- function (x){
        PPI= split_string_vector(x)
        unique_gene= unique(c(PPI[,1], PPI[,2]))
        return (unique_gene)
}
### extract the PPI containing the one protein
extract_repeat_PPI= function(x){
        repeat_PPI= vector(mode= "character", length = 0)
        for (i in 1:length(x)){
                a= split_string(x[i])
                if (a[1] == a[2]){
                        repeat_PPI= c(repeat_PPI, x[i])
                }
        }
        return(repeat_PPI)
}


Combine <- function(R1, R2){
        overlap = intersect(R1[,1], R2[,1])
        R1_match = R1[match(overlap, R1[,1]),]
        R2_match = R2[match(overlap, R2[,1]),]
        counts_overlap= as.numeric(R1_match[,2]) + as.numeric(R2_match[,2])
        counts_overlap_dedup= as.numeric(R1_match[,3]) + as.numeric(R2_match[,3])
        R_overlap= cbind(overlap, counts_overlap, counts_overlap_dedup)
        R1_non_match= R1[which(!R1[,1] %in% overlap),]
        R2_non_match= R2[which(!R2[,1] %in% overlap),]
        R_final= rbind(R_overlap, R1_non_match, R2_non_match)
        colnames(R_final) = c("BC1cluster_BC2cluster", "all_reads", "deduped_reads")
        return(R_final)
        
} #### Combine two cluster results from two Hiseq runs

lineage_combine= function(G0, G3, G6, G9, G12, lineage){
        barcodes= unique(c(G0[,1], G3[,1], G6[,1], G9[,1], G12[,1]))
        G0_matched= G0[match(barcodes, G0[,1]), 3]
        G3_matched= G3[match(barcodes, G3[,1]), 3]
        G6_matched= G6[match(barcodes, G6[,1]), 3]
        G9_matched= G9[match(barcodes, G9[,1]), 3]
        G12_matched= G12[match(barcodes, G12[,1]), 3]
        matrix_final= cbind(barcodes, G0_matched, G3_matched, G6_matched, G9_matched, G12_matched)
        matrix_final[is.na(matrix_final)] = 0
        colnames(matrix_final)= c("Barcodes", "G0", "G3", "G6", "G9", "G12")
        csvWriter(matrix_final, lineage)
}

lineage_combine_4= function(G0, G6, G9, G12, lineage){
        barcodes= unique(c(G0[,1], G6[,1], G9[,1], G12[,1]))
        G0_matched= G0[match(barcodes, G0[,1]), 3]
        G6_matched= G6[match(barcodes, G6[,1]), 3]
        G9_matched= G9[match(barcodes, G9[,1]), 3]
        G12_matched= G12[match(barcodes, G12[,1]), 3]
        matrix_final= cbind(barcodes, G0_matched, G6_matched, G9_matched, G12_matched)
        matrix_final[is.na(matrix_final)] = 0
        colnames(matrix_final)= c("Barcodes", "G0", "G6", "G9", "G12")
        csvWriter(matrix_final, lineage)
}



# 1st column of barcode_matrix and 2nd column of PPI_matrix are barcodes
barcodes_match_PPI <- function(barcodes_matrix, PPI_matrix){
        barcodes_matrix_matched= barcodes_matrix[match(PPI_matrix[,3], barcodes_matrix[,1]),]
        PPI_barcodes= cbind(PPI_matrix[,1:2], barcodes_matrix_matched)
        PPI_barcodes= PPI_barcodes[!is.na(PPI_barcodes[,3]),]
        return(PPI_barcodes)
}

### Match barcode to known PPIs. The sepcific folder contain all the PPIs and their corresponding barcodes
PPI_match= function(lineage_trajectories, known_PPI_barcode, unmatched_barcode_lineage, unmatched_PPI_barcodes){
        lineage_combine= csvReader_T(lineage_trajectories)
        PPI_barcodes_all = csvReader_T("/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/PPI_barcodes/PPiseq_all_barcodes.csv")
        DBC_known_counts = barcodes_match_PPI(lineage_combine, PPI_barcodes_all)#6126314
        environment= rep("DMSO", nrow(DBC_known_counts))
        DBC_known_final= cbind(DBC_known_counts[,1:3], environment, DBC_known_counts[, 4:8])
        colnames(DBC_known_final)= c("Systematic_PPI", "Standard_PPI", "barcode", "environment", "G0", "G3", "G6", "G9", "G12")
        csvWriter(DBC_known_final, known_PPI_barcode)
        unmatched_DBC = lineage_combine[which(!lineage_combine[,1] %in% DBC_known_final[,3]),]
        csvWriter(unmatched_DBC, unmatched_barcode_lineage)
        unmatched_PPI_barcodes = PPI_barcodes_all[which(!PPI_barcodes_all[,3] %in% DBC_known_final[,3]),]
        csvWriter(unmatched_PPI_barcodes, unmatched_PPI_barcodes)
        
}


### Match barcodes to known PPIs with allowing 2 mismatches for each double barcodes; take too much time

mismatches = function(seq1, seq2, max){
        if (seq1 != seq2){
                a= unlist(strsplit(seq1,""))
                b= unlist(strsplit(seq2,""))
                mismatch= abs(length(a)- length(b))
                for(i in 1:min(length(a),length(b))){
                        if(a[i] != b[i]){
                                mismatch = mismatch + 1
                        }
                        if (mismatch >= max){
                                break
                        }
                }
                return(mismatch)
        }
        else{
                return (mismatch)}
}

barcodes_mismatch_PPI <- function(barcodes_matrix, PPI_matrix){
        new_matrix = matrix(0, nrow(barcodes_matrix), ncol(barcodes_matrix)+2)
        for(i in 1:nrow(barcodes_matrix)){
                dist = rep(10, nrow(PPI_matrix))
                for (j in 1: nrow(PPI_matrix)){
                        dist[j] = mismatches(barcodes_matrix[i,1], PPI_matrix[j, 3], 3)
                }
                min_dist_index = which(dist == min(dist))
                if (min(dist) < 3 & length(min_dist_index) == 1){
                        new_matrix[i, 1:2] = PPI_matrix[min_dist_index, 1:2]
                        new_matrix[i, 3: ncol(new_matrix)] = barcodes_matrix[i,]
                }
                else{new_matrix[i,3] = "NO"}
        }
        new_matrix= new_matrix[which(new_matrix[,3] != "NO"),]
        return(new_matrix)
}

#### Replace counts with frequencies and write into a frequency matrix
Freq_convert = function(input, column_of_count, output){
    input_matrix = csvReader_T(input)
    out_matrix = matrix(0, nrow(input_matrix), ncol(input_matrix))
    out_matrix[,1:(column_of_count -1)] = input_matrix[,1:(column_of_count - 1)]
    #out_matrix[,(column_of_count + 5) : ncol(input_matrix)] = input_matrix[, (column_of_count + 5) : ncol(input_matrix)]
    for (i in column_of_count: (column_of_count +4)){ # only calculate the frequency with the columns of counts (5:9)
        out_matrix[,i] = frequency(input_matrix[,i])
    }
    colnames(out_matrix) = colnames(input_matrix)
    csvWriter(out_matrix, output)
    return(out_matrix)
}

##### combine lineage trajectory counts to each PPI or barcode 
add_trajectory_PPI_fitness <- function(PPI_estimated_fitness, known_PPI_barcode_counts, output){
        PPI_barcode= csvReader_F(PPI_estimated_fitness)
        barcode_counts= csvReader_T(known_PPI_barcode_counts)
        matched_trajectory= barcode_counts[match(PPI_barcode[,3], barcode_counts[,3]), 5:9]
        PPI_barcode_trajectory= cbind(PPI_barcode, matched_trajectory)
        colnames(PPI_barcode_trajectory) = c("PPI", "barcode_kinds", "barcode", "mean_fitness", "likelihoods", "distance", "G0_counts", "G3_counts", "G6_counts", "G9_counts", "G12_counts")
        csvWriter(PPI_barcode_trajectory, output)
        
} 

#### calculate the mean fitness of control strains (Negative or positive strains)
mean_fitness_calculation <- function(DHFR_negative_PPI){
        neg_freq= frequency(as.numeric(DHFR_negative_PPI[,7]))
        neg_fitness= as.numeric(DHFR_negative_PPI[,4])
        individual_fitness = rep(0, length(neg_freq))
        for (i in 1:length(neg_freq)){
                individual_fitness[i]= neg_freq[i] * neg_fitness[i]
        }
        mean_fitness = sum(individual_fitness)
        return(mean_fitness)
}

### Calculate the P-values for PPIs marked with different barcodes. 
# Input: relative fitness file; Output 1 barcode PPI data and multiple barcodes PPI data with pvalues
P_value <- function(PPI_relative_fitness, output_01, output_02){
        PPI_indiv= csvReader_T(PPI_relative_fitness)
        PPI_negative= PPI_indiv[which(PPI_indiv[,1] == "negative_non_DHFR"),]
        PPI_negative_fitness= as.numeric(PPI_negative[,4])
        PPI_unique= unique(PPI_indiv[,1])
        matrix_summary= matrix(0, length(PPI_unique), 6)
        matrix_summary[,1]= PPI_unique
        colnames(matrix_summary)= c("PPI", "barcode_counts", "mean", "sd", "p.value", "FDR_adjusted_value")
        b=1
        for (i in 1: nrow(matrix_summary)){
                a= as.numeric(PPI_indiv[b,2])
                if(a > 1){
                        matrix_summary[i,2]= a
                        matrix_summary[i,3]= mean(as.numeric(PPI_indiv[b:(b+a-1), 4]))
                        matrix_summary[i,4]= sd(as.numeric(PPI_indiv[b:(b+a-1),4]))
                        matrix_summary[i,5]= t.test(as.numeric(PPI_indiv[(b:(b+a-1)),4]), PPI_negative_fitness, alternative = "greater")$p.value
                        b= b+a
                }
                else{
                        matrix_summary[i, 2]= a
                        matrix_summary[i, 3]= as.numeric(PPI_indiv[b,4])
                        matrix_summary[i, 4]= "NA"
                        matrix_summary[i, 5]= "NA"
                        b= b+a
                }
                
        }
        matrix_PPI_01= matrix_summary[which(matrix_summary[,5] == "NA"),]
        matrix_filtered= matrix_summary[which(!matrix_summary[,5] == "NA"),]
        matrix_filtered[,6]= p.adjust(as.numeric(matrix_filtered[,5]), method= "fdr")
        
        csvWriter(matrix_PPI_01, output_01)
        csvWriter(matrix_filtered, output_02)
}

P_value_all_negative <- function(PPI_relative_fitness, output_01, output_02){
        PPI_indiv= csvReader_T(PPI_relative_fitness)
        MATa_DHFR12 = csvReader_T("/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/PPI_barcodes/MATa_genome_combine.csv")
        MATalpha_DHFR3 = csvReader_T("/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/PPI_barcodes/MATalpha_genome_combine.csv")
        PPI_DHFR12 = PPI_indiv[which(PPI_indiv[,3] %in% MATa_DHFR12[,3]),]
        PPI_DHFR3 = PPI_indiv[which(PPI_indiv[,3] %in% MATalpha_DHFR3[,3]),]
        PPI_negative_DHFR= PPI_indiv[which(PPI_indiv[,1] == "negative_non_DHFR"),]
        PPI_negative = rbind(PPI_DHFR12, PPI_DHFR3, PPI_negative_DHFR)
        PPI_negative_fitness= as.numeric(PPI_negative[,4])
        
        PPI_unique= unique(PPI_indiv[,1])
        matrix_summary= matrix(0, length(PPI_unique), 6)
        matrix_summary[,1]= PPI_unique
        colnames(matrix_summary)= c("PPI", "barcode_counts", "mean", "sd", "p.value", "FDR_adjusted_value")
        b=1
        for (i in 1: nrow(matrix_summary)){
                a= as.numeric(PPI_indiv[b,2])
                if(a > 1){
                        matrix_summary[i,2]= a
                        matrix_summary[i,3]= mean(as.numeric(PPI_indiv[b:(b+a-1), 4]))
                        matrix_summary[i,4]= sd(as.numeric(PPI_indiv[b:(b+a-1),4]))
                        matrix_summary[i,5]= t.test(as.numeric(PPI_indiv[(b:(b+a-1)),4]), PPI_negative_fitness, alternative = "greater")$p.value
                        b= b+a
                }
                else{
                        matrix_summary[i, 2]= a
                        matrix_summary[i, 3]= as.numeric(PPI_indiv[b,4])
                        matrix_summary[i, 4]= "NA"
                        matrix_summary[i, 5]= "NA"
                        b= b+a
                }
                
        }
        matrix_PPI_01= matrix_summary[which(matrix_summary[,5] == "NA"),]
        matrix_filtered= matrix_summary[which(!matrix_summary[,5] == "NA"),]
        matrix_filtered[,6]= p.adjust(as.numeric(matrix_filtered[,5]), method= "fdr")
        
        csvWriter(matrix_PPI_01, output_01)
        csvWriter(matrix_filtered, output_02)
}
### Input: a vector of PPI; output: a matrix in which the first column is unique PPI, 
### and the second column is the oppsite PPI of the first column. If no opposite PPI, the second column is Zero 
mark_duplicates <- function(x){
        matrix_all= matrix(0, length(x), 2)
        # first label PPIs that have duplicates of opposite direction (two columns)
        for (i in 1:length(x)){
                a= split_string(x[i])
                b= paste(a[2], a[1], sep="_") # generate the PPI with opposite direction
                PPI_remaining= x[which(x != x[i])]
                c= which(PPI_remaining == b) # check the duplicates in the list
                if (length(c) == 0){
                        matrix_all[i,1] = x[i]
                }
                else{
                        matrix_all[i,1] = x[i]
                        matrix_all[i,2] = PPI_remaining[c]
                }
        }
        if(length(unique(matrix_all[,2])) > 1){
                matrix_all_duplicate= matrix_all[which(matrix_all[,2] !=0),]
                PPI_unique= matrix_all_duplicate[1,1] # store the first PPI in the matrix
                matrix_all_duplicate_half= matrix(0, nrow(matrix_all_duplicate), ncol(matrix_all_duplicate))
                matrix_all_duplicate_half[1,] = matrix_all_duplicate[1,]
                for (i in 2: nrow(matrix_all_duplicate_half)){
                        if (!matrix_all_duplicate[i,2] %in% PPI_unique){
                                matrix_all_duplicate_half[i,] = matrix_all_duplicate[i,]
                                PPI_unique= c(PPI_unique, matrix_all_duplicate[i,1])
                        }
                        else{
                                matrix_all_duplicate_half[i,]= c(0,0)}
                }
                matrix_all_duplicate_half= matrix_all_duplicate_half[which(matrix_all_duplicate_half[,1] != 0),]
                matrix_all_non_duplicate= matrix_all[which(matrix_all[,2] == 0),]
                matrix_duplicate_marked = rbind(matrix_all_duplicate_half, matrix_all_non_duplicate)
        }
        else{
                matrix_duplicate_marked= matrix_all
        }
        
        return (matrix_duplicate_marked)
        
}

# extract the self-interacting PPIs
extract_repeat_PPI= function(x){
        repeat_PPI= vector(mode= "character", length = 0)
        for (i in 1:length(x)){
                a= split_string(x[i])
                if (a[1] == a[2]){
                        repeat_PPI= c(repeat_PPI, x[i])
                }
        }
        return(repeat_PPI)
}

### Input: a vector of PPI; output: a matrix in which the first column is unique PPI, 
### and the second column is the oppsite PPI of the first column. If no opposite PPI, the second column is Zero

mark_duplicates_fast <- function(x){
        repeat_PPI= extract_repeat_PPI(x)
        repeat_PPI_matrix= cbind(repeat_PPI, rep(0, length(repeat_PPI)))
        x= x[which(!x %in% repeat_PPI)]
        x_reverse= rep(0, length(x))
        for (i in 1:length(x)){
                a= split_string(x[i])
                b= paste(a[2], a[1], sep="_") # generate the PPI with opposite direction
                x_reverse[i] = b
        }
        matrix_all= matrix(0, length(x), 2)
        matrix_all[,1]= x
        matrix_all[,2]= x[match(x, x_reverse)]
        matrix_all[is.na(matrix_all)] = 0
        
        if(length(unique(matrix_all[,2])) > 1){
                matrix_all_duplicate= matrix_all[which(matrix_all[,2] !="0"),]
                all_duplicate = as.vector(t(matrix_all_duplicate))
                all_duplicate_unique= unique(all_duplicate)
                matrix_all_duplicate_order= matrix(all_duplicate_unique, length(all_duplicate_unique)/2, byrow = TRUE)
                matrix_all_non_duplicate= matrix_all[which(matrix_all[,2] == "0"),]
                matrix_duplicate_marked = rbind(matrix_all_duplicate_order, matrix_all_non_duplicate, repeat_PPI_matrix)
        }
        else{
                matrix_duplicate_marked= rbind(matrix_all, repeat_PPI_matrix)
        }
        return (matrix_duplicate_marked)
}

# Match the known PPI matrix with PPI sets, output the PPI matrix that is found in the PPI sets (two directions)
match_both_direction <- function (PPiseq_matrix, known_PPI){
        a= rep(0, length(known_PPI))
        b= rep(0, length(known_PPI))
        for (i in 1:length(known_PPI)){
                a[i]= split_string(known_PPI[i])[1]
                b[i]= split_string(known_PPI[i])[2]
        }
        known_PPI_02= paste(b, a, sep="_")
        matched_PPI = PPiseq_matrix[which(PPiseq_matrix[,1] %in% known_PPI | PPiseq_matrix[,1] %in% known_PPI_02),]
        return (matched_PPI)
}

#### Check double GO terms enriched for each PPI. 
## Input: PPI, GO terms result (each GO terms: genes), name conversion; Output: counts of PPI for each of double GO terms enriched
map_PPI_GO <- function (PPI_not_reported, Function_slim_result, standard_name_convert, output_matrix, out_put_matrix_ratio){
        PPI_overlap= csvReader_T(PPI_not_reported)[,1]
        result_overlap = csvReader_T(Function_slim_result)
        cell_component_overlap= result_overlap[,2]
        cell_component_overlap = gsub(",", "-", cell_component_overlap)
        protein_overlap= list()
        for (i in 1:length(cell_component_overlap)){
                protein_overlap[[i]]= standard_convert_systematic(split_string_comma(result_overlap[i,5]), standard_name_convert)
        }
        matrix_overlap= matrix(0, length(cell_component_overlap), length(cell_component_overlap))
        colnames(matrix_overlap) = cell_component_overlap
        rownames(matrix_overlap) = cell_component_overlap
        matrix_overlap_ratio= matrix(0, length(cell_component_overlap), length(cell_component_overlap))
        colnames(matrix_overlap_ratio) = cell_component_overlap
        rownames(matrix_overlap_ratio) = cell_component_overlap
        for (i in 1:nrow(matrix_overlap)){
                for(j in 1:ncol(matrix_overlap)){
                        genes_1= protein_overlap[[i]]
                        genes_2= protein_overlap[[j]]
                        a= 0
                        for(k in 1:length(PPI_overlap)){
                                PPI= unlist(split_string(PPI_overlap[k]))
                                if (PPI[1] %in% genes_1 & PPI[2] %in% genes_2){
                                        a= a + 1
                                }
                                else if (PPI[1] %in% genes_2 & PPI[2] %in% genes_1){
                                        a= a +1
                                }
                                else {a= a}
                        }
                        matrix_overlap[i,j] = a
                        matrix_overlap_ratio[i,j] = a/length(PPI_overlap)
                        
                }
        }
        csvWriter_rownames(matrix_overlap, output_matrix)
        csvWriter_rownames(matrix_overlap_ratio, output_matrix_ratio)
}

##### Transform lineage matrix to data frame which can be plotted in ggplot2
transform_lineages <- function(x, PPI_pos, barcode_pos, fitness_pos, G0_pos){
        time_points = rep(c(0,3,6,9,12), nrow(x))
        lineages = c(t(x[, G0_pos:(G0_pos + 4)]))
        PPI = rep("0", length(time_points))
        barcode = rep("0", length(time_points))
        fitness = rep(0, length(time_points))
        for(i in 1:nrow(x)){
                PPI[(i*5 -4): (i*5)] = as.character(rep(x[i,PPI_pos], 5))
                barcode[(i*5 -4): (i*5)] = as.character(rep(x[i,barcode_pos], 5))
                fitness[(i*5 -4): (i*5)] = rep(x[i,fitness_pos], 5)
        }
        matrix_lineages = data.frame(PPI, barcode, fitness, time_points, lineages)
        csvWriter(matrix_lineages, "Lineage_trajectories_plot_data_frame.csv")
        return (matrix_lineages)
}

#### Check the counts of index hopping

#index_to_be_checked file format (4 column matrix):
#Forward_primer_name Reverse_primer_name Forward_primer_sequence Reverse_primer_seq
#p204 P216 CATAGG GTATTG

IndexHopping <- function(unmatched_multitag, index_to_be_checked, output_seq, output_counts){
        primer_seq = csvReader_T(index_to_be_checked)
        matrix_multitag = as.data.frame(table(read.delim(unmatched_multitag)))
        matrix = matrix(0, nrow(primer_seq), nrow(primer_seq))
        rownames(matrix) = primer_seq[,1]
        colnames(matrix) = primer_seq[,2]
        matrix_counts = matrix(0, nrow(primer_seq), nrow(primer_seq))
        rownames(matrix_counts) = primer_seq[,1]
        colnames(matrix_counts) = primer_seq[,2]
        for (i in 1:nrow(matrix)){
                for(j in 1: nrow(matrix)){
                        matrix[i,j] = paste(primer_seq[i,3], primer_seq[j,4], sep = "")
                        matrix_counts[i,j] = matrix_multitag[match(matrix[i,j], as.character(matrix_multitag[,1])),2]
                }
        }
        
        csvWriter_rownames(matrix, output_seq)
        csvWriter_rownames(matrix_counts, output_counts)
        
}

PPV_matrix = function(matrix_ref, output){
        number_Pos = length(which(matrix_ref[,1] == "Positive")) 
        number_Neg = length(which(matrix_ref[,1] == "Negative")) 
        Fitness = seq(0, 1, by= 0.01) # 51
        Q_values = seq(-20, -1.5, by = 0.5) # 38
        PPV_matrix = matrix(0, length(Fitness), length(Q_values))
        rownames(PPV_matrix) = as.character(Fitness)
        colnames(PPV_matrix) = as.character(Q_values)
        for (i in 1:length(Fitness)){
                for (j in 1:length(Q_values)){
                        positive_threshold = matrix_ref[which(matrix_ref[,3] >= Fitness[i] &
                                                                      log10(matrix_ref[,6]) <= Q_values[j]),1]
                        true_positive = length(which(positive_threshold == "Positive"))
                        false_positive = length(which(positive_threshold == "Negative"))
                        PPV_matrix[i, j] = true_positive/(true_positive + false_positive)
                }
        }
        
        csvWriter_rownames(PPV_matrix, output )
        return(PPV_matrix)
}

Random_reference <- function(PPI_multiple, all_PPI_filtered, size, sampling_number){
        PPI_PPiseq = PPI_multiple[,1]
        RRS_size= 0
        RRS= character(length=0)
        while (RRS_size < size) {
                yeast_PPI_random= all_PPI_filtered[sample(1: length(all_PPI_filtered), sampling_number, replace=F)]
                RRS= unique(c(RRS, yeast_PPI_random))
                RRS_overlap= intersect(RRS, PPI_PPiseq)
                RRS_size= length(RRS_overlap)
        } 
        RRS_duplicate_marked = mark_duplicates_fast(RRS_overlap)
        return (RRS_duplicate_marked[,1])
}

Random_reference_generation = function(PPI_multiple, all_PPI_filtered, number_PPI, sample_number, data_set_number){
        neg_set = vector("list", data_set_number)
        for(i in 1:data_set_number){
                PPI_neg = Random_reference(PPI_multiple, all_PPI_filtered, number_PPI, sample_number)
                PPI_neg_matrix = PPI_multiple[which(PPI_multiple[,1] %in% PPI_neg),]
                neg_set[[i]] = PPI_neg_matrix
        }
        return (neg_set)
}

PPV_matrix = function(matrix_ref, output){
        number_Pos = length(which(matrix_ref[,1] == "Positive")) 
        number_Neg = length(which(matrix_ref[,1] == "Negative")) 
        Fitness = seq(0, 1, by= 0.01) # 51
        Q_values = seq(-20, -1.5, by = 0.5) # 38
        PPV_matrix = matrix(0, length(Fitness), length(Q_values))
        rownames(PPV_matrix) = as.character(Fitness)
        colnames(PPV_matrix) = as.character(Q_values)
        for (i in 1:length(Fitness)){
                for (j in 1:length(Q_values)){
                        positive_threshold = matrix_ref[which(matrix_ref[,3] >= Fitness[i] &
                                                                      log10(matrix_ref[,6]) <= Q_values[j]),1]
                        true_positive = length(which(positive_threshold == "Positive"))
                        false_positive = length(which(positive_threshold == "Negative"))
                        PPV_matrix[i, j] = true_positive/(true_positive + false_positive)
                }
        }
        
        csvWriter_rownames(PPV_matrix, output )
        return(PPV_matrix)
}

PPV_threshold = function(PPV_matrix_ref, PPV_threshold, index){
        matrix_PPV_threshold = matrix(NA, ncol(PPV_matrix_ref), 3)
        for (i in 1: ncol(PPV_matrix_ref)){
                matrix_PPV_threshold[i,2] = colnames(PPV_matrix_ref)[i]
                for (j in 1: nrow(PPV_matrix_ref)){
                        if (is.na(PPV_matrix_ref[j,i])){
                                next
                        }
                        else {
                                if (PPV_matrix_ref[j,i] >= PPV_threshold){
                                        matrix_PPV_threshold[i,3] = rownames(PPV_matrix_ref_01)[j]
                                        break
                                }  
                        }
                }
        }
        matrix_PPV_threshold [,1]= paste("data",as.character(index), sep = "_")
        return(matrix_PPV_threshold)
}

ROC_matrix = function(matrix_ref, selected_fitness, selected_Q_value, pos_neg_ratio, index){
        number_Pos = length(which(matrix_ref[,1] == "Positive")) 
        number_Neg = length(which(matrix_ref[,1] == "Negative")) 
        Fitness = seq(0, 1, by= 0.02) # 51
        Q_values = seq(-20, -1.5, by = 0.5) # 38
        l_f = length(Fitness)
        l_q = length(Q_values)
        row_number = l_f * l_q
        ROC_matrix = matrix(0, row_number, 7)
        colnames(ROC_matrix) = c("Index", "Ratio", "Fitness", "Q_value", "TPR", "FPR", "PPV")
        if (selected_fitness != 0 & selected_Q_value == 0){
                for(i in 1:l_f){
                        ROC_matrix[(i *l_q - (l_q-1)): (i*l_q),3] = rep(Fitness[i], l_q)
                        ROC_matrix[(i *l_q - (l_q-1)): (i*l_q),4] = Q_values
                }
        }else if (selected_fitness == 0 & selected_Q_value != 0){
                for(i in 1:l_q){
                        ROC_matrix[(i *l_f - (l_f-1)): (i*l_f),3] = Fitness
                        ROC_matrix[(i *l_f - (l_f-1)): (i*l_f),4] = rep(Q_values[i], l_f)
                }
        }
        
        for(i in 1:nrow(ROC_matrix)){
                positive_threshold = matrix_ref[which(matrix_ref[,3] >= ROC_matrix[i,3] &
                                                              log10(matrix_ref[,6]) <= ROC_matrix[i,4]),1]
                true_positive = length(which(positive_threshold == "Positive"))
                false_positive = length(which(positive_threshold == "Negative"))
                ROC_matrix[i, 5] = true_positive/number_Pos
                specificity = (number_Neg - false_positive)/number_Neg
                ROC_matrix[i, 6] = 1 - specificity
                ROC_matrix[i, 7] = true_positive/(true_positive + false_positive)
                
        }
        if (selected_fitness != 0 & selected_Q_value == 0){
                ROC_matrix_selected = ROC_matrix[which(ROC_matrix[,3] == selected_fitness),]
        }else if (selected_fitness == 0 & selected_Q_value != 0){
                ROC_matrix_selected = ROC_matrix[which(ROC_matrix[,4] == selected_Q_value),]
        }
        
        ROC_matrix_selected[,1] = index
        ROC_matrix_selected[,2] = pos_neg_ratio
        return (ROC_matrix_selected)
}

protein_degree_count = function(PPI){
        source("/Volumes/zmliu_02/PPiseq/HU/R_code/function.R")
        all_PPI_gene = split_string_vector(PPI)
        protein_degree = as.data.frame(table(as.character(c(all_PPI_gene[,1], all_PPI_gene[,2]))))
        protein_degree_order= protein_degree[order(protein_degree[,2], decreasing = T),]
        return(protein_degree_order)
}

####### Function to calculate normalized fitness values for each replicated PPI, mean fitness value, Q-value

### First extract positive PPIs, DHFR positive and negative controls from lineage data
dynamic_PPI_residual = function(DMSO_lineage, H2O2_lineage, DMSO_real, H2O2_real, output_file){
        DMSO_DHFR_Pos= DMSO_lineage[which(DMSO_lineage[,1] == "positive_DHFR"),]
        DMSO_DHFR_Neg= DMSO_lineage[which(DMSO_lineage[,1] == "negative_non_DHFR"),]
        
        DMSO_DHFR_Pos_mean = mean(DMSO_DHFR_Pos[,4]) # 0.97388
        DMSO_DHFR_Neg_mean = mean(DMSO_DHFR_Neg[,4]) # 0.0725669
        
        H2O2_DHFR_Pos= H2O2_lineage[which(H2O2_lineage[,1] == "positive_DHFR"),]
        H2O2_DHFR_Neg= H2O2_lineage[which(H2O2_lineage[,1] == "negative_non_DHFR"),]
        
        H2O2_DHFR_Pos_mean = mean(H2O2_DHFR_Pos[,4]) # 0.753049
        H2O2_DHFR_Neg_mean = mean(H2O2_DHFR_Neg[,4]) # 0.01723175
        
        DMSO_H2O2_PPI = unique(c(as.character(DMSO_real[,1]), as.character(H2O2_real[,1]))) # 6654
        
        ### extract lineage data of positive PPIs
        DMSO_lineage_Pos = DMSO_lineage[which(DMSO_lineage[,1] %in% DMSO_H2O2_PPI),c(1,2,4)]
        H2O2_lineage_Pos = H2O2_lineage[which(H2O2_lineage[,1] %in% DMSO_H2O2_PPI),c(1,2,4)]
        ## normalize the fitness with two controls
        DMSO_lineage_Pos[,3] = (DMSO_lineage_Pos[,3] - DMSO_DHFR_Neg_mean)/(DMSO_DHFR_Pos_mean - DMSO_DHFR_Neg_mean)
        H2O2_lineage_Pos[,3] = (H2O2_lineage_Pos[,3] - H2O2_DHFR_Neg_mean)/(H2O2_DHFR_Pos_mean - H2O2_DHFR_Neg_mean)
        ## Consider all the negative fitness values as 0
        #length(which(DMSO_lineage_Pos[,3] < 0)) # 296
        #DMSO_lineage_Pos[which(DMSO_lineage_Pos[,3] < 0),3] = 0
        #length(which(H2O2_lineage_Pos[,3] < 0)) # 591
        #H2O2_lineage_Pos[which(H2O2_lineage_Pos[,3] < 0),3] = 0
        
        DMSO_H2O2_matrix = matrix(0, length(DMSO_H2O2_PPI), 15)
        DMSO_H2O2_matrix[,1] = DMSO_H2O2_PPI
        for(i in 1:nrow(DMSO_H2O2_matrix)){
                index_1 = which(as.character(DMSO_lineage_Pos[,1]) == DMSO_H2O2_matrix[i,1])
                index_2 = which(as.character(H2O2_lineage_Pos[,1]) == DMSO_H2O2_matrix[i,1])
                
                a = mean(DMSO_lineage_Pos[index_1,3])
                b = mean(H2O2_lineage_Pos[index_2,3])
                if(length(index_1) > 0 & a < 0){ a = 0}
                if(length(index_2) > 0 & b < 0){ b = 0}
                DMSO_H2O2_matrix[i,2] = a
                DMSO_H2O2_matrix[i,3] = b
                
                DMSO_H2O2_matrix[i,4] = length(index_1)
                DMSO_H2O2_matrix[i,5] = length(index_2)
                DMSO_H2O2_matrix[i,6] = b - a
                if(length(index_1) > 1 & length(index_2) > 1){
                        DMSO_H2O2_matrix[i,7] = t.test(DMSO_lineage_Pos[index_1,3], H2O2_lineage_Pos[index_2,3], alternative = "two.sided")$p.value  
                }else{
                        DMSO_H2O2_matrix[i,7] = 10
                }
                
                if(length(index_1) > 0){
                        DMSO_H2O2_matrix[i, 8:(7 + length(index_1))] = DMSO_lineage_Pos[index_1, 3]
                }
                if(length(index_2) > 0){
                        DMSO_H2O2_matrix[i, 12:(11 + length(index_2))] = H2O2_lineage_Pos[index_2, 3]
                }
                
                
        }
        DMSO_H2O2_matrix[,7]= p.adjust(DMSO_H2O2_matrix[,7], method= "fdr")
        colnames(DMSO_H2O2_matrix) = c("PPI", "mean_DMSO", "mean_other", "Num_replicate_DMSO",
                                       "Num_replicate_other","difference", "Q-value","DMSO_1", "DMSO_2", "DMSO_3", "DMSO_4",
                                       "Other_1", "Other_2", "Other_3", "Other_4")
        csvWriter(DMSO_H2O2_matrix, output_file)
}

***
##### Calculate the p-values for the same PPI in two environments by Cyber-T test. The input matrix should be output of dynamic_PPI_residual (above function)
source("/Volumes/zmliu_02/PPiseq/Combine_environments/one_binary_PPV/dynamic_PPI/cyberT_test/bayesreg.R")
Cal_Bayes_q = function(DMSO_H2O2_all, output_name){
        replicate_number = paste(as.character(DMSO_H2O2_all[,4]), as.character(DMSO_H2O2_all[,5]), sep = "_")
        DMSO_H2O2_all = cbind(DMSO_H2O2_all, replicate_number)
        group = unique(replicate_number)
        data_all = rep(0, (ncol(DMSO_H2O2_all)+1))
        for(i in 1:length(group)){
                group[i]
                data = DMSO_H2O2_all[which(replicate_number == group[i]),]
                if (data[1,4]== 1 | data[1,4] == 0 | data[1,5] == 1 | data[1,5] == 0){
                        bayes_Q= rep(1, nrow(data))
                        data_final = as.matrix(cbind(data, bayes_Q))
                }else{
                        num_rep_01 = data[1,4]
                        num_rep_02 = data[1,5]
                        data_cyber = data[,c(8:(7 + num_rep_01), 12:(11 + num_rep_02))]
                        if (nrow(data_cyber) >= 2000){
                                winSize = 101
                        }else if (nrow(data_cyber) >= 1000 & nrow(data_cyber) < 2000){
                                winSize = 51
                        }else if (nrow(data_cyber < 1000)){
                                temp_size = floor(nrow(data_cyber)/20)
                                if (temp_size == 1 | temp_size == 0) { winSize = 3}
                                else if(temp_size %% 2 == 0){winSize = temp_size + 1}
                                else {winSize = temp_size}
                        }        
                        bayes_Q = bayesT(aData= data_cyber, numC= num_rep_01, numE= num_rep_02, ppde=FALSE, 
                                         betaFit=1, bayes=TRUE, winSize, conf=12, 
                                         doMulttest=TRUE, bayesIntC=FALSE, bayesIntE=FALSE)$BH
                        data_final= as.matrix(cbind(data, bayes_Q))
                }
                data_all = rbind(data_all, data_final)
        }
        data_all = data_all[2:nrow(data_all),]
        csvWriter(data_all, output_name)
        return(data_all)
}
***
#### To show the number in y axis in a scientific way
fancy_scientific <- function(l) {
        # turn in to character string in scientific notation
        l <- format(l, scientific = TRUE)
        # quote the part before the exponent to keep all the digits
        l <- gsub("^(.*)e", "'\\1'e", l)
        # turn the 'e+' into plotmath format
        l <- gsub("e", "%*%10^", l)
        # return this as an expression
        parse(text=l)
}