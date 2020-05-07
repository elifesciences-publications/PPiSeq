def pairwise_mismatch(directory = "/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files_final/counts/",
                      cluster_threshold = 1, # only take the cluster with counts larger than threshold
                      known_barcodes_1 ='/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/PPI_barcodes/protein_barcodes_first.csv',
                      known_barcodes_2 ='/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/PPI_barcodes/protein_barcodes_second.csv',
                      old_index_1 = '/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/Unmatched_barcodes/Unmatched_barcodes_protein_index_01.csv',
                      old_index_2 = '/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/Unmatched_barcodes/Unmatched_barcodes_protein_index_02.csv',
                      old_BC1_cluster = '/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/Unmatched_barcodes/All_time_01_BC1_cluster.csv',
                      old_BC2_cluster = '/Volumes/Zhimin/PPiseq/DMSO/all_lintag_files/Unmatched_barcodes/All_time_01_BC2_cluster.csv',
                      unmatched_lintag1 = 'G0_all_lintag1.txt',
                      unmatched_lintag2 = 'G0_all_lintag2.txt',
                      unmatched_seqtag = 'G0_all_seqtag.txt',
                      bartender_BC1_barcode = 'All_time_BC1_barcode.csv',
                      bartender_BC2_barcode = 'All_time_BC2_barcode.csv',
                      bartender_BC1_cluster = 'All_time_BC1_cluster.csv',
                      bartender_BC2_cluster = 'All_time_BC2_cluster.csv',
                      output= "G0_all_PPIcounts.csv",
                      index_construct = 0):
    import csv
    if index_construct ==1:
        ### First check the overlap with the map of old-cluster and known barcodes, and update the index for new cluster
        with open(directory + bartender_BC1_cluster) as f:
        	reader = csv.reader(f)
        	new_bartender_BC1_cluster=list(reader)
        del new_bartender_BC1_cluster[0]
        new_bartender_BC1_cluster = [i for i in new_bartender_BC1_cluster if int(i[3]) > cluster_threshold]
        with open(directory + 'bartender_cluster_01_filtered.csv', 'w') as file: # output the filtered cluster
        	file.writelines(','.join([a[0], a[1]]) + '\n' for a in new_bartender_BC1_cluster)
    	
    	with open(old_BC1_cluster) as f: # make a dictionary of old cluster id and old cluster sequence (key: id; value: sequence)
    		reader = csv.reader(f)
    		old_BC1_cluster = list(reader)
    	del old_BC1_cluster[0]
    	dict_old_cluster_01 = dict(zip([a for a,b,c,d in old_BC1_cluster], [b for a,b,c,d in old_BC1_cluster]))
    	
    	
    	with open(old_index_1) as f: ### Replace keys in the above dictionary with sequences
    		reader = csv.reader(f)
    		old_index_1 = list(reader)
    	del old_index_1[0]
    	dict_old_index_1 = dict(zip([dict_old_cluster_01[a] for a, b, c in old_index_1], [b for a,b,c in old_index_1])) # key: sequence, value: index of known-PPI-barcode
    	dict_old_mismatch_1 = dict(zip([dict_old_cluster_01[a] for a, b, c in old_index_1], [c for a,b,c in old_index_1])) # key: sequence, value: number of mismatches between cluster sequence and specific known-PPI-barcode
    	
    	#Take the overlap clusters and update the cluster indices to be the new ones
    	new_bartender_BC1_cluster_overlap = [i for i in new_bartender_BC1_cluster if i[1] in dict_old_index_1.keys()]
    	new_index_1 = zip([a for a, b, c, d in new_bartender_BC1_cluster_overlap],[dict_old_index_1[b] for a, b, c, d in new_bartender_BC1_cluster_overlap], [dict_old_mismatch_1[b] for a, b, c, d in new_bartender_BC1_cluster_overlap])
    	with open(directory + 'Updated_matched_index_01.csv', 'w') as file: # output the filtered cluster
    		file.writelines(','.join(["cluster_index", "known_barcode_index", "Distance"]) + '\n')
    		file.writelines(','.join([a[0], a[1], a[2]]) + '\n' for a in new_index_1)

    	#### extract the non_overlap clusters and follow the previous code to map to the known barcodes
    	unmatched_1 = [i for i in new_bartender_BC1_cluster if i[1] not in dict_old_index_1.keys()]
    	
    	#######################################################
    	#Do the same thing for the second cluster
    	with open(directory + bartender_BC2_cluster) as f:
        	reader = csv.reader(f)
        	new_bartender_BC2_cluster=list(reader)
        del new_bartender_BC2_cluster[0]
        new_bartender_BC2_cluster = [i for i in new_bartender_BC2_cluster if int(i[3]) > cluster_threshold]
        with open(directory + 'bartender_cluster_02_filtered.csv', 'w') as file: # output the filtered cluster
        	file.writelines(','.join([a[0], a[1]]) + '\n' for a in new_bartender_BC2_cluster)
    	
    	with open(old_BC2_cluster) as f: # make a dictionary of old cluster id and old cluster sequence (key: id; value: sequence)
    		reader = csv.reader(f)
    		old_BC2_cluster = list(reader)
    	del old_BC2_cluster[0]
    	dict_old_cluster_02 = dict(zip([a for a,b,c,d in old_BC2_cluster], [b for a,b,c,d in old_BC2_cluster]))
    	
    	
    	with open(old_index_2) as f: ### Replace keys in the above dictionary with sequences
    		reader = csv.reader(f)
    		old_index_2 = list(reader)
    	del old_index_2[0]
    	dict_old_index_2 = dict(zip([dict_old_cluster_02[a] for a, b, c in old_index_2], [b for a,b,c in old_index_2])) # key: sequence, value: index of known-PPI-barcode
    	dict_old_mismatch_2 = dict(zip([dict_old_cluster_02[a] for a, b, c in old_index_2], [c for a,b,c in old_index_2])) # key: sequence, value: number of mismatches between cluster sequence and specific known-PPI-barcode
    	
    	#Take the overlap clusters and update the cluster indices to be the new ones
    	new_bartender_BC2_cluster_overlap = [i for i in new_bartender_BC2_cluster if i[1] in dict_old_index_2.keys()]
    	new_index_2 = zip([a for a, b, c, d in new_bartender_BC2_cluster_overlap],[dict_old_index_2[b] for a, b, c, d in new_bartender_BC2_cluster_overlap], [dict_old_mismatch_2[b] for a, b, c, d in new_bartender_BC2_cluster_overlap])
    	with open(directory + 'Updated_matched_index_02.csv', 'w') as file: # output the filtered cluster
    		file.writelines(','.join(["cluster_index", "known_barcode_index", "Distance"]) + '\n')
    		file.writelines(','.join([a[0], a[1], a[2]]) + '\n' for a in new_index_2)
    		
    	#### extract the non_overlap clusters and follow the previous code to map to the known barcodes
    	unmatched_2 = [i for i in new_bartender_BC2_cluster if i[1] not in dict_old_index_2.keys()]
    	
    	with open(known_barcodes_1) as f:
			reader = csv.reader(f)
			protein_barcodes_1=list(reader)
    	with open(known_barcodes_2) as f:
        	reader = csv.reader(f)
        	protein_barcodes_2=list(reader)
    	
    	del protein_barcodes_1[0]
    	del protein_barcodes_2[0]
    	
    	index_match_1 = []
    	import time
    	print "First barcodes"
    	print time.ctime()
    	
    	# get the exactly matched indices 
    	unmatched_1_barcodes = [i[1] for i in unmatched_1]
    	protein_1_barcodes = [i[2] for i in protein_barcodes_1]
    	overlap_1 = intersection(unmatched_1_barcodes, protein_1_barcodes)
    	for i in range(len(overlap_1)):
    		unmatched_index = unmatched_1[unmatched_1_barcodes.index(list_to_str(overlap_1[i]))][0]
    		known_index = protein_barcodes_1[protein_1_barcodes.index(list_to_str(overlap_1[i]))][0]
    		index_match_1.append([unmatched_index, known_index, str(0)])
    	
    	unmatched_1_mismatch = [a for a in unmatched_1 if a[1] not in overlap_1]
		
  	# pairwise comparison between unmatched barcodes(first half) with known barcodes(first half)
    	
    	print time.ctime()
    	
        tracker =0
    	for i in range(len(unmatched_1_mismatch)):
    		for j in range(len(protein_barcodes_1)):
				a= unmatched_1_mismatch[i][1]
				b= protein_barcodes_1[j][2]
				mismatch = levenshtein_distance(a,b)
				if mismatch <= 2:
					tracker += 1
					index_match_1.append([unmatched_1_mismatch[i][0], str(j), str(mismatch)])
        	if i%100 ==0:
				print "Finish %d of %d iterations; %d match" %(i, len(unmatched_1_mismatch), tracker)
				print time.ctime()

       	with open(directory + 'Updated_unmatched_index_01.csv', 'w') as file:
       		file.writelines(','.join(["cluster_index", "known_barcode_index", "Distance"]) + '\n')
        	file.writelines(','.join(a) + '\n' for a in index_match_1)
        
        ##### Merge the updated matched index_01 with the updated unmatched index_01
        for i in range(len(new_index_1)):
        	index_match_1.append(list(new_index_1[i]))
        with open(directory + 'Unmatched_barcodes_protein_index_01.csv', 'w') as file:
       		file.writelines(','.join(["cluster_index", "known_barcode_index", "Distance"]) + '\n')
        	file.writelines(','.join(a) + '\n' for a in index_match_1)
       
       
		index_match_2 = []
    	print "second barcodes"
    	print time.ctime()
    	# first append the exactly matched indices
    	unmatched_2_barcodes = [i[1] for i in unmatched_2]
    	protein_2_barcodes = [i[2] for i in protein_barcodes_2]
    	overlap_2 = intersection(unmatched_2_barcodes, protein_2_barcodes)
    	for i in range(len(overlap_2)):
    		unmatched_index = unmatched_2[unmatched_2_barcodes.index(list_to_str(overlap_2[i]))][0]
    		known_index = protein_barcodes_2[protein_2_barcodes.index(list_to_str(overlap_2[i]))][0]
    		index_match_2.append([unmatched_index, known_index, str(0)])
    	unmatched_2_mismatch = [a for a in unmatched_2 if a[1] not in overlap_2]
    	 
    	 # pairwise comparison between unmatched barcodes(second half) with known barcodes(second half) 
    	print time.ctime()
    	tracker =0
    	for i in range(len(unmatched_2_mismatch)):
			for j in range(len(protein_barcodes_2)):
				a= unmatched_2_mismatch[i][1]
				b= protein_barcodes_2[j][2]
				mismatch = levenshtein_distance(a,b)
				if mismatch <= 2:
					tracker += 1
					index_match_2.append([unmatched_2_mismatch[i][0], str(j), str(mismatch)])
			if i%100 ==0:
				print "Finish %d of %d iterations; %d match" %(i, len(unmatched_2_mismatch), tracker)
				print time.ctime()
   
    	with open(directory + 'Updated_unmatched_index_02.csv', 'w') as file:
    		file.writelines(','.join(["cluster_index", "known_barcode_index", "Distance"]) + '\n')
        	file.writelines(','.join(a) + '\n' for a in index_match_2)
        
        for i in range(len(new_index_2)):
        	index_match_2.append(list(new_index_2[i]))
        with open(directory + 'Unmatched_barcodes_protein_index_02.csv', 'w') as file:
       		file.writelines(','.join(["cluster_index", "known_barcode_index", "Distance"]) + '\n')
        	file.writelines(','.join(a) + '\n' for a in index_match_2)
       
        

        unmatched_barcodes_combine(directory, known_barcodes_1, known_barcodes_2, unmatched_lintag1, unmatched_lintag2, unmatched_seqtag,  
                                   bartender_BC1_barcode, bartender_BC2_barcode, bartender_BC1_cluster, bartender_BC2_cluster, output)
    
    else:
		unmatched_barcodes_combine(directory, known_barcodes_1, known_barcodes_2, unmatched_lintag1, unmatched_lintag2, unmatched_seqtag, 
                                   bartender_BC1_barcode, bartender_BC2_barcode, bartender_BC1_cluster, bartender_BC2_cluster, output)

def intersection(list1, list2):
	return list(set(list1) & set(list2))  

def list_to_str(list): # one string
	return ''.join(list)                           

def levenshtein_distance(s1, s2):
  
    l1 = len(s1)
    l2 = len(s2)
 
    matrix = [range(l1 + 1)] * (l2 + 1)
    for zz in xrange(l2 + 1):
        matrix[zz] = range(zz, zz + l1 + 1)
    for zz in range(0, l2):
        for sz in range(0, l1):
            if s1[sz] == s2[zz]:
                matrix[zz + 1][sz + 1] = min(matrix[zz + 1][sz] + 1,
                                    matrix[zz][sz + 1] + 1, matrix[zz][sz])
            else:
                matrix[zz + 1][sz + 1] = min(matrix[zz + 1][sz] + 1,
                                    matrix[zz][sz + 1] + 1, matrix[zz][sz] + 1)
    return matrix[l2][l1]

def dictionary_index(input_index_file):
    import csv
    with open(input_index_file) as f:
        reader = csv.reader(f)
        index_1=list(reader)
    dict_match_1 = {}
    for i in index_1:
        unknown_index= i[0]
        known_index= i[1]
    
        if unknown_index in dict_match_1:
            dict_match_1[unknown_index].append(known_index)
        else:
            dict_match_1[unknown_index]={}
            dict_match_1[unknown_index]=[known_index]
    print "%d all possible matches" %len(dict_match_1.keys())
    duplicate_keys = 0
    duplicate_values = 0
    keys_remove = []
    for key, value in dict_match_1.iteritems():
        if len(value) > 1:
            duplicate_keys += 1
            duplicate_values += len(value)
            keys_remove.append(key)
    for key in keys_remove:
        del dict_match_1[key]
    print "Remove %d of unamtched barcodes which have been matched to multiple known barcodes" % duplicate_keys
    print "%d duplicated known barcodes" %duplicate_values
    print "%d unique matches" %len(dict_match_1.keys())
    
    return dict_match_1
    
        
def unmatched_barcodes_combine(directory,
                               known_barcodes_1, 
                               known_barcodes_2,
                               unmatched_lintag1,
                               unmatched_lintag2,
                               unmatched_seqtag,
                               bartender_BC1_barcode,
                               bartender_BC2_barcode,
                               bartender_BC1_cluster,
                               bartender_BC2_cluster,
                               output):
    import csv
    d1 = dictionary_index(directory + 'Unmatched_barcodes_protein_index_01.csv')
    d2 = dictionary_index(directory + 'Unmatched_barcodes_protein_index_02.csv')
    
    with open(known_barcodes_1, 'rb') as f:
        reader = csv.reader(f)
        b5= list(reader)
    with open(known_barcodes_2, 'rb') as f:
        reader = csv.reader(f)
        b6= list(reader)
#make dictionaries
    
    del b5[0]
    del b6[0]
    
    d5 = dict(zip([a for a,b,c in b5], [c for a,b,c in b5]))
    d6 = dict(zip([a for a,b,c in b6], [c for a,b,c in b6])) 
    #d7 = dict(zip([a for a,b,c in b5], [b for a,b,c in b5]))
    #d8 = dict(zip([a for a,b,c in b6], [b for a,b,c in b6])) 
    
    with open(directory + bartender_BC1_barcode, 'rb') as f:
    	reader = csv.reader(f)
    	c1 = list(reader)
	
	with open(directory + bartender_BC2_barcode, 'rb') as f:
		reader = csv.reader(f)
		c2 = list(reader)
	#with open(directory + bartender_BC1_cluster, 'rb') as f:
		#reader = csv.reader(f)
		#c3= list(reader)
	#with open(directory + bartender_BC2_cluster, 'rb') as f:
		#reader = csv.reader(f)
		#c4= list(reader)
    
    d9 = dict(zip([a for a,b,c in c1], [c for a,b,c in c1]))
    d10 = dict(zip([a for a,b,c in c2], [c for a,b,c in c2]))
    #d11 = dict(zip([a for a,b,c,d in c3], [b for a,b,c,d in c3]))
    #d12 = dict(zip([a for a,b,c,d in c4], [b for a,b,c,d in c4]))
    
    bc1_file = open(directory + unmatched_lintag1, "r")
    bc2_file = open(directory + unmatched_lintag2, "r")
    umi_file = open(directory + unmatched_seqtag, "r")
    
    d = dict(); e = dict(); #d is the count of each DBC, e is a list of UMIs for each DBC
    #PPI_barcodes=[]
    
    from itertools import izip
        
    for bc1, bc2, umi in izip(bc1_file, bc2_file, umi_file):
		bc1 = bc1.strip()
		bc2 = bc2.strip()
		umi = umi.strip()
		if bc1 in d9 and bc2 in d10:
			cluster_id_BC1 = d9[bc1]
			cluster_id_BC2 = d10[bc2]
			if cluster_id_BC1 in d1 and cluster_id_BC2 in d2:
				#reads = bc1 + '_' + bc2
				#cluster = d11[cluster_id_BC1] + '_' + d12[cluster_id_BC2]
				seq = d5[d1[cluster_id_BC1][0]] + '_' + d6[d2[cluster_id_BC2][0]]
				#PPI = d7[d1[cluster_id_BC1][0]] + '_' + d8[d2[cluster_id_BC2][0]]
				#PPI_barcodes.append([PPI, seq, cluster, reads])
				d[seq] = d.get(seq,0)+1
				j = e.get(seq)
				if j == None:
					e[seq] = [umi]
				else:
					j.append(umi)
					e[seq] = j
	
    #with open(directory + "Sequences_matched.csv", 'w') as file:
		#file.writelines(','.join(["PPI", "known_BC", "cluster_BC","Unmatched_BC"]) + '\n')
		#file.writelines(','.join(a) + '\n' for a in PPI_barcodes)
	
    m = list([['knownBC01_knownBC02', 'all_reads', 'deduped_reads']])
    for k, v in d.iteritems():
		temp = [k,str(v), str(len(set(e[k])))]
		m.append(temp)
	
    with open(directory + output, 'w') as file:
		file.writelines(','.join(a) + '\n' for a in m)	
		