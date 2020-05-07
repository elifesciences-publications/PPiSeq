### sort barcode strains together according to their PPIs
def sort_PPI_count(barcode_fitness, genotype_pos, barcode_pos, fitness_pos, distance_pos, G0_pos, output):
    import csv
    with open(barcode_fitness) as f:
        reader = csv.reader(f)
        d=list(reader)
    del(d[0]) # function if the input file has header
    ### change this according to the input
    genotype_dict= {}
    for i in d:
        genotype= i[genotype_pos - 1]
        barcode= i[barcode_pos - 1]
        fitness= i[fitness_pos - 1]
        distance= i[distance_pos - 1]
        G0= i[G0_pos -1 ]
        G6= i[G0_pos]
        G9= i[G0_pos + 1]
        G12= i[G0_pos + 2]
        
       
    
        if genotype in genotype_dict:
            genotype_dict[genotype]['barcode'].append(barcode)
            genotype_dict[genotype]['fitness'].append(fitness)
            genotype_dict[genotype]['distance'].append(distance)
            genotype_dict[genotype]['G0'].append(G0)
            genotype_dict[genotype]['G6'].append(G6)
            genotype_dict[genotype]['G9'].append(G9)
            genotype_dict[genotype]['G12'].append(G12)
        
        else:
            genotype_dict[genotype]={}
            genotype_dict[genotype]['barcode']=[barcode]
            genotype_dict[genotype]['fitness']= [fitness]
            genotype_dict[genotype]['distance']= [distance]
            genotype_dict[genotype]['G0']=[G0]
            genotype_dict[genotype]['G6']=[G6]
            genotype_dict[genotype]['G9']= [G9]
            genotype_dict[genotype]['G12']=[G12]
    
    import numpy as np
    PPI_name= []
    counts= []
    barcode= []
    fitness= []
    distance= []
    G0= []
    G6= []
    G9= []
    G12= []
    
    for i in genotype_dict.keys():
        kinds= len(genotype_dict[i]['barcode'])
        for j in range(kinds):
            PPI_name.append(i)
            counts.append(kinds)
            barcode.append(genotype_dict[i]['barcode'][j])
            fitness.append(genotype_dict[i]['fitness'][j])
            distance.append(genotype_dict[i]['distance'][j])
            G0.append(genotype_dict[i]['G0'][j])
            G6.append(genotype_dict[i]['G6'][j])
            G9.append(genotype_dict[i]['G9'][j])  
            G12.append(genotype_dict[i]['G12'][j])  
    
    fitness_individual=open(output, 'w')
    writer = csv.writer(fitness_individual)
    for PPI, count, barcode, fitness, distance, G0, G6, G9, G12 in zip(PPI_name, counts, barcode, fitness, distance, G0, G6, G9, G12):
        writer.writerow([PPI, count, barcode, fitness, distance, G0, G6, G9, G12])
    fitness_individual.close()
