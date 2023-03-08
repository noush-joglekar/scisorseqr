### By Anoushka Joglekar

import sys
import pandas as pd
import time

start_time = time.time()

print("Creating Isoform X cell matrix")

input_file = sys.argv[1]
all_info = [x.strip('\n').split('\t') for x in open(input_file).readlines()][1:]

print("Processing file with ", len(all_info), " entries")

gene_iso_names = [x[0] for x in all_info]
cluster_names = [x[1] for x in all_info]
num_iso = [int(x[2]) for x in all_info]

cluster_list = list(set(cluster_names))
iso_list = list(set(gene_iso_names))

# indexes = [gene_iso_names.index(iso) for iso in iso_list]

numIso_frame = pd.DataFrame(index=iso_list, columns=cluster_list)
numIso_frame = numIso_frame.fillna(0)

for i in range(len(all_info)):
    numIso_frame.loc[gene_iso_names[i], cluster_names[i]] += num_iso[i]

elapsed_time = time.time() - start_time
print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

print("Writing output file")
numIso_frame.to_csv('Matrices/IsoformXCluster_matrix.csv', sep="\t")
