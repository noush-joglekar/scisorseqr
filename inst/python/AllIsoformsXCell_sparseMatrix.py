### By Anoushka Joglekar
### Modified 2019_02_27

import sys
import pandas as pd
import time
from itertools import chain

start_time = time.time()

input_file = sys.argv[1]
all_info = [x.strip('\n').split('\t') for x in open(input_file).readlines()]
iso_names = [x[0] for x in all_info]
cellsPerIso = [x[1::2] for x in all_info]
numIsoPerCell = [[int(y) for y in x[2::2]] for x in all_info]

print("Creating Isoform X cell matrix")
print("Processing file with ",len(all_info)," entries")


cell_list = list(set(list(chain.from_iterable([x for x in cellsPerIso]))))
lc = len(cell_list)

def MakeListToAppend(index):
	toAppend = [0]*lc
	ixes = [cell_list.index(x) for x in cellsPerIso[i]]
	for ix in ixes:
		toAppend[ix] = numIsoPerCell[i][ixes.index(ix)]
	return(toAppend)

tA = []
for i in range(len(all_info)):
	tA.append(MakeListToAppend(i))
	if i%10000 == 0:
		print(i/1000,'k done')

print("Writing to dataframe")
numIso_frame = pd.DataFrame(tA,columns = cell_list, index = iso_names)

elapsed_time = time.time() - start_time
print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

print("Writing output file")
numIso_frame.to_csv('Matrices/AllIsoformsXAllCells_matrix.csv',sep="\t")
