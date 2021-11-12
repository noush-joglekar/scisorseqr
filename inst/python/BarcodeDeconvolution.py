## By Anoushka Joglekar
## Created 2018. Edited AS 2019. Edited AJ 2020

import sys
import argparse
import gzip
import pandas as pd
import os
import multiprocessing as mp
import re
from traceback import print_exc

def revComp(my_seq):            ## obtain reverse complement of a sequence
        base_comp = {'A':'T', 'C':'G','G':'C', 'T':'A', 'N':'N', " ":" "}
        lms = list(my_seq[:-1]) ## parse string into list of components
        lms.reverse()
        try:
                lms = [base_comp[base] for base in lms]
        except TypeError:
                pass
        lms = ''.join(lms)      ## get string from reversed and complemented list
        return(lms)

def BC_scan(barcodes,start_scan,end_scan,seq):  ## Given start and end position of scanning sequence, look for barcodes
        str_list = [seq[i:i+16] for i in range(start_scan,end_scan)]    ## barcode file has "\n" in encoding
        bc_intersect = list(set(str_list).intersection(barcodes))
        if bc_intersect:
                bc_intersect = [element for element in bc_intersect]    ## allow for possibility of multiple barcodes matching
                bc_pos = [seq.index(element) for element in bc_intersect]
        else:
                bc_intersect = 'no_bc'
                bc_pos='-'
        return{'bc_pos':bc_pos,'bc_intersect':bc_intersect}



def TDetectFwd(seq):    ## Looks for T9 in the forward strand
        try:
                fwd_ix = seq.index('TTTTTTTTT')         ### Only detects the first occurrence of T9
                return(fwd_ix)
        except ValueError:
                fwd_ix = -1
                return(fwd_ix)

def TDetectRev(rev_seq):        ## If it was actually the reverse strand then it looks at T9 in the reverse complement
        try:
                rev_ix = rev_seq.index('TTTTTTTTT')
                return(rev_ix)
        except ValueError:
                rev_ix = -1
                return(rev_ix)


def TSODetectFwd(end_seq):      ## Looks for middle 15 bases of TSO in the forward strand
        tso = 'TGGTATCAACGCAGA'
        try:
                fwdTSO_ix = len(end_seq) - end_seq.index(revComp(tso))
                return(fwdTSO_ix)
        except ValueError:
                fwdTSO_ix = -1
                return(fwdTSO_ix)

def TSODetectRev(revEnd_seq):   ## If it was actually the reverse strand then it looks at TSO in the reverse complement
        tso = 'AAGCAGTGGTATCAACGCAGAGTACAT'
        try:
                revTSO_ix = len(revEnd_seq) - revEnd_seq.index(revComp(tso))
                return(revTSO_ix)
        except ValueError:
                revTSO_ix = -1
                return(revTSO_ix)

def prelim(args):
	global barcodes
	global cluster_id
	global umiLength

	file_name = re.split('/|.fq.gz|.fastq.gz',args.fq)[-2]
	print(file_name)

	# If output file already exists, delete it.
	if os.path.isfile('%sPolyT_BCDetection_%s_.csv' %(args.outDir, file_name)):
		os.system('rm %sPolyT_BCDetection_%s_.csv' %(args.outDir,file_name))

        # If secondary output file for reads that are too short already exists, delete it.
	if os.path.isfile('%sTooShort.csv' %(args.outDir)):
		os.system('rm %sTooShort.csv' %(args.outDir))

	bc_file = args.bcClust.replace(u'\xa0', u'')

	barcodes = [x.strip('\n').split('\t')[0] for x in open(bc_file).readlines()]
	cluster_id = [x.strip('\n').split('\t')[1] for x in open(bc_file).readlines()]
	
	if args.chemistry == "v2":
	  umiLength = 10
	elif args.chemistry == "v3":
	  umiLength = 12
	  
	return()

def addToDict(d, line, rn):
	seq = line[:-1][:200]
	rev_seq = revComp(line)[:200]
	end_seq = line[-200:-1]
	revEnd_seq = revComp(line)[-200:]

	d['Readname'].append(rn)
	d['length'].append(len(line))

	## Getting TSO stats for read
	fwdTSO_ix = TSODetectFwd(end_seq)
	revTSO_ix = TSODetectRev(revEnd_seq)

	if fwdTSO_ix == revTSO_ix == -1:
		d['tso_status'].append('-')
		d['tso_position'].append('-')
	elif fwdTSO_ix == -1 and revTSO_ix != -1:
		d['tso_status'].append('TSO_found')
		d['tso_position'].append([revTSO_ix])
	elif fwdTSO_ix != -1 and revTSO_ix == -1:
		d['tso_status'].append('TSO_found')
		d['tso_position'].append([fwdTSO_ix])
	elif fwdTSO_ix != -1 and revTSO_ix != -1:
		d['tso_status'].append('DoubleTSO')
		d['tso_position'].append('-')

	## Getting T9 position
	fwd_ix = TDetectFwd(seq)
	rev_ix = TDetectRev(rev_seq)

	if fwd_ix == rev_ix == -1:              ## If valuError, output is -1, so no polyT found
		d['T9_status'].append('poly_T_not_found')
		d['position'].append('-')
		d['bc_position'].append('-')
		d['BarcodeFound'].append('no_bc')
		d['Cluster'].append('no_clust')
		d['Strand_info'].append('none')
		d['UMIs'].append('-')

	elif fwd_ix == -1 and rev_ix != -1:     ## PolyT found in reverse complement only
		d['T9_status'].append('poly_T_found')
		d['position'].append([rev_ix])
		d['Strand_info'].append('rev')
		start_scan = rev_ix-36
		end_scan = rev_ix-6
		if start_scan >= 0 and end_scan > 0:
			bc_found = BC_scan(barcodes,start_scan,end_scan,rev_seq)
		elif start_scan < 0 and end_scan >0:
			start_scan = 0
			bc_found = BC_scan(barcodes,start_scan,end_scan,rev_seq)
		elif start_scan <0 and end_scan <= 0:
			bc_found = {'bc_pos':'-','bc_intersect':'no_bc'}
		if bc_found:
			d['BarcodeFound'].append(bc_found.get('bc_intersect'))
			d['bc_position'].append(bc_found.get('bc_pos'))
			if 'no_bc' in bc_found.get('bc_intersect'):
				d['Cluster'].append('no_clust')
				d['UMIs'].append('-')
			else:
				d['Cluster'].append([cluster_id[x] for x in [barcodes.index(item) for item in bc_found.get('bc_intersect')]])
				UMI_start = int(bc_found.get('bc_pos')[0])+16
				UMI_end = int(bc_found.get('bc_pos')[0])+16+umiLength
				d['UMIs'].append(rev_seq[UMI_start:UMI_end])
				#bc_count += 1
		else:
			d['BarcodeFound'].append('no_bc')
			d['Cluster'].append('no_clust')
			d['bc_position'].append('-')
			d['UMIs'].append('-')

	elif fwd_ix != -1 and rev_ix == -1:     ## PolyT found in sequence but NOT the reverse complement
		d['T9_status'].append('poly_T_found')
		d['position'].append([fwd_ix])
		d['Strand_info'].append('fwd')
		start_scan = fwd_ix-36
		end_scan = fwd_ix-16
		if start_scan >= 0 and end_scan > 0:
			bc_found = BC_scan(barcodes,start_scan,end_scan,seq)
		elif start_scan < 0 and end_scan >0:
			start_scan = 0
			bc_found = BC_scan(barcodes,start_scan,end_scan,seq)
		elif start_scan <0 and end_scan <= 0:
			bc_found = {'bc_pos':'-','bc_intersect':'no_bc'}
		else:
			print("wtf",fwd_ix,rev_ix,start_scan,end_scan)
		if bc_found:
			d['BarcodeFound'].append(bc_found.get('bc_intersect'))
			d['bc_position'].append(bc_found.get('bc_pos'))
			if 'no_bc' in bc_found.get('bc_intersect'):
				d['Cluster'].append('no_clust')
				d['UMIs'].append('-')
			else:
				d['Cluster'].append([cluster_id[x] for x in [barcodes.index(item) for item in bc_found.get('bc_intersect')]])
				UMI_start = int(bc_found.get('bc_pos')[0])+16
				UMI_end = int(bc_found.get('bc_pos')[0])+16+umiLength
				d['UMIs'].append(seq[UMI_start:UMI_end])
				#bc_count += 1
		else:
			d['BarcodeFound'].append('no_bc')
			d['Cluster'].append('no_clust')
			d['bc_position'].append('-')
			d['UMIs'].append('-')

	else:                   ## PolyT found in both. Could mean one of three things
		start_scan_f = fwd_ix-36
		end_scan_f = fwd_ix-16
		if start_scan_f >= 0 and end_scan_f > 0:
			bc_found_f = BC_scan(barcodes,start_scan_f,end_scan_f,seq)
		elif start_scan_f < 0 and end_scan_f >0:
			start_scan_f = 0
			bc_found_f = BC_scan(barcodes,start_scan_f,end_scan_f,seq)
		elif start_scan_f <0 and end_scan_f <= 0:
			bc_found_f = {'bc_pos':'-','bc_intersect':'no_bc'}

		start_scan_r = rev_ix-36
		end_scan_r = rev_ix-6
		if start_scan_r >= 0 and end_scan_r > 0:
			bc_found_r = BC_scan(barcodes,start_scan_r,end_scan_r,rev_seq)
		elif start_scan_r < 0 and end_scan_r >0:
			start_scan_r = 0
			bc_found_r = BC_scan(barcodes,start_scan_r,end_scan_r,rev_seq)
		elif start_scan_r <0 and end_scan_r <= 0:
			bc_found_r = {'bc_pos':'-','bc_intersect':'no_bc'}

		if bc_found_f and bc_found_r and bc_found_r.get('bc_intersect') != 'no_bc' and bc_found_f.get('bc_intersect') != 'no_bc':
			## BC found in forward AND reverse strand implies something is wrong, discard the read
			d['T9_status'].append('poly_T_found')
			d['position'].append([fwd_ix,rev_ix])
			d['BarcodeFound'].append('DoubleBC')
			d['bc_position'].append('-')
			d['Cluster'].append('no_clust')
			d['Strand_info'].append('both')
			d['UMIs'].append('-')
			#chimera += 1
		elif bc_found_f and not bc_found_r:     ## Barcode found in fwd strand, reverse strand polyT was a false positive
			d['T9_status'].append('poly_T_found')
			d['position'].append([fwd_ix])
			d['Strand_info'].append('fwd')
			d['BarcodeFound'].append(bc_found_f.get('bc_intersect'))
			d['bc_position'].append(bc_found_f.get('bc_pos'))
			if 'no_bc' in bc_found.get('bc_intersect'):
				d['Cluster'].append('no_clust')
				d['UMIs'].append('-')
			else:
				d['Cluster'].append([cluster_id[x] for x in [barcodes.index(item) for item in bc_found_f.get('bc_intersect')]])
				UMI_start = int(bc_found.get('bc_pos')[0])+16
				UMI_end = int(bc_found.get('bc_pos')[0])+16+umiLength
				d['UMIs'].append(seq[UMI_start:UMI_end])
				#bc_count += 1
		elif bc_found_r and not bc_found_f:     ## Barcode found in reverse strand, fwd T9 was a false positive
			d['T9_status'].append('poly_T_found')
			d['position'].append([rev_ix])
			d['Strand_info'].append('rev')
			d['BarcodeFound'].append(bc_found_f.get('bc_intersect'))
			d['bc_position'].append(bc_found_f.get('bc_pos'))
			if 'no_bc' in bc_found.get('bc_intersect'):
				d['Cluster'].append('no_clust')
				d['UMIs'].append('-')
			else:
				d['Cluster'].append([cluster_id[x] for x in [barcodes.index(item) for item in bc_found_f.get('bc_intersect')]])
				UMI_start = int(bc_found.get('bc_pos')[0])+16
				UMI_end = int(bc_found.get('bc_pos')[0])+16+umiLength
				d['UMIs'].append(rev_seq[UMI_start:UMI_end])
				#bc_count += 1
		else:
			d['T9_status'].append('poly_T_found')
			d['position'].append('-')
			d['Strand_info'].append('both')
			d['bc_position'].append('-')
			d['BarcodeFound'].append('no_bc')
			d['Cluster'].append('no_clust')
			d['UMIs'].append('-')

	return(d)


def chunkAndProcess(args,d,tooShort):
	step = 4	## read lines in file with a step size of 4
	startLine = 0
	sim_processes = args.numProc
	print("Number of simultaneous processes: ", sim_processes)

	readCount = sum(1 for i in gzip.open(args.fq, 'rb')) // 4
#	print("Number of reads in fastq file: ", readCount)

	toFork = (readCount // sim_processes) + 1
#	print("Number of reads per process: ", toFork) ## test

        # Set childNo to 1 to give first child that childNo.
	childNo = 1

	while readCount >= 1:
		isChild = os.fork()
		if isChild == 0:
			break
		else:
			if childNo == sim_processes:
				break
			else:
				startLine = startLine + (toFork * 4)
				readCount = readCount - toFork
				childNo += 1

	if isChild != 0:
		os.waitpid(isChild, 0)
		sys.exit()

	with gzip.open(args.fq,"rt",encoding ='utf-8') as f:
		for _ in zip(range(startLine), f):
			pass
		for lineno, line in enumerate(f, start = startLine):
			if lineno < startLine + (toFork * 4):
				if lineno % step == 0:
					rn = line[:-1].split(' ')[0]
				if lineno % step == 1 and len(line) >= 201:
					d = addToDict(d, line, rn)
				if  lineno % step == 1 and len(line) < 201:
                                        tooShort.append(rn)
			else:
				break
	writeTooShort(args, tooShort, childNo)
	df = pd.DataFrame(data=d)
	df.to_csv('%stmp%d' %(args.outDir, childNo), sep = "\t",index=False, header=False)

	return()


def writeTooShort(args, tooShort, childNo):
	tmp_d = {'Name':tooShort}
	tmp_df = pd.DataFrame(data=tmp_d)
	if childNo == 1:
		tmp_df.to_csv('%stoo_short%d' %(args.outDir, childNo), sep = "\t", index=False)
	elif childNo > 1:
		tmp_df.to_csv('%stoo_short%d' %(args.outDir, childNo), sep = "\t", index=False, header=False)

	# Append tmp file to output file. Appending (>>) is necessary because we
	# cannot know when each process will write to the output file.

	os.system('cat %stoo_short%d | tr -d "[ ]\'" >> %sTooShort.csv' %(args.outDir, childNo, args.outDir))
	#print("Writing too short %d data to output file" %(childNo)) #test

	# Remove unnecessary tmp file because it has already been written to the output file.
	if os.path.isfile("%stoo_short%d" %(args.outDir, childNo)):
	#print("Removing too_short%d" %(childNo)) #test
		os.system('rm %stoo_short%d' %(args.outDir, childNo))
	return()



def parse_args():
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
        parser.add_argument('fq', type=str, help='fastq.gz file')
        parser.add_argument('bcClust', type=str, help='single cell barcode-cluster assignment')
        parser.add_argument('--outDir', default="./", type=str, help='directory to put output in')
        parser.add_argument('--numProc', default=10,type=int, help='number of simultaneous processes')
        parser.add_argument('--chemistry', default="v2", type=str, help="10x chemistry version - v2 or v3")
        args = parser.parse_args()
        return args

def main():
	d = {'Readname':[],'T9_status':[],'Strand_info':[],'position':[],'BarcodeFound':[],'bc_position':[],\
'Cluster':[],'UMIs':[],'tso_status':[],'tso_position':[],'length':[]}
	tooShort = []
	args = parse_args()
	prelim(args)
	chunkAndProcess(args,d,tooShort)


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
