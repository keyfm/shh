

import argparse,sys,re
from Bio import SeqIO
from multiprocessing.dummy import Pool as ThreadPool
import numpy as np
import pandas as pd


''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Produces two matrices: (A) Hamming Distance Matrix, and (B) Count valid pairwise base pairs, on MultitFasta Alignemnt file (e.g. output multivcfanalyzer).
    If pairwise comparison contains N, site is ignored.
    Multithreading optional.
                                ''',
                                epilog="Questions or comments? --> key@shh.mpg.de")
parser.add_argument("-f", dest="multifa_file", help="Input MultitFasta. Each seq has to have same length.", type=argparse.FileType('rt'))
parser.add_argument('-t', dest="threads", action="store", type=int, default=1, help="Number of threads. Default: 1")
parser.add_argument("-g", "--genome", dest="specGenome", help="Specific genome ID. Output matrix only contains data for requested genome. Currently works only by supplying single genome. Default: All",action="store",default=["all"] )
parser.add_argument("-o", "--out", dest="outTag", help="outputPathFileNameTag",action="store")
args = parser.parse_args()



''' FUNCTIONS'''

# def hamdist(pairList, recordDict):
def hamdist(pairID):
	"""Count the # of differences between equal length strings str1 and str2. Return count and number of valid pairwise comparisons"""
	diffs = 0
	ctrBase = 0
	# record_dict global, thus can be used in function! Consider to replace re.search with explicit ch1/ch2 != N check. Might be faster.
	for ch1, ch2 in zip( record_dict[ pairID[0] ] , record_dict[ pairID[1] ] ):
		if not re.search('N', ''.join([ch1,ch2])):
			ctrBase += 1
			if ch1 != ch2:
				diffs += 1
	return [diffs , ctrBase]

def list2pairs(longlist):
	'''all id pairs we want. Note: Same ID pairs kept!'''
	halfpairs = []
	fullpairs = []
	for i in longlist:
		for j in longlist:
			fullpairs.append([i,j])
			if [j,i] not in halfpairs:
				halfpairs.append([i,j])
	return [halfpairs , fullpairs]

def results2separateList(ham_res_lol,list_pairs_half_full):
	'''Separate listoflist elements into separate list. Extend list when pair not calculated, due to redundancies'''
	sub = list()
	base = list()
	ctr = 0
	for fullPair in list_pairs_half_full[1]:
		if fullPair in list_pairs_half_full[0]:
			sub.append( ham_res_lol[ctr][0] )
			base.append( ham_res_lol[ctr][1] )
			ctr += 1
		else:
			sub.append(0)
			base.append(0)
	return [sub,base]



'''main'''

if __name__ == "__main__":
	# record_dict = SeqIO.to_dict(SeqIO.parse("/Users/fmk/Documents/shh/sese/mvcf/cap12/ana_v2_5x_c1tr4b_c2r1234_seseRefsParaC/snpAlignment.fasta", "fasta"))
	# record_dict = SeqIO.to_dict(SeqIO.parse("/projects1/enterica/neolithic/mvcf_ana.backup/cap12/ana_v1_3x_c1tr4b_c2r1234_seseRefsParaC/snpAlignment.fasta", "fasta"))
	record_dict = SeqIO.to_dict(SeqIO.parse( args.multifa_file , "fasta"))
	keys_multifa = record_dict.keys() # sve so we know order! even though py3.6+ saves ordered dictionary keys (just to keep code recyclable)
	list_pairs_half_full = list2pairs(keys_multifa)
	# modify halfpair list to represent only desired genome-pairs (if any requested)
	if not args.specGenome == 'all':
		list_pairs_half_full[0] = [ list_pairs_half_full[0][j] for j, item in enumerate(list_pairs_half_full[0]) if re.search(args.specGenome , item) ]
	# estraxt hamming distance and num valid base-pairs in parallel
	pool = ThreadPool(args.threads) 
	results = pool.map(hamdist, list_pairs_half_full[0])
	BaseCtrLs = results2separateList(results,list_pairs_half_full) # separates and adds 0 for the otherwise redundant pairs: i,j == X; j,i == 0
	# turn lists into np.array necessary for turning it into matrix, name header/row and write as tab separated file
	SubArray = np.array(BaseCtrLs[0])
	matrixSub = SubArray.reshape(len(keys_multifa),len(keys_multifa))
	dfSub = pd.DataFrame(matrixSub, index=keys_multifa, columns=keys_multifa)
	dfSub.to_csv(args.outTag + "_substitution_mtx.tsv", index=True, header=True, sep='\t')
	CtrArray = np.array(BaseCtrLs[1])
	matrixCtr = CtrArray.reshape(len(keys_multifa),len(keys_multifa))
	dfCtr = pd.DataFrame(matrixCtr, index=keys_multifa, columns=keys_multifa)
	dfCtr.to_csv( args.outTag + "_positionCount_mtx.tsv", index=True, header=True, sep='\t')

exit(0)