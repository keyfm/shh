from Bio import SeqIO
import argparse,sys,re



''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Infer pseudogene from reference and query gene sequence (bot in fasta format).
    Fasta files need to contain only single sequence.
    Header: ['name','pseudo','refLenBP','refNumCodons','queryLenBP','queryNumCodons','propN','startCodon','prematureStop','frameshift','BaseSubRate','CodonSubRate']
    "# LEGEND: propN...proportion of N-calls in query sequence."
	"# LEGEND: start_codon: Y...first codon is a start codon; N...first codon contains N(s); O...first codon non-start codon (ATG|GTG|TTG)"
	"# LEGEND: prematurestop: R...reference faulty: no stop at all, first stop prior last codon, or seq not devisible of three; -...stop of query and ref at same codon; numeric...premature-stop-codon ('TAA|TGA|TAG') position; N...no stop-codon found in query (e.g. due to N's); E extension...query codon sequence stop later than reference codon stop"
	"# LEGEND: frameshift: F...frameshift present (not divisible by three); -...no framehsift existent"

                                ''',
                                epilog="Questions or comments? --> key@shh.mpg.de")
parser.add_argument("-q", "--query_fasta",dest="query_seq", help="Query fasta sequence ", type=argparse.FileType('rt'))
parser.add_argument("-r", "--reference_fasta", dest="ref_seq", help="Reference fasta sequence", type=argparse.FileType('rt'))
parser.add_argument('-e', "--extra_header", dest="header", action='store_true',default=False,help="Add header line to output")
parser.add_argument("-o", "--out", dest="out", help="Output table. If not specified stdout", type=argparse.FileType('w'),default=sys.stdout)
args = parser.parse_args()


## testing part input:
# qRec = SeqIO.read("/Users/fmk/Documents/shh/sese/tmp/PCL0050/00_3856/00_3856_consensus.66pA1T1.fasta", "fasta")
# rRec = SeqIO.read("/Users/fmk/Documents/shh/sese/zhou/preprint_warwick/WRAP85593%2020171102/Additional_Data_11-pan_genes/single_fa/PCL0050.fa", "fasta")
# qRec = SeqIO.read("/projects1/enterica/neolithic/tmp/mem/v1_ancBranchFull_etal/00_3856/gencons/00_3856_PCL0034_consensus.66pA1T1.fasta", "fasta")
# rRec = SeqIO.read("/projects1/enterica/neolithic/references.backup/custom/PCL_single_fa/PCL0034.fa", "fasta")

# q_dict = SeqIO.to_dict(SeqIO.parse("/Users/fmk/Documents/shh/sese/tmp/PCsxd
''' FUNCTIONS'''
def prop_N(record):
	return str(record.seq).count('N') / len(record.seq)

def start_codon(record):
	'''function reads fa seq and tests the first codon of this seq'''
	first_codon = str(record.seq[0:3])
	if re.match('(ATG|GTG|TTG)', first_codon):
		return "Y" # yes, statt codon present
	elif re.search('N', first_codon):
		return "N" # meaning it contains N; however, we do not test if possible non-N in the start codon make up any typically start codon (which should not matter). 
	else:
		return "O"

def stop_codon(qSeq,rSeq):
	'''function infers stop codon [TAA,TGA,TAG] position in query and reference. If no stop at all: NA '''
	qStop = [j for j, item in enumerate([qSeq[i : i + 3] for i in range(0, len(qSeq), 3)]) if re.search("TAA|TGA|TAG", item)] # get ordered list of stop codon indices
	if not qStop: # test if stop available at all: if not stop (aka empty list) we put NA
		qStop = ['NA']
	else:
		qStop[0] = qStop[0]+1
	rStop = [j for j, item in enumerate([rSeq[i : i + 3] for i in range(0, len(rSeq), 3)]) if re.search("TAA|TGA|TAG", item)]
	if not rStop:
		rStop = ['NA']
	else:
		rStop[0] = rStop[0]+1
	return [qStop[0],rStop[0]] # list of first stops, with 1-based positional index

def test_premature_stop(idxQR,rLen,rCodonLen):
	'''test q-stop before r-stop or no stop at all'''
	# sanity check ref: being divisible by three and w/o premature-stop
	if ((rLen % 3) > 0) or (idxQR[1] != rCodonLen) or idxQR[1] == "NA":
		return "R" # synonym for reference faulty. Not called as pseudo per se.
	else:
		if idxQR[0] == idxQR[1]:
			return "-" # when query and ref agree in stop-codon position
		elif idxQR[0] == "NA":
			return "N"
		elif idxQR[0] < idxQR[1]:
			return idxQR[0] # return position of stop codon ; extension not considered
		elif idxQR[0] > idxQR[1]:
			return "E"

def test_frameshift(qLen):
	'''test if query seq not devisible by 3'''
	if ((qLen % 3) > 0):
		return "F" # framehshift present
	else:
		return "-" # no frameshift

def estimate_substitution_rates(qSeq,rSeq):
	''' estimate base and codon substituion rate for all bases until first stop in query/ref; NOTE: N ignored; codon_substitution != AA change (not tested) '''
	rTrip = [rSeq[i : i + 3] for i in range(0, len(rSeq), 3)]
	qTrip = [qSeq[i : i + 3] for i in range(0, len(qSeq), 3)]
	tripMinCount = min(map(len, [qTrip, rTrip] )) # estimate rates only for shared length. Overhangs of ref/query do not contribute to codon / base substitution rates
	base_subst_count = 0
	codon_subst_count = 0
	for j, qCod in enumerate(qTrip[:tripMinCount]):
		codon_switch = False
		if qCod != rTrip[j]:
			# print(j,qCod,rTrip[j]) 
			for index,base in enumerate(qCod):
				if base != "N" and base != rTrip[j][index]:
					codon_switch = True
					base_subst_count = base_subst_count + 1
					# print(j,index,base,rTrip[j][index])
		if codon_switch == True:
			codon_subst_count = codon_subst_count + 1
	rate_qBase_subst = base_subst_count/len(qSeq)
	rate_qCodon_subst = codon_subst_count/len(qTrip)
	return [ "%.5f" % round(rate_qCodon_subst,5) , "%.5f" % round(rate_qBase_subst,5) ]

def is_pseudogene(premature_stop,frameshift):
	''' return "yes" when premature_stop, no stop at all, or frameshift detected '''
	if isinstance(premature_stop, (int, float, complex)) or (frameshift == "F"):
		return "yes"
	else:
		return "no"



'''MAIN'''
if __name__ == '__main__':
	qRec = SeqIO.read(args.query_seq, "fasta")
	rRec = SeqIO.read(args.ref_seq, "fasta")
	args.query_seq.close()
	args.ref_seq.close()
	## get stats
	qNprop = float("%.5f" % round( prop_N(qRec) ,5) )
	qLen = len(qRec.seq)
	rLen = len(rRec.seq)
	rCodonLen = len(rRec.seq)/3
	qCodonLen = len(qRec.seq)/3
	## infer start codon state
	startCod = start_codon(qRec)
	## first stop codon position in q and r (no stop > NA), and test for premature stop
	stopIdxQR = stop_codon(str(qRec.seq),str(rRec.seq))
	## prematurestop?
	premature_stop = test_premature_stop(stopIdxQR,rLen,rCodonLen)
	## frameshift?
	frameshift = test_frameshift(qLen)
	## codon/base substitution rates
	subRatesListCodBase = estimate_substitution_rates(str(qRec.seq),str(rRec.seq))
	## test if pseudogene
	pseudogene = is_pseudogene(premature_stop,frameshift)
	## write output table
	if args.header == True:
		args.out.write('\t'.join(['name','pseudo','refLenBP','refNumCodons','queryLenBP','queryNumCodons','propN','startCodon','prematureStop','frameshift','BaseSubRate','CodonSubRate']) + '\n')
	args.out.write('\t'.join( str(x) for x in [ rRec.name , pseudogene , rLen , rCodonLen , qLen , qCodonLen , qNprop, startCod , premature_stop , frameshift , subRatesListCodBase[1] , subRatesListCodBase[0] ]) + '\n')
	args.out.close()

exit()
