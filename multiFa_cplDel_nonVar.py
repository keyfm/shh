#!/usr/bin/env python3

###USAGE: python3 completedeletioncode.py <fasta_filename> <percentage_of_deletion> <output_filename> <0 for variant sites only; anyother number above 0 for including invariant sites> '####
'''
Script filters sites that contain N in any seq in multifasta
4th argument set to 0, will filter in addition for variable sites only (i.e. >1 non-N character!)
'''

import argparse,sys

''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Filters MultiFasta for CompleteDeletion (i.e. removes sites w/ "N"), and for invariable sites (if specified)
    NOTE: Seq order not maintained.
                                ''',
                                epilog="Questions or comments to Adytia ;) ")
parser.add_argument("-f", dest="fasta_file", help="Input multifasta", type=argparse.FileType('rt'))
parser.add_argument("-o", "--out_fasta", dest="out", help="Output multifasta", type=argparse.FileType('w'),default=sys.stdout)
parser.add_argument("-l", "--seq_length", dest="seq_length", help="Number of bases to be printed per line in output fasta [100]", type=int,default=100)
parser.add_argument('-n', dest="allSites", help="Set flag to print non-varaible sites too [False]", action='store_true', default=False) # default= not necessary as implied by action
args = parser.parse_args()


out = args.out


## read the file and splits it and make a hash table####
with args.fasta_file as f:
    g=f.read().strip().split('>')

Alignment={}
for i in g[1:]:
    Alignment[i.split('\n')[0]]=[''.join(i.split('\n')[1:])]

### Calculate the length of Alignment ###
length=[]
for i in Alignment:
    length.append(len(Alignment[str(i)][0][0:]))
###select for seqeunces within the percentage of deletion ####
p=list(set(length))
nucbases=['A','T','G','C']

## build alignment matrix (lol) and filter for N (i.e. complete deletion)
m=[]
for j in range(0,p[0]):
    s=[]
    for i in Alignment:
        s.append(Alignment[i][0][j])
    if 'N' in s:
        pass
    else:
        m.append(s)

## filter for variable sites
m1=[]
if args.allSites==False:
    for o in m:
        if (len(list(set(o))))>1:
            m1.append(o)
else:
    m1=m

del m # remove m to save resources (in case alignment is large!)

## re-build fasta
Alignment1={}
a=0
for k in Alignment:
    n=[]
    for l in m1:
        n.append(l[a])
    Alignment1[k]=[''.join(n)]
    a=a+1
### Write to the file ###
seq_length = args.seq_length
o = args.out
for i in Alignment1:
    # o.write('>'+str(i)+'\n'+Alignment1[str(i)][0]+'\n')
    o.write('>'+str(i)+'\n')
    sequence = Alignment1[str(i)][0]
    while len(sequence) > 0:
        o.write(sequence[:seq_length]+'\n')
        sequence = sequence[seq_length:]

o.close()
### END #####





