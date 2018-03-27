import argparse,sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    Extract PCL ID "\t" geneName_productInfo
                                ''',
                                epilog="Questions or comments? --> key@shh.mpg.de")
parser.add_argument("-g", dest="genbank_file", help="Input genbank", type=argparse.FileType('rt'))
parser.add_argument("-o", "--out_fasta", dest="out", help="Output file", type=argparse.FileType('w'),default=sys.stdout)
args = parser.parse_args()


out = args.out



''' FUNCTIONS'''
def genbankHeader(seqIOrecord,outFile,feature='CDS'):
	"""
	Exrtact FA sequence for all desired features from genbank file
	"""
	for feat in seqIOrecord.features:
		if feat.type == feature:
			# turn Seqfeature to Seqrecord to use SeqIO.write function
			# tmp_record = SeqRecord(feat.extract(record.seq), id=feat.qualifiers['gene'][0]+'_'+feat.qualifiers['product'][0].title().replace(" ", ""), name=feat.qualifiers['gene'][0], description="")
			if 'pseudo' not in feat.qualifiers:
				header_id = str(feat.qualifiers.get('locus_tag')[0])+'\t'+feat.qualifiers.get('gene','NA')[0]+'_'+feat.qualifiers.get('product','NA')[0].title().replace(" ", "")+'_'+feat.qualifiers.get('protein_id','NA')[0]+"\n"
				out.write(header_id)
				# SeqIO.write(tmp_record, outFile, "fasta")

			

'''main'''

if __name__ == "__main__":
	record = SeqIO.read(args.genbank_file, "genbank")
	genbankHeader(record,out)
