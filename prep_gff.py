#!/usr/bin/env python
""" prep_gff.py """
__author__ = "Matthias Hirsch-Hoffmann"
__copyright__ = "Copyright 2018, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Prototype"

"""
This script needs 2 input files: 
1. gff annotation file
2. fasta file with sequences of above annotated genes

The script will create protein sequence based on CDS 
annotation and save CDS information, strand and protein-sequence
in the output as JSON for later usage in the script which
determines the position of modifications.

Depending on gff (no clear rules) you 
have to modify the script to correctly extract
protein names to match input file for peptide position determination
"""
import sys,os
import click
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
#==========================================================================================
def get_ioseq(genomefasta):
	#read fasta file and create dict with sequence
	ioseq={}
	for record in SeqIO.parse(genomefasta, 'fasta'):
		ioseq[record.id]=record.seq
	return ioseq
#==========================================================================================
def seqN(sequence):
	#fills sequence with N to be a multiple of 3
	mod = len(sequence)%3
	if mod==1:
		sequence+='NN'
	elif mod==2:
		sequence+='N'
	return sequence	
#==========================================================================================
@click.command()
@click.argument('gff', type=click.Path(exists=True,readable=True))
@click.argument('genomefasta', type=click.Path(exists=True,readable=True))
def main(gff,genomefasta):
	basegff=os.path.basename(gff)
	proteinexons={}
	ioseq=get_ioseq(genomefasta)
	"""
	read gff file and create dict for each protein.
	"""
	with open(gff) as f:
		for line in f:
			#split the line
			cols=line.replace('\n','').split('\t')
			if not cols[0][0:1]=='#' and cols[2]=='CDS':
				
				if cols[0] not in ioseq:
					#the chromosome is missing
					sys.exit('Sequence missing in fasta '+genomefasta+' for entry '+cols[0])
				
				#split the tags
				tags=cols[8].split(';') # may produce empty tag
				for tag in tags:
					#split the tag
					if len(tag) > 0:
						key,value=tag.split('=')
						if key=='Parent':
							""" check the content of value """
							#print key,value
							"""
							adjust here if protein name does not match 
							protein names in xlsx file
							"""
							parents=value.split(',')
							protname=parents[0]
							
							if protname not in proteinexons:
								proteinexons[protname]={}
								proteinexons[protname]['strand']=cols[6]
								proteinexons[protname]['seqio']=cols[0]
								proteinexons[protname]['exons']=[]
								proteinexons[protname]['CDS']=''
							#append CDS exon to list
							proteinexons[protname]['exons'].append([cols[3],cols[4]])
	"""
	the exons are sorted, pieces are cut out of the sequece of 
	the chromosome fasta file, joined, reverse-complemented if 
	needed and translated into protein sequence.
	"""
	for protein in proteinexons:
		#sort exons based on strand
		if proteinexons[protein]['strand']=='+':
			#sort asc on start pos for forward strand
			proteinexons[protein]['exons']=sorted(proteinexons[protein]['exons'], key=lambda x: x[0])
		else:
			#sort desc on end pos pro reverse strand
			proteinexons[protein]['exons']=sorted(proteinexons[protein]['exons'], key=lambda x: x[1], reverse=True)
		#build complete CDS Sequence
		for exons in proteinexons[protein]['exons']:
			#cut sequence out, python 0-index reduce start-pos by 1!!!
			seqpiece = ioseq[proteinexons[protein]['seqio']][(int(exons[0])-1):int(exons[1])]
			if proteinexons[protein]['strand']=='+':
				#append to then end for forward
				proteinexons[protein]['CDS']+=seqpiece
			else:
				#append before for reverse
				proteinexons[protein]['CDS']=seqpiece+proteinexons[protein]['CDS']
			
		coding_dna = Seq(str(proteinexons[protein]['CDS']), generic_dna)
		#check strand... create reverse if necessary
		if proteinexons[protein]['strand']=='-':
			#reverse-complement CDS
			coding_dna=coding_dna.reverse_complement()
		#fill sequence with N to complete last codon
		coding_dna = Seq(seqN(str(coding_dna)),generic_dna)
		#translate the sequence
		proteinexons[protein]['protseq'] = str(coding_dna.translate(stop_symbol=''))
		"""
		CDS, seqio are not needed anymore and key(s) are removed
		from dictionary to keep it as small as possible
		"""
		del proteinexons[protein]['CDS'] 
		del proteinexons[protein]['seqio'] 
		#print protein,proteinexons[protein]['strand']
		#print proteinexons[protein]['CDS']
		#print proteinexons[protein]['protseq']
		#print proteinexons[protein]
	
	with open(basegff+'.json', 'w') as outfile:
		#json.dump(proteinexons,outfile, indent=4, sort_keys=True, ensure_ascii=True)
		json.dump(proteinexons,outfile, indent=None, sort_keys=True, ensure_ascii=True)

	return 1

#==========================================================================================
if __name__ == '__main__':
	sys.exit(main())
		
