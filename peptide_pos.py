#!/usr/bin/env python
""" peptide_pos.py """
__author__ = "Matthias Hirsch-Hoffmann"
__copyright__ = "Copyright 2018, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Prototype"
#==========================================================================================
"""
script to read XLSX file 
xlsx file has a header
geneid in Col C
unmodified peptide sequence in Col E
modified peptide sequence in Col F
"""
#==========================================================================================
# define XLSX value columns
cols={}
cols['geneid']='C'
cols['peptide']='E'
cols['peptide_mod']='F'
#==========================================================================================
import sys,os
import click
import json
import openpyxl
#==========================================================================================
@click.command()
@click.argument('workbook', type=click.Path(exists=True,readable=True))
@click.argument('worksheet', type=click.STRING)
@click.argument('jsongff', type=click.Path(exists=True,readable=True))
def main(workbook,worksheet,jsongff):
	#basefilename 
	basewb=os.path.basename(workbook)
	outfile=open(basewb+'.txt', 'w')
	outfile.write("geneid\tstrand\tpeptide\tpeptide_mod\tpeptide_org\tmodification\tprotein position\tgenomic position\tstatus\n")
	#read json gff file
	proteins=read_jsongff(jsongff)
	uniquer={} #variable to keep all analysed entries to reducre runtime 
	#print "reading xlsx"
	W = openpyxl.load_workbook(workbook)
	""" 
	modify here if header row is missing
	"""
	for row in range(2,W[worksheet].max_row+1):
		geneid = W[worksheet][cols['geneid']+str(row)].value
		#next row if geneid is not annotated
		if geneid not in proteins:
			""" PRINT RESULT """
			print_result(outfile,geneid,[],'','','','',[],"not annotated ")
			continue #next row
		
		peptide         = W[worksheet][cols['peptide']+str(row)].value
		peptide_mod_org = W[worksheet][cols['peptide_mod']+str(row)].value
		protseq         = proteins[geneid]['protseq']
		#find unmodified peptide in protein-sequence
		peppos = protseq.lower().find(peptide.lower())
		if peppos == -1:
			""" PRINT RESULT """
			print_result(outfile,geneid,[],peptide,'',peptide_mod_org,'',[],"peptide could not be located in protein sequence")		
			#print_result(geneid,[],'','','','',[]," not annotated ")
			continue #next row
		"""
		create modified peptide where positions of modification
		are represented as .
		"""
		peptide_mod=peptide_cleanup(peptide_mod_org)
		"""
		check if geneid/peptide_mod already seen
		"""
		if geneid+"-"+peptide_mod in uniquer:
			""" PRINT RESULT """
			print_result(outfile,geneid,[],peptide,peptide_mod,peptide_mod_org,'',[],"already analyzed")
			#print "already analyzed: "+geneid+"/"+peptide+"="*60		
			continue #next row
		#add to uniquer
		uniquer[geneid+"-"+peptide_mod]=1
		"""
		do we have a modification?
		"""
		peptide_mod_pos=peptide_mod.find('.')		
		if peptide_mod_pos == -1:
			""" PRINT RESULT """
			print_result(outfile,geneid,[],peptide,peptide_mod,peptide_mod_org,'',[],"peptide does not have a modification of interrest")
			continue #next row
		
		#dna position for each amino acid
		aas = get_aas(geneid,proteins)
		#print aas
		mymodpos=1
		while peppos >=0:
			#print geneid+" ("+proteins[geneid]['strand']+")"+" hit at "+str(peppos)+" "+peptide_mod
			ppos=0 #position in the peptide
			for aapos in range(peppos,peppos+len(peptide)):
				if peptide_mod[ppos:ppos+1]=='.':
					if proteins[geneid]['strand']=='+' and aas[aapos][0]+1 == aas[aapos][1] and aas[aapos][0]+2 == aas[aapos][2]:
						#print peptide_mod[ppos:ppos+1],aas[aapos]
						""" PRINT RESULT """
						print_result(outfile,geneid,proteins,peptide,peptide_mod,peptide_mod_org,aapos,aas,"ok",mymodpos)
					elif proteins[geneid]['strand']=='-' and aas[aapos][0]-1 == aas[aapos][1] and aas[aapos][0]-2 == aas[aapos][2]:
						#print peptide_mod[ppos:ppos+1],aas[aapos]
						""" PRINT RESULT """
						print_result(outfile,geneid,proteins,peptide,peptide_mod,peptide_mod_org,aapos,aas,"ok",mymodpos)
					else:
						""" PRINT RESULT """
						print_result(outfile,geneid,[],peptide,peptide_mod,peptide_mod_org,aapos,[],"modification in splitted codon",mymodpos)
						#print peptide_mod[ppos:ppos+1],aas[aapos]
					mymodpos+=1
						
				ppos+=1 #increase peptide position
			peppos = protseq.lower().find(peptide.lower(),(peppos+1)) #check if peptide exists multiple times in protein
	outfile.close()
	return 1			
#==========================================================================================
def print_result(outfile,geneid,proteins,peptide,peptide_mod,peptide_mod_org,aapos,aas,msg,modpos=''):
	result=geneid
	if proteins:
		result+="\t"+proteins[geneid]['strand']
	else:
		result+="\t"
		
	result+="\t"+peptide+"\t"+peptide_mod+"\t"+peptide_mod_org
	result+="\t"+str(modpos)
	
	if not aapos=='':
		result+="\t"+str(aapos+1) #position of the modification in protein
	else:
		result+="\t"

	if aas:
		result+="\t"+str(aas[aapos][1]) # middle position of codon
	else:
		result+="\t"
		
	result+="\t"+msg	
	
	outfile.write(result+"\n")
#	print proteins[geneid]['protseq']
#	print proteins[geneid]['exons']
#	for i in range(0,len(proteins[geneid]['protseq'])):
#		print i,proteins[geneid]['protseq'][i:i+1],aas[i]
#	sys.exit(1)
#==========================================================================================
def peptide_cleanup(peptide):
	""" 
	replace all modifications by its original character
	""" 
	peptide_org=peptide
	peptide=peptide.replace('(t)','T').replace('(T)','T')
	peptide=peptide.replace('(s)','S').replace('(S)','S')
	peptide=peptide.replace('(y)','Y').replace('(Y)','Y')
	peptide=peptide.replace('(oxM)','M').replace('(m)','M')
	peptide=peptide.replace('(c)','C').replace('(C)','C').replace('(C*)','C').replace('(c*)','C')
	peptide=peptide.replace('(ac)','')
	peptide=peptide.replace('(K)','K')
	""" 
	replace all wanted modifications into . (dot)
	"""
	peptide=peptide.replace('(pS)','.').replace('(pT)','.').replace('(pY)','.')
	#print peptide
	if not peptide.find('(') == -1:
		sys.exit('*** ERROR *** unhandled modification in peptide '+peptide+" ("+peptide_org+") at position "+str(peptide.find('(')))
	return peptide
#==========================================================================================	
def get_aas(geneid,proteins):
	#create aa-genome-position array
	aas={}
	if geneid in proteins:
		#print proteins[geneid]
		p=0
		aa=0
		if proteins[geneid]['strand']=='+':
			for exon in proteins[geneid]['exons']:
				for i in range(int(exon[0]),int(exon[1])+1):
					if not aa in aas:
						aas[aa]={}
					aas[aa][p]=i
					p+=1
					if p==3:
						aa+=1
						p=0
		else:
			for exon in proteins[geneid]['exons']:
				for i in range(int(exon[1]),int(exon[0])-1,-1):
					if not aa in aas:
						aas[aa]={}
					aas[aa][p]=i
					p+=1
					if p==3:
						aa+=1
						p=0
	return aas
#==========================================================================================
#read json file with CDS information
def read_jsongff(jsongff):
	proteins={}
	#print "reading json gff"
	with open(jsongff) as data_file:    
		proteins = json.load(data_file)
	return proteins

#==========================================================================================
if __name__ == '__main__':
	sys.exit(main())


