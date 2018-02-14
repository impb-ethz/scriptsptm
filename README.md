# scriptsptm
Various script related to analyse/process PTM (post translational modification) proteomics data

## pref_gff.py
Script to preprocess a genome file with its GFF annotation file.

python pref_gff.py <gff-file> <genome-fasta-file>

Outputs a json file containing for each annotated protein its protein sequence, strand and all CDS feature positions. 
Filename of the new created json file will be <gff-file>.json. The Output file is needed as input file for script 

### Prerequisites
* click
* json
* Bio
* Bio.Seq
* Bio.Alphabet

**peptide_pos.py**

## peptide_pos.py
Script to read XLSX file containing proteinid, unmodified peptide sequence and modified peptide_sequence to determine position of PTM in protein sequence and on the genome.

python peptide_pos.py <XLSX Workbook-file> <Worksheet> <gff-json-file>
  
Outputs a tab-delimited file containing 
* geneid
* strand
* unmodified peptide
* modified peptide (all PTM's of interest are replaced by a dot (.))
* original modified peptide 
* modification number
* position of modification on protein
* position of modification on genome
* status






