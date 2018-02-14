# PTM Scripts
Various script related to analyse/process PTM (post translational modification) proteomics data

## prep_gff.py
Script to preprocess a genome file with its GFF annotation file.
```
python prep_gff.py <gff-file> <genome-fasta-file>
```
Outputs a json file containing for each annotated protein its protein sequence, strand and all CDS feature positions. 
Filename of the new created json file will be <gff-file>.json. The Output file is needed as input file for script [**peptide_pos.py**](#peptide_pos.py)

### Prerequisites
* click
* json
* Bio
* Bio.Seq
* Bio.Alphabet

## peptide_pos.py
Script to read XLSX file containing geneid, unmodified peptide sequence and modified peptide sequence to determine position of PTM in protein sequence and on the genome. Further input is the name of the worksheet within XLSX and the gff-json-file create with script [**prep_gff.py**](#prep_gff.py).

Development was performed with a XLSX file having geneid in column C, unmodified peptide in column E and modified peptide in column F. Please adapt script accoring your needs (line 20 ff). 

Depending on the input modified peptide sequence the peptide\_cleanup()-function has to be adapted, e.g. TVNGG<b>(oxM)</b>S<b>(pS)</b>PSLQK vs. TVNGG**M147**S**S166**PSLQK (line 145 ff).

```
python peptide_pos.py <XLSX Workbook-file> <Worksheet> <gff-json-file>
``` 
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

### Prerequisites
* click
* json
* openpyxl

## License
Copyright &copy; 2018, Matthias Hirsch-Hoffmann, ETH Zurich. Licensed under the GNU General Public License, version 3 (https://www.gnu.org/licenses/gpl-3.0.en.html) or see [LICENSE.txt](LICENSE.txt).

## Disclamer
This software is supplied 'as is' without any warranty or guarantee of support. ETH Zurich is not responsible for its use, misuse, or functionality. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability arising from, out of, or in connection with this software.

