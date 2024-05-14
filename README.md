# Pymol-pdb index aligner
This repository is a series of python scripts to clean, fetch sequences, align, and report adjusted indexes given a .csv of GUs of interest.

#more specifically...
When pulling sequence information from pymol and its pdb accession, often times things don't match up.

Pymol is often missing residues that could not be resolved, but it indexes match the measurements you make
within the program. When extracting the pymol sequence using XXXX, the sequence will be reformatted in fasta
format with Ns such that the length matches the pymol indexes.

When fetching the sequence from a pdb accession, the fasta file is often more complete than in pymol. However,
the first nucleotide doesn't always correspond to index 1 of the pymol file. After the pymol and pdb fastas are
both fetched, they can be aligned with XXXX that uses a needle aligned to calculate the "pdb_adjusted_index" from
the pymol index.

We also ran into trouble that not all molecules are labelled with the same index. For this reason, the XXX script
will align pdb fastas from the same organism for a final "adjusted index" so you can directly compare indexes within
a species.

An example file is included to compare XXXX INFORMATION. However, these scripts are compatible with an input csv with
any number of measurement coumns, as long as it includes X,Y,Z.



##SCRIPT 1: 0_molecule_format.py
This script was written to refine molecule names for more consistency. Users can go in an replace labels as they see fit.

#example to run
python 0_molecule_format.py input.csv




##SCRIPT 2: 1_get_fastas.py
This script fetches pdb, pymol, or reference sequences necessary for the following alignment step and saves them
in a fasta format. PDB sequences are from each entry's pdb access. Pymol sequences are extracted and reformatted for each entry using pymol. Reference sequences include only one example of each molecule type per organism.
*** this script will report errors for the pdb sequence that might need to be filled in by hand by downloading the full fasta and finding the sequence of interest by eye**

#dependencies:
pymol

#how to run
##for Pymol fasta
python 1_get_fastas.py -i input.csv -ft pymol -o pymol_out_name.fasta
##for pdb fasta
python 1_get_fastas.py -i input.csv -ft pdb -o pdb_out_name.fasta
##for reference fasta
pyton 1_get_fastas.py -i input.csv -ft ref -pdb pdb.fasta -o ref_out_name.fasta





##SCRIPT 3: 2_adjust_index.py
This script takes each of the fasta files from the last script and returns a csv file that has all the original information as the input csv along with pdb_adjusted_indexes and ref_adjusted_indexes

#dependencies
BioPython

#how to run
python 2_adjust_index.py -i input.csv -py pymol.fasta -pdb pdb.fasta -ref reference.fasata -o output_name.csv
