#this script takes the 3 fasta files from 1_get_fastas.py as well as the original csv 
#it adjust the residue index to one common fasta (reference.fasta) within species for rRNAs
#it's written to easily be combined with 1_get_fastas.py but since they are installed weird on my computer I kept them separate

#the general flow of the steps are:
#1. align pymol fasta to pdb fasta to fill in the gaps (this can be later modified downstream to get flanking sequences) and get pdb adjusted index
#2. align pdb fasta to reference fasta to get reference adjusted index
#3. output a .csv file with pdb, chain, res1 index, res2 index, ref pdb_chain, res1 ref adjusted index and res2 ref ajusted index

import argparse
import pandas as pd
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Application import ApplicationError
#from Bio.Align.Applications import NeedleCommandline


#set arguments
parser = argparse.ArgumentParser(description='This script will make three fasta files needed to align and remove redundancy')
parser.add_argument('-i', '--inputcsv', help='non-redundant csv file', required=True)
parser.add_argument('-py','--pymolfasta', help='fasta file from pymol', required=True)
parser.add_argument('-pdb','--pdbfasta', help='fasta file from pdb', required=True)
parser.add_argument('-ref','--referencefasta', help='fasta file from reference', required=True)
parser.add_argument('-o','--output', help='output fasta file', required=True)
args= parser.parse_args()

#read in csv file
df = pd.read_csv(args.inputcsv)

#define function to read in fasta file as dictionary
def fasta_to_dict(fasta_file):
    fasta_dict = {}
    header = None
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line
                fasta_dict[header] = ""
            elif header != None:
                fasta_dict[header] += line
    return fasta_dict

pymol_dict = fasta_to_dict(args.pymolfasta)
pdb_dict = fasta_to_dict(args.pdbfasta)
reference_dict = fasta_to_dict(args.referencefasta)



#define function to run needle alignment and reformat output to be useable
#def run_alignment(fasta1,fasta2):

#run needle command
def run_needle(asequence, bsequence):
    try:
        needle_cline = NeedleCommandline(
            asequence="asis::" + asequence,
            bsequence="asis::" + bsequence,
            gapopen=10,
            gapextend=0.5,
            outfile="align.txt")
        stdout, stderr = needle_cline()
        with open("align.txt", "r") as alignment_file:
            alignment_text = alignment_file.read()
        return alignment_text
    except ApplicationError as e:
        return "alignment failed: " + str(e)

#reformat alignment output
def alignment_format(alignment_file):
    #read in alignment output
    with open(alignment_file, "r") as f:
        alignment_text = f.read()
    align_lines = []
    asequence = []
    bsequence = []
    #split alignment output into lines
    for line in alignment_text.splitlines():
        #if line starts with "asis"
        if line.startswith("asis"):
            #split line by space
            line = line.split()
            #concatonate 3rd element of line to asequence
            align_lines.append(line[2])
    #store lines in appropriate list
    for i in align_lines:
        if len(asequence) == len(bsequence):
            asequence.append(i)
        elif len(asequence) > len(bsequence):
            bsequence.append(i)
    #reformat asequence/bsequence into one giant string
    asequence = "".join(asequence)
    bsequence = "".join(bsequence)
    #make a pandas dataframe with one column for asequence, one column for b sequence, and one row for each character
    seq_align = pd.DataFrame({"asequence":list(asequence), "bsequence":list(bsequence)})
    #make a list starting at 1 with the first not dash in asequence, and increasing by 1 for each not dash in asequence, add - if asequence is a dash using pandas
    def counter_add(row, column_name):
        nonlocal counter
        if row[column_name] != "-":
            counter += 1
            return counter
        else:
            return counter
    counter = 0
    seq_align["asequence_index"] = seq_align.apply(counter_add, args=("asequence",), axis=1)
    counter = 0
    seq_align["bsequence_index"] = seq_align.apply(counter_add, args=("bsequence",), axis=1)
    return(seq_align)



#this function is getting replaced 
def mini_adjust_index(combined_list,index):
    for i in combined_list:
        if i[1] == index:
            return i[3]



def loop_flank2(res_info, pymoldict, pdbdict, refdict):
    #store reference molecule and organism as ref_org
    res_info["ref_org"] = ">" + res_info["Molecule_res1"] + "_" + res_info["Source_Organism_chain_res1"]
    #make a list of unique ref_org
    #make a list of unique pdb_id and chain (res1_data_ext)
    unique_pdb_chain = res_info["res1_data_ext"].unique()
    #align pymoldict and pdbdict with run_needle using pandas
    for i in unique_pdb_chain:
        #print progress with pdb_id and chain_id
        print(f"progress {unique_pdb_chain.tolist().index(i)+1}/{len(unique_pdb_chain)} {i}")
        #get pdb_id and chain_id
        pdb_id, chain_id = i.split("_")
        #get pymol sequence
        head1 = ">"+pdb_id + "_" + chain_id
        pymolsequence = pymoldict[head1]
        #get pdb sequence
        pdbsequence = pdbdict[head1]
        #run needle
        run_needle(pymolsequence,pdbsequence)
        #reformat alignment output
        align_df = alignment_format("align.txt")
        #get pdb adjusted index for each res_index_res1 and res_index_res2 where res1_data_ext == i
        for index, row in res_info.iterrows():
            if row["res1_data_ext"] == i:
        # Find matching indices in align_df
                matched_res1 = align_df[align_df["asequence_index"] == row["res_index_res1"]]["bsequence_index"].values
                matched_res2 = align_df[align_df["asequence_index"] == row["res_index_res2"]]["bsequence_index"].values
        # Check if matching indices are found
            if len(matched_res1) > 0:
                res_info.at[index, "pdb_adj_res1"] = matched_res1[0]
                res_info.at[index, "flank1"] = pdbsequence[int(matched_res1[0])-6:int(matched_res1[0])+5]
            else:
                res_info.at[index, "pdb_adj_res1"] = "Index Error: Matching index not found"
                res_info.at[index, "flank1"] = "Index Error: Matching index not found"
            if len(matched_res2) > 0:
                res_info.at[index, "pdb_adj_res2"] = matched_res2[0]
                res_info.at[index, "flank2"] = pdbsequence[int(matched_res2[0])-6:int(matched_res2[0])+5]
            else:
                res_info.at[index, "pdb_adj_res2"] = "Index Error: Matching index not found"
                res_info.at[index, "flank2"] = "Index Error: Matching index not found"
        #store reference sequence as the key of refdict that starts with res_info["ref_org"]
        #store the ref_org associated with each unique_pdb_chain as refdictstart
        #set refdictstart to the value in res_info where res_info["res1_data_ext"] == i
        refdictstart = res_info[res_info["res1_data_ext"] == i]["ref_org"].values[0]
        for key in refdict:
            if key.startswith(refdictstart):
                refsequence = refdict[key]
                #save the key after the "|" character as ref_pdb_chain
                ref_pdb_chain = key.split("|")[1]
        #add ref_pdb_chain to res_info when res1_data_ext == i
        res_info.loc[res_info["res1_data_ext"] == i, "ref_pdb_chain"] = ref_pdb_chain
        #run needle
        run_needle(pdbsequence,refsequence)
        #reformat alignment output
        align_df = alignment_format("align.txt")
        #get reference adjusted index for each res_index_res1 and res_index_res2 where res1_data_ext == i
        for index, row in res_info.iterrows():
            if row["res1_data_ext"] == i:
                #set ref_adj_res1 to the value in b_sequence index where asequence_index == pdb_adj_res1
                res_info.at[index, "adj_res1"] = align_df[align_df["asequence_index"] == row["pdb_adj_res1"]]["bsequence_index"].values[0]
                #set ref_adj_res2 to the value in b_sequence index where asequence_index == pdb_adj_res2
                res_info.at[index, "adj_res2"] = align_df[align_df["asequence_index"] == row["pdb_adj_res2"]]["bsequence_index"].values[0]
    return res_info

#make a new test_col equal to 1 if res_index_res1 is equal to adj_res1 and res_index_res2 is equal to adj_res2 using pandas
#res_info["test_col"] = res_info.apply(lambda row: 1 if row["res_index_res1"] == row["adj_res1"] and row["res_index_res2"] == row["adj_res2"] else 0, axis=1)        

        
final_resinfo = loop_flank2(df, pymol_dict, pdb_dict, reference_dict)


#write final_resinfo to csv


final_resinfo.to_csv(args.output, index=False)
