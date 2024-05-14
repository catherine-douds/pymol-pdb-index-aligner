#this script reads in a csv and checks the formatting of the molecule and organism names

import argparse
import pandas as pd
import csv

#set arguments
parser = argparse.ArgumentParser(description='This script will check the formatting of the molecule and organism names in the GU csv from Sharear')
parser.add_argument('i', help='input csv file')
args= parser.parse_args()

#read in csv file
df = pd.read_csv(args.i)

#print column names
print(df.columns)

#check formatting of molecule and organism names
def get_molecules(df):
    molecules = df["Molecule"].unique()
    #sort molecules alphabetically
    molecules = sorted(molecules)
    return molecules



#replace "23 S rRNA with 23S rRNA" in df
#df["Molecule_res1"] = df["Molecule_res1"].str.replace("23 S rRNA", "23S rRNA")



#within the same PDB_ID, sort by chain_ID_res1
#df = df.sort_values(by = ["PDB_ID", "chain_ID_res1"])
#sort dataframe so that anything that contains rRNA in "molecule_res1" is at the end
#df_else = df[~df["Molecule_res1"].str.contains("rRNA")]
#df_rRNA = df[df["Molecule_res1"].str.contains("rRNA")]
#df_sorted = pd.concat([df_else, df_rRNA])

#save df_sorted as csv
#df_sorted.to_csv("r1_GU_sorted.csv", index = False)

#replace "18GAAA (52-MER)"" with "18GAAA(52-MER)"
df["Molecule"] = df["Molecule"].str.replace("23 S rRNA", "23S rRNA")
df["Molecule"] = df["Molecule"].str.replace("5.8 S rRNA", "5.8S rRNA")
df["Molecule"] = df["Molecule"].str.replace("C-di-GMP riboswitch", "c-di-GMP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("c-di-GMP Riboswitch", "c-di-GMP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("E. coli yybP-ykoY riboswitch", "yybP-ykoY riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Core double ENE RNA (Xtal construct) from Oryza sativa transposon,Core double ENE RNA (Xtal construct) from Oryza sativa transposon", "Core double ENE RNA (Xtal construct)")
df["Molecule"] = df["Molecule"].str.replace("Bacillus subtilis small cytoplasmic RNA (scRNA),RNA", "small cytoplasmic RNA (scRNA)")
df["Molecule"] = df["Molecule"].str.replace("Bacillus subtilis xpt", "xpt")
df["Molecule"] = df["Molecule"].str.replace("Escherichia coli strain K-12 substr. MG1655_TMP32XR1 chromosome, complete genome", "EMG1655_TMP32XR1 chromosome")
df["Molecule"] = df["Molecule"].str.replace("C92U mutant c-di-GMP riboswitch", "c-di-GMP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("E.coli RNase P RNA", "RNase P RNA")
df["Molecule"] = df["Molecule"].str.replace("G20A mutant c-di-GMP Riboswitch", "c-di-GMP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("G22A mutant c-di-GMP Riboswitch", "c-di-GMP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("GLMS RIBOZYME", "GlmS ribozyme")
df["Molecule"] = df["Molecule"].str.replace("GLMS ribozyme", "GlmS ribozyme")
df["Molecule"] = df["Molecule"].str.replace("Group I intron, 5 prime fragment", "Group I intron 5 prime fragment")
df["Molecule"] = df["Molecule"].str.replace("Group I intron, 3 prime fragment", "Group I intron 3 prime fragment")
df["Molecule"] = df["Molecule"].str.replace("Guanine Riboswitch", "Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Guanine riboswitch", "Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Guanine riboswitch C74U mutant", "Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("GGuanine riboswitch aptamer domain", "Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Guanine riboswitch aptamer domain", "Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("L-glutamine riboswitch (58-MER)", "L-glutamine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("L-glutamine riboswitch RNA (61-MER)", "L-glutamine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Lactococcus lactis yybP-ykoY riboswitch", "yybP-ykoY riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Lactococcus lactis yybP-ykoY riboswitch", "yybP-ykoY riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Lysine riboswitch  RNA", "Lysine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Lysine Riboswitch RNA", "Lysine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Lysine riboswitch RNA", "Lysine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("M-Box RNA, ykoK riboswitch aptamer", "M-Box Riboswitch")
df["Molecule"] = df["Molecule"].str.replace("PRPP Riboswitch", "PRPP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("PRPP riboswitch", "PRPP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("RIBONUCLEASE P", "RNase P RNA")
df["Molecule"] = df["Molecule"].str.replace("RIBONUCLEASE P RNA", "RNase P RNA")
df["Molecule"] = df["Molecule"].str.replace("RNA (83-MER),RNA (83-MER)", "RNA (83-MER)")
df["Molecule"] = df["Molecule"].str.replace("SAM-I RIBOSWITCH", "SAM-I riboswitch")
df["Molecule"] = df["Molecule"].str.replace("TPP-specific riboswitch", "TPP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("TPP riboswitch", "TPP riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Vibrio Vulnificus Adenine Riboswitch", "Adenine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Vibrio vulnificus A-riboswitch", "Adenine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Vibrio vulnificus Adenine Riboswitch", "Adenine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Vibrio vulnificus strain 93U204 chromosome II, adenine riboswitch aptamer domain", "Adenine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("adenine riboswitch aptamer variant", "Adenine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("X. oryzae Mn riboswitch optimized construct", "Mn riboswitch optimized construct")
df["Molecule"] = df["Molecule"].str.replace("X. oryzae Mn riboswitch optimized construct", "Mn riboswitch optimized construct")
df["Molecule"] = df["Molecule"].str.replace("glmS glucosamine-6-phosphate activated ribozyme", "GlmS ribozyme")
df["Molecule"] = df["Molecule"].str.replace("glmS ribozyme RNA", "GlmS ribozyme")
df["Molecule"] = df["Molecule"].str.replace("group II intron, domain 1", "group II intron domain 1")
df["Molecule"] = df["Molecule"].str.replace("guanine riboswitch", "Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("guanine riboswitch aptamer", "Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("tetrahydrofolate riboswitch aptamer domain", "Tetrahydrofolate riboswitch")
df["Molecule"] = df["Molecule"].str.replace("tetrahydrofolate riboswitch aptamer", "Tetrahydrofolate riboswitch")
df["Molecule"] = df["Molecule"].str.replace("The Tet-S1 state molecule of co-transcriptional folded G264A mutated Tetrahymena group I intron with 6nt 3'/5'-exon and 2-aminopurine nucleoside", "Group I self-splicing intron")
df["Molecule"] = df["Molecule"].str.replace("The Tet-S2 state molecule of co-transcriptional folded wild-type Tetrahymena group I intron with 30nt 3'/5'-exon (intron and 3'-exon)", "Group I self-splicing intron")
df["Molecule"] = df["Molecule"].str.replace("The Tet-S2 state with a pseudoknotted 4-way junction molecule of co-transcriptional folded wild-type Tetrahymena group I intron with 30nt 3'/5'-exon (intron and 3'-exon)", "Group I self-splicing intron")
df["Molecule"] = df["Molecule"].str.replace("The pre-Tet-C state molecule of co-transcriptional folded wild-type Tetrahymena group I intron with 30nt 3'/5'-exon", "Group I self-splicing intron")
df["Molecule"] = df["Molecule"].str.replace("RNase P RNA RNA","RNase P RNA")
df["Molecule"] = df["Molecule"].str.replace("TRANSFER-RNA, PHE","tRNA phe")
df["Molecule"] = df["Molecule"].str.replace("SRS2 fragment of Rgs4 3' UTR, RNA (31-MER)","SRS2 fragment of Rgs4 3' UTR")
df["Molecule"] = df["Molecule"].str.replace("A24U mutant of the B. subtilis xpt-pbuX Guanine riboswitch","Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("A24U/U25A/A46G mutant of the B. subtilis xpt-pbuX Guanine riboswitch","Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("A24U/U25A/A46G/C74U mutant of the B. subtilis xpt-pbuX Guanine riboswitch","Guanine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Adenosine RIboswitch","Adenosine riboswitch")
df["Molecule"] = df["Molecule"].str.replace("Lysine riboswitch RNA","Lysine riboswitch")


molecules = get_molecules(df)
#save molecules as csv where each molecule is new line
with open("molecules.csv", "w") as f:
    for molecule in molecules:
        f.write(molecule + "\n")

#sort df to have all rRNA at the front
df_rRNA = df[df["Molecule"].str.contains("rRNA")]
df_else = df[~df["Molecule"].str.contains("rRNA")]
df_sorted = pd.concat([df_rRNA, df_else])
df_sorted.to_csv("r1_red_sorted.csv", index = False)
df_rRNA.to_csv("r1_red_rRNA.csv", index = False)
df_else.to_csv("r1_red_else.csv", index = False)