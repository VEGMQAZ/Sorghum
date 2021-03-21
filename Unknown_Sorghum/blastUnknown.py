from Bio import SearchIO
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast import NCBIXML

import os
import ntpath
import argparse
import numpy as np
import pandas as pd
#####################################################################

# Custom database contains 17453 LTR sequences selected using NCBI database query commands.
# Sequences are specific to green plants. Database was created using the following bash commands:
# makeblastdb -in transposon_old_17453.fasta -parse_seqids -dbtype nucl
# (transposon) OR (LTR) AND "green plants"[porgn:__txid33090]

#####################################################################

cwd = os.getcwd()
parser = argparse.ArgumentParser(description='blastUnknown.py -i <file>.fasta')
parser.add_argument('mfilename', metavar='<filename>.fasta', type=str)
args = parser.parse_args()
input = ntpath.basename(args.mfilename)
db_name = "transposon_old_17453.fasta"
sample=open('output_' + ntpath.basename(args.mfilename) + '.txt' , 'a')

LTR_array = ['LTR','Retrotransposon','gag', 'gypsy','Ty3', ' copia-like ', 'copia', 'gypsy-like', 'centromeric retrotranposon']
DNA_array = [' transposon','Transposon', 'Transposable element','transposon-like', 'Rim2', 'CACTA', 'transposase', 'y1-rr', 'hAT', 'Anac', 'TNP2-like','transposable element', 'Mu tranposon']
HEL_array = ['Helitron', 'helitron']

blastx_cline = NcbiblastxCommandline(cmd='blastn', query=input, out=ntpath.basename(input) +'_blastout.txt', outfmt='6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore', db=db_name, evalue=0.01)
result_handle = blastx_cline()

count_Unknown=0
count_LTR=0
count_DNA=0
count_Helitron=0

fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(ntpath.basename(args.mfilename)),'fasta'))
# print(fasta_sequences.keys())

for blast_record in open(ntpath.basename(args.mfilename) + '_blastout.txt'):
    df = (pd.DataFrame(blast_record.split('\t')))
    if any(c in df[0][2] for c in LTR_array):
        if "#Unknown" in df[0][0]:
            count_LTR+=1
    elif any(d in df[0][2] for d in DNA_array):
        if "#Unknown" in df[0][0]:
            count_DNA+=1
    elif any(e in df[0][2] for e in HEL_array):
        if "#Unknown" in df[0][0]:
            count_Helitron+=1
    else:
        count_Unknown+=1
        continue  # Breaking here gives you only the best HSP.

sumTE = count_LTR + count_DNA  + count_Helitron
print("Total reclassified : " + str(sumTE))
print("count_LTR" +"\t"+ "count_DNA" +"\t"+ "count_Helitron" +"\t"+ "count_Unknown")
print(str(count_LTR) +"\t"+ str(count_DNA) +"\t"+ str(count_Helitron) +"\t"+ str(count_Unknown))
