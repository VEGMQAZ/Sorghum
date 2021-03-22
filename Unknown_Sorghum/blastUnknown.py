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


LTR_array = ['LTR','Retrotransposon','gag', 'gypsy','Ty3', ' copia-like ', 'copia', 'gypsy-like', 'centromeric retrotranposon']
DNA_array = ['transposon','Transposon', 'Transposable element','transposon-like', 'Rim2', 'CACTA', 'transposase', 'y1-rr', 'hAT', 'Anac', 'TNP2-like','transposable element', 'Mu tranposon']
HEL_array = ['Helitron', 'helitron']

blastx_cline = NcbiblastxCommandline(cmd='blastn', query=input, out=ntpath.basename(input) +'_blastout.txt', outfmt='6 qseqid sseqid salltitles pident length mismatch gapopen qstart qend sstart send evalue bitscore', db=db_name, evalue=0.01)
result_handle = blastx_cline()

count_Unknown=0
count_LTR=0
count_DNA=0
count_Helitron=0

sample=open('reclassified_' + ntpath.basename(args.mfilename) + '.txt' , 'a')

with open('reclassified_' + ntpath.basename(args.mfilename) + '.txt') as fasta_file:  # Will close handle cleanly
    identifiers = []
    lengths = []
    for seq_record in SeqIO.parse(ntpath.basename(args.mfilename), 'fasta'):  # (generator)
        # print(seq_record.id)
        for blast_record in open(ntpath.basename(args.mfilename) + '_blastout.txt'):
            df = (pd.DataFrame(blast_record.split('\t')))
        #     # print(df[0][0])
            if "#Unknown" in df[0][0] and df[0][0] == seq_record.id and any(c in df[0][2] for c in LTR_array):
                count_LTR+=1
                print(">" + seq_record.id.split("#")[0] + "#LTR", file=sample)
                print(seq_record.seq, file=sample)
            elif "#Unknown" in df[0][0] and df[0][0] == seq_record.id and any(d in df[0][2] for d in DNA_array):
                count_DNA+=1
                print(">" + seq_record.id.split("#")[0] + "#DNA", file=sample)
                print(seq_record.seq, file=sample)
            elif "#Unknown" in df[0][0] and df[0][0] == seq_record.id and any(e in df[0][2] for e in HEL_array):
                count_Helitron+=1
                print(">" + seq_record.id.split("#")[0]  + "#RC/Helitron", file=sample)
                print(seq_record.seq, file=sample)
            else:
                # count_Unknown+=1
                continue  # Breaking here gives you only the best HSP.

sample.close()

sumTE = count_LTR + count_DNA  + count_Helitron
print("Total reclassified : " + str(sumTE))
print("count_LTR" +"\t"+ "count_DNA" +"\t"+ "count_Helitron")
print(str(count_LTR) +"\t"+ str(count_DNA) +"\t"+ str(count_Helitron))
