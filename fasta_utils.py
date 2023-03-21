# FastaIO
# For writing and parsing Fasta files

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import os
tqdm.pandas()
script_path = os.getcwd()

def write_fasta(writer_df, output_path):
    # Recieves a DataFrame format
    # Each row in the fasta is as follows:
    # '>' index_row
    # row[1]
    # Headlines of each gene will be in the index
    # Sequence will be in in column[0]    

    with open(output_path, 'w') as out:
        for row in writer_df.iterrows():
            out.write('>' +  str(row[0]) + '\n')
            out.write(str(row[1][0]) + '\n')
            
def parse_fasta_file(input_file, sepe=None, take_list=0):
    # Recieves a path to a fasta file
    # Returns as a DataFrame with name and seq columns
    
    fasta_df = pd.DataFrame(columns=['name', 'seq'])
    fasta_parse = SeqIO.parse(open(input_file),'fasta')
    i = 0
    for fasta in fasta_parse:
        
        fasta_id = fasta.id.split(sepe)
        fasta_df.loc[i] = [fasta_id[take_list], str(fasta.seq)]
        i += 1    
    return fasta_df
