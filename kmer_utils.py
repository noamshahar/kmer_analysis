# K-mer utilities

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import os
tqdm.pandas()
script_path = os.getcwd()


def k_mer_split(seq, k):
    # splits a sequence into k-mer running windows
    
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

def fasta_k_mer_split(input_file, k, output_file):
    # Recieves a fasta input_file, k and an output_file path
    # Writes kmer sequences as a fasta file
    
    fasta_parse = SeqIO.parse(open(input_file),'fasta')
    i = 0
    
    with open(output_file, 'w') as outfile:
        for fasta in tqdm(fasta_parse):
            
            curr_name = fasta.id
            curr_seq = str(fasta.seq)
            curr_kmers = k_mer_split(curr_seq, k)
            k_count = 0
            for curr_kmer in curr_kmers:
                outfile.write('>' + curr_name + ':' + str(i) + ':' + str(k_count) + '\n')
                outfile.write(curr_kmer + '\n')
                i += 1
                k_count += 1


if __name__ == '__main__':

    ref_fastq_df = pd.read_csv('/home/noamshahar/kmer_analysis/fastq_as_csvs/021523_BnS_dPPR2x8swap_200nM_BIOO_trim_Adapt_Size20.csv')
    outpath = '/home/noamshahar/kmer_analysis/kmer_csvs/200nm'
    for i in range(4,17):

        print (i)
        curr_kmer_series = ref_fastq_df.seq.apply(k_mer_split, args=(i,))
        curr_kmer_series_explode = curr_kmer_series.explode()
        curr_outpath = os.path.join(outpath, '021523_BnS_dPPR2x8swap_200nM_BIOO_trim_Adapt_Size20_kmer_' + str(i) + '.csv')

        curr_kmer_series_explode.to_csv(curr_outpath)
