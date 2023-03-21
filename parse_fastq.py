from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
from main import is_N_in_seq
import os
import fasta_utils
tqdm.pandas()

def parse_fastq(fastq_path):
    """Parse a fastq file and return a Pandas dataframe with two columns:
    one for read sequences and one for read headers.
    Sequences are filtered to a given length.
    """
    records = SeqIO.parse(fastq_path, "fastq")
    sequences = []
    for record in tqdm(records):
        seq = str(record.seq)
        sequences.append(seq)
        
    df = pd.DataFrame({'seq': sequences})

    return df

if __name__ == '__main__':

    fastq_list = ['/home/noamshahar/kmer_analysis/fastqs/021523_BnS_dPPR2x8swap_100nM_BIOO_trim_Adapt_Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/021523_BnS_dPPR2x8swap_200nM_BIOO_trim_Adapt_Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_50nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_100nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_200nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/122021_MMR_Input20mer_BIOO_trim_4at5and3Size20.fastq']
    outpath = '/home/noamshahar/kmer_analysis/fastq_as_fastas'

    for curr_fastq in fastq_list:

        curr_fasta_name = os.path.basename(curr_fastq).replace('.fastq', '.fasta')
        fasta_outpath = os.path.join(outpath, curr_fasta_name)

        fastq_df = parse_fastq(curr_fastq)

        # filter sequences with N nt
        print ('filter sequences with N nt...')
        seq_N_series = fastq_df.seq.progress_apply(is_N_in_seq)
        N_index = seq_N_series[seq_N_series].index.tolist()
        fastq_df = fastq_df[~fastq_df.index.isin(N_index)]
        
        fastq_df_sample = fastq_df.sample(1000000)

        fasta_utils.write_fasta(fastq_df_sample, fasta_outpath)


