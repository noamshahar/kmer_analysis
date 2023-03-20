from Bio import SeqIO
import pandas as pd
from tqdm import tqdm
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

    ref_fastq = '/home/noamshahar/kmer_analysis/fastqs/021523_BnS_dPPR2x8swap_200nM_BIOO_trim_Adapt_Size20.fastq'
    ref_fastq_df = parse_fastq(ref_fastq)

