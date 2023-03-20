import pandas as pd
from tqdm import tqdm
import get_ext_files
import os
import numpy as np
tqdm.pandas()


def calculate_pssm(sequences):
    """
    Calculates a PSSM based on a list of sequences with the same length.
    Returns a pandas DataFrame representing the PSSM.
    """
    sequences = [seq for seq in sequences if 'N' not in seq]
    # calculate the frequency matrix
    freq_matrix = pd.DataFrame(index=['A', 'C', 'G', 'T'], columns=range(len(sequences[0])), dtype=float)
    for i in tqdm(range(len(sequences[0]))):
        col = [seq[i] for seq in sequences]
        freq_matrix[i] = pd.Series(col).value_counts(normalize=True)
    
    
    return freq_matrix

def get_prob_per_pos(curr_row, pssm):

    pos = curr_row['pos']
    nt = curr_row['nt']

    return pssm.loc[nt, pos]

def calc_sequence_probability(sequence, pssm):
    """
    Calculate the probability of a given DNA sequence occurring in the PSSM table.

    Args:
        sequence (str): A DNA sequence with the same length as the columns in the PSSM table.
        pssm (pandas.DataFrame): A PSSM table with nucleotides as the rows and positions as the columns.

    Returns:
        float: The probability of the given sequence occurring in the PSSM table.
    """
    # Get the column labels of the PSSM table as integers
    pssm.columns = pssm.columns.astype(int)

    # Convert the DNA sequence to a NumPy array of nucleotides
    seq_df = pd.DataFrame({'pos':range(len(sequence)),
                           'nt':list(sequence)})

    prob_series = seq_df.apply(get_prob_per_pos, args=(pssm,), axis=1)

    # Calculate the total probability by multiplying the probabilities for each position
    total_probability = prob_series.prod()

    return total_probability


if __name__ == '__main__':

    kmer_csv_dir = '/home/noamshahar/kmer_analysis/kmer_csvs/ref_library'
    outpath = '/home/noamshahar/kmer_analysis/pssm_csvs'
    csv_list = get_ext_files.find_ext_files(kmer_csv_dir, 'csv')
    """
    for curr_csv in csv_list:

        print (curr_csv)
        curr_df = pd.read_csv(curr_csv)

        curr_pssm = calculate_pssm(curr_df['seq'].tolist())

        file_name = os.path.basename(curr_csv).replace('.csv', '_pssm.csv')
        curr_outpath = os.path.join(outpath, file_name)

        curr_pssm.to_csv(curr_outpath)
    """
    kmer_csv_cannon = '/home/noamshahar/kmer_analysis/kmer_csvs/200nm/021523_BnS_dPPR2x8swap_200nM_BIOO_trim_Adapt_Size20_kmer_'
    pssm_csv_cannon = '/home/noamshahar/kmer_analysis/pssm_csvs/122021_MMR_Input20mer_BIOO_trim_4at5and3Size20_kmer_'
    outpath = '/home/noamshahar/kmer_analysis/enrichment_csvs'

    for i in range(4,11):

        print (i)
        curr_kmer_path = kmer_csv_cannon + str(i) + '.csv'
        curr_pssm_path = pssm_csv_cannon + str(i) + '_pssm.csv'

        curr_kmer_df = pd.read_csv(curr_kmer_path, index_col=0)
        curr_pssm_df = pd.read_csv(curr_pssm_path, index_col=0)

        sequences = [seq for seq in curr_kmer_df['seq'].tolist() if 'N' not in seq]

        seq_unique = pd.Series(sequences).value_counts(normalize=True)

        prob_in_ref = pd.Series(seq_unique.index.tolist()).progress_apply(calc_sequence_probability, args=(curr_pssm_df,))

        enrichment_df = pd.DataFrame({'seq':seq_unique.index.tolist(),
                                      'prob_in_lib':seq_unique.tolist(),
                                      'prob_in_ref': prob_in_ref.tolist()})

        enrichment_df['enrichment_index'] = enrichment_df['prob_in_lib'] / enrichment_df['prob_in_ref']

        enrichment_df = enrichment_df.sort_values('enrichment_index', ascending=False)

        output_file = os.path.join(outpath, kmer_csv_cannon + 'enrichment_' + str(i) + '.csv')

        enrichment_df.to_csv(output_file)