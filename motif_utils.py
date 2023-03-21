# Motif utilites

import fasta_utils
import bash_utils
import os
import pandas as pd
import parse_fastq
from main import is_N_in_seq
from tqdm import tqdm
tqdm.pandas()
script_path = os.getcwd()

def motif_STREME(target_list, ref_list, DNA=True):
    # Recieves a target and ref list of sequences
    # Returns STREME motif results as PSSM DataFrames.
    
    if DNA:
        nc_type = ' --dna '
    else:
        nc_type = ' --rna '
    
    # Write all as temp fasta
    target_list_fasta = os.path.join(script_path, 'STREME_target_temp.fasta')
    ref_list_fasta = os.path.join(script_path, 'STREME_ref_temp.fasta')
    fasta_utils.write_fasta(pd.DataFrame(target_list), target_list_fasta)
    fasta_utils.write_fasta(pd.DataFrame(ref_list), ref_list_fasta)
    
    # Execute command
    result_dir_path = os.path.join(script_path, 'STREME_result_temp')
    bash_utils.write_cmd('streme -p ' + target_list_fasta + ' -n ' + ref_list_fasta + nc_type + '-oc ' + result_dir_path)
    
    # Read results and parse
    pssm_file = os.path.join(result_dir_path, 'streme.txt')
    final_pssm_df = pd.DataFrame(columns=['motif', 'length', 'nsites',
                                          'e-value', 'pssm'])
    with open(pssm_file) as f:
        i = 0
        lines_list = f.readlines()
    
    motif_i = 0
    for i in range(len(lines_list)):
        
        if lines_list[i].startswith('MOTIF'): # New motif PSSM
            motif_line = lines_list[i].split(' ')
            motif_i = int(motif_line[1].split('-')[0])
            motif_seq = motif_line[1].split('-')[1]
            summary_line = lines_list[i+1] # Next line from MOTIF is a summary line of the PSSM
            motif_len = int(summary_line.split(' ')[5]) # Motif length
            motif_nsites = int(summary_line.split(' ')[7]) # Motif nsites
            motif_evalue = float(summary_line.split(' ')[9]) # Motif evalue
            
            
            pssm_lists = lines_list[i+2:i+motif_len+1] # pssm table as a list
            pssm_df = pd.DataFrame(pssm_lists) # Transform pssm list to a DataFrame
            pssm_df = pssm_df[0].str.split(expand=True) # Split into columns
            pssm_df.columns = ['A', 'C', 'G', 'T']
            pssm_df = pssm_df.transpose() # Transpose the pssm table
            final_pssm_df.loc[motif_i] = [motif_seq, motif_len, motif_nsites,
                             motif_evalue, pssm_df] # append current motif information the final pssm df
            i += motif_len + 1 # Skips lines to the n
            
    return final_pssm_df
                

def _get_index_motif_score(pssm, curr_nt, curr_i):
    
    return pssm[curr_i][curr_nt]

def get_motif_score(seq, motif_pssm):
    # Recieves a sequence and a motif pssm
    # Returns the geometric mean of the motif in the sequence
    # Sequence and motif must be at the same length!

    if len(seq) != motif_pssm.shape[1]:
        return 'Sequence and motif are not at the same length!'
    
    seq_split_df = pd.DataFrame(list(seq))
    seq_split_df[1] = range(seq_split_df.shape[0])
    print (seq_split_df)
    a = seq_split_df.apply(lambda row: _get_index_motif_score(motif_pssm, row[0], row[1]))
    return a
    

if __name__ == '__main__':

    fastq_list = ['/home/noamshahar/kmer_analysis/fastqs/021523_BnS_dPPR2x8swap_100nM_BIOO_trim_Adapt_Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/021523_BnS_dPPR2x8swap_200nM_BIOO_trim_Adapt_Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_50nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_100nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_200nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq']
    fastq_ref = '/home/noamshahar/kmer_analysis/fastqs/122021_MMR_Input20mer_BIOO_trim_4at5and3Size20.fastq'
    outpath = '/home/noamshahar/kmer_analysis/streme_files'
    
    fastq_ref_df = parse_fastq.parse_fastq(fastq_ref)
    seq_N_series = fastq_ref_df.seq.progress_apply(is_N_in_seq)
    N_index = seq_N_series[seq_N_series].index.tolist()
    fastq_ref_df = fastq_ref_df[~fastq_ref_df.index.isin(N_index)]

    for curr_fastq in fastq_list:

        curr_csv_name = os.path.basename(curr_fastq).replace('.fastq', '.csv')
        csv_outpath = os.path.join(outpath, curr_csv_name)

        fastq_df = parse_fastq.parse_fastq(curr_fastq)

        # filter sequences with N nt
        print ('filter sequences with N nt...')
        seq_N_series = fastq_df.seq.progress_apply(is_N_in_seq)
        N_index = seq_N_series[seq_N_series].index.tolist()
        fastq_df = fastq_df[~fastq_df.index.isin(N_index)]

        curr_motif_df = motif_STREME(fastq_df.seq.tolist(), fastq_ref_df.seq.tolist(), False)

        curr_motif_df.to_csv(csv_outpath)