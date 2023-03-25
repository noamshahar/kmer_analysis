import pandas as pd
import get_ext_files
import fasta_utils
import bash_utils
import calculate_pssm
import kmer_utils
import parse_fastq
import os
import glob_param
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
tqdm.pandas()

def is_N_in_seq(seq:str) -> bool:
    """returns True if N is found in seq, else False

    Args:
        seq (str): provide DNA string

    Returns:
        bool: True if N found, otherwise False
    """
    # if N in sequence (address as lowercase)
    if 'n' in seq.lower():
        bool_n = True
    else:
        bool_n = False

    return bool_n


def get_ref_pssm(ref_fastq_path:str, k:int) -> pd.DataFrame:
    """main function for obtaining reference library PSSM table

    Args:
        ref_fastq_path (str): path to reference (unbound) library fastq
        k (int): k used for analysis

    Returns:
        pd.DataFrame: PSSM table ref kmer
    """    
    # parse reference fastq to a pandas dataframe
    print ('parse reference fastq to a pandas dataframe...')
    ref_fastq_df = parse_fastq.parse_fastq(ref_fastq_path)

    # filter sequences with N nt
    print ('filter sequences with N nt...')
    seq_N_series = ref_fastq_df.seq.progress_apply(is_N_in_seq)
    N_index = seq_N_series[seq_N_series].index.tolist()
    ref_fastq_df = ref_fastq_df[~ref_fastq_df.index.isin(N_index)]

    print ('create kmer for reference fastq...')
    # create kmer series for ref_fastq_df
    ref_kmer_series = ref_fastq_df.seq.progress_apply(kmer_utils.k_mer_split, args=(k,))
    # explode all lists in all rows to a single pandas series
    ref_kmer_series_explode = ref_kmer_series.explode()

    print ('calculate pssm table for ref_kmer_series_explode')
    # calculate pssm table for ref_kmer_series_explode
    pssm_table = calculate_pssm.calculate_pssm(ref_kmer_series_explode.tolist())

    return pssm_table

def get_enrichment_index(fastq_path:str, k:int, pssm_table:pd.DataFrame) -> pd.DataFrame:
    """calculates enrichmet indices on a single positive fastq file
    providing k and pssm table of reference library

    Args:
        fastq_path (str): path to a single fastq file
        k (int): kmer used in analysis
        pssm_table (pd.DataFrame): pssm table of reference library

    Returns:
        pd.DataFrame: enrichment_df
    """
    # parse fastq to a pandas dataframe
    print ('parse fastq to a pandas dataframe...')
    fastq_df = parse_fastq.parse_fastq(fastq_path)

    # filter sequences with N nt
    print ('filter sequences with N nt...')
    seq_N_series = fastq_df.seq.progress_apply(is_N_in_seq)
    N_index = seq_N_series[seq_N_series].index.tolist()
    fastq_df = fastq_df[~fastq_df.index.isin(N_index)]

    # slice seqs
    fastq_df['seq'] = fastq_df['seq'].apply(lambda x:x[5:16])

    print ('create kmer for fastq...')
    # create kmer series for ref_fastq_df
    kmer_series = fastq_df.seq.progress_apply(kmer_utils.k_mer_split, args=(k,))
    # explode all lists in all rows to a single pandas series
    kmer_series_explode = kmer_series.explode()
    
    print ('calculate probabilities based on pssm table...')
    # get unique kmers and their ratio in the data
    kmer_unique = kmer_series_explode.value_counts(normalize=True)

    # check the probability to find each unique kmer in the ref library, based on the pssm table
    prob_in_ref = pd.Series(kmer_unique.index.tolist()).progress_apply(calculate_pssm.calc_sequence_probability,
                                                                       args=(pssm_table,))

    print ('calculate enrichment index...')
    # add all data to enrichment_df
    enrichment_df = pd.DataFrame({'seq':kmer_unique.index.tolist(),
                                    'prob_in_lib':kmer_unique.tolist(),
                                    'prob_in_ref': prob_in_ref.tolist()})

    # calculate enrichment index per each seq
    enrichment_df['enrichment_index'] = enrichment_df['prob_in_lib'] / enrichment_df['prob_in_ref']

    # sort enrichment_df according to the enrichment_index, top to bottom
    enrichment_df = enrichment_df.sort_values('enrichment_index', ascending=False)

    return enrichment_df

def multiply_seq_by_enrichment(curr_row):

    curr_seq = [curr_row['seq']]
    curr_enrichment = int(curr_row['enrichment_index'])

    return curr_seq * curr_enrichment


def save_enrichment_fasta(enrichment_df:pd.DataFrame, sd_above_mean:float,
                          k:int, outpath:str) -> str:
    """saves enrichment fasta sequences based on sd above threshold param
    also saved histogram of enrichment

    Args:
        enrichment_df (pd.DataFrame): dataframe of kmers with enrichment index
        sd_above_mean (float): threshold of STD above the mean to filter kmers
        outpath (str): output directory of temp fasta

        Returns:
        str: path to enrichment fasta file
    """
    # define histplot figure outpath
    fig_outpath = os.path.join(outpath, 'k_' + str(k) + '_sd_' + str(sd_above_mean) + '_enrichment_hist.png')
   
    # define temp fasta outpath
    fasta_outpath = os.path.join(outpath, 'k_' + str(k) + '_sd_' + str(sd_above_mean) + '_enrichment_fasta.fasta')
   
    # calculate the enrichment_index threshold value
    std_thres = enrichment_df.enrichment_index.mean() + (sd_above_mean * enrichment_df.enrichment_index.std())

    # filter kmers with higher enrichment index than std_thres
    enrichment_filter_df = enrichment_df[enrichment_df.enrichment_index >= std_thres]

    # multiply each by its cognate enrichment index
    multiple_seq = enrichment_filter_df.apply(multiply_seq_by_enrichment, axis=1)
    if multiple_seq.empty:
        enrichment_filter_explode = pd.Series()
    else:
        # explode
        enrichment_filter_explode = multiple_seq.explode()
    # save as a global variable current number of kmers
    global N_KMERS
    N_KMERS = enrichment_filter_explode.shape[0]

    # use round enrichment_index values in the histplot
    ###??###
    bin_edges = range(int(min(enrichment_df['enrichment_index'])), int(max(enrichment_df['enrichment_index'])) + 2)
    # plot enrichmet hist plot and save it
    sns.histplot(data=enrichment_df, x='enrichment_index', bins=bin_edges)
    # add threshold line for std thresold
    plt.axvline(std_thres, linestyle='--', color='red')  
    # add a legend with meta data and number of sequences
    legend_str = f'k={k}\nSD={sd_above_mean}\nvalue={std_thres:.2f}\nn={N_KMERS}'
    plt.legend([legend_str])
    # save figure as png
    plt.savefig(fig_outpath, format='png', dpi=300)
    # close figure
    plt.close()

    # replace 'T' with 'U' in each sequence
    enrichment_filter_explode = enrichment_filter_explode.apply(lambda sequence: sequence.replace('T', 'U'))

    # export fasta to fasta_outpath
    fasta_utils.write_fasta(pd.DataFrame(enrichment_filter_explode), fasta_outpath)

    return fasta_outpath

def get_muscle_fasta(enrichment_fasta_path:str, outpath:str) -> str:
    """performs multiple sequence alignment on enrichment fasta kmers
    uses MUSCLE

    Args:
        enrichment_fasta_path (str): path to enriched kmers fasta
        outpath (str): output directory for muscle fasta

    Returns:
        str: the outpath of the muscle fasta
    """
    # define the muscle outpath fasta file
    muscle_fasta_path = os.path.basename(enrichment_fasta_path).replace('.fasta', '_muscle.fasta')

    # define the full path given outpath
    muscle_fasta_outpath = os.path.join(outpath, muscle_fasta_path)

    # define the command as a list ready for subprocess
    muscle_cmd = [glob_param.MUSCLE_PATH, '-in', enrichment_fasta_path,
                  '-out', muscle_fasta_outpath]

    # run command using subprocess
    subprocess.run(muscle_cmd)

    return muscle_fasta_outpath

def create_logo(muscle_fasta_path:str, outpath:str, k:int,
                sd_above_mean:float):
    """creates a weblogo on muscle fasta

    Args:
        muscle_fasta_path (str): path to muscle fasta of enriched kmers
        outpath (str): output directory
        k (int): kmer used in analysis
        sd_above_mean (float): threshold of STD above the mean to filter kmers

    """
    global N_KMERS

    # define meta data of analysis in the figure
    fineprint = f'k={k}; SD={sd_above_mean}; n={N_KMERS}'

    # define the weblogo outpath file
    logo_path = os.path.basename(muscle_fasta_path).replace('_muscle.fasta', '_logo.eps')

    # define the full path given outpath
    logo_outpath = os.path.join(outpath, logo_path)

    # define the command as a list ready for subprocess
    logo_cmd = [glob_param.WEBLOGO_PATH, '-f', muscle_fasta_path,
                '-o', logo_outpath, '--fineprint', fineprint,
                '--size', 'large', '--errorbars', 'NO',
                 '--scale-width', 'NO',
                 '--resolution', '600',
                 '-c', 'classic',
                 '--sequence-type', 'rna',
                 '--number-interval', '1']

    # run command using subprocess
    subprocess.run(logo_cmd)
   

if __name__ == '__main__':

    positive_fastq_list = ['/home/noamshahar/kmer_analysis/fastqs/021523_BnS_dPPR2x8swap_100nM_BIOO_trim_Adapt_Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/021523_BnS_dPPR2x8swap_200nM_BIOO_trim_Adapt_Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_50nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_100nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq',
                           '/home/noamshahar/kmer_analysis/fastqs/082621_MMR_200nM_dPPR8x2_BIOO_trim_4at5and3Size20.fastq']
    ref_fastq_path = '/home/noamshahar/kmer_analysis/fastqs/122021_MMR_Input20mer_BIOO_trim_4at5and3Size20.fastq'
    outpath = '/home/noamshahar/kmer_analysis/k_results_middle_part'
    k_list = [7]
    sd_above_mean_list = [5,7,10,15,20,30]


    for k in k_list:

        print (k)
        # define current dir as current iterated k
        curr_dir_name = 'run_k_' + str(k)
        # define full dir path and create dir
        curr_outpath_dir = os.path.join(outpath, curr_dir_name)
        print (curr_outpath_dir)
        os.makedirs(curr_outpath_dir, exist_ok=True)

        # calculate pssm table for current iterated k
        pssm_table = get_ref_pssm(ref_fastq_path, k)

        # iterate through all positive fastqs in current k
        for positive_fastq_path in positive_fastq_list:

            # create sub-folder for current fastq
            curr_fastq_dir = os.path.join(curr_outpath_dir, os.path.basename(positive_fastq_path))
            os.makedirs(curr_fastq_dir)
            # calculate enrichment index
            enrichment_df = get_enrichment_index(positive_fastq_path, k, pssm_table)

            # iterate through the list of sd_above_mean_list
            for curr_sd_above_mean in sd_above_mean_list:
                
                # define and create subdir for each sd_above_mean
                curr_sd_k_dir = os.path.join(curr_fastq_dir, 'sd_' + str(curr_sd_above_mean))
                os.makedirs(curr_sd_k_dir, exist_ok=True)

                # save current enrichment fasta
                enrichment_fasta_path = save_enrichment_fasta(enrichment_df, curr_sd_above_mean, k, curr_sd_k_dir)
                # run muscle and save its fasta
                muscle_fasta_path = get_muscle_fasta(enrichment_fasta_path, curr_sd_k_dir)
                # create logo for current muscle fasta
                create_logo(muscle_fasta_path, curr_sd_k_dir, k, curr_sd_above_mean)