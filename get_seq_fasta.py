import pandas as pd
import fasta_utils
import get_ext_files
import os
import subprocess

sd_thres = 2
muscle_path = '/home/noamshahar/kmer_analysis/scripts/muscle3.8.31_i86linux64'
outpath = '/home/noamshahar/kmer_analysis/enrichment_fastas'
enrichment_path = '/home/noamshahar/kmer_analysis/enrichment_csvs/200nm'
enrichment_files = get_ext_files.find_ext_files(enrichment_path, 'csv')
outpath_muscle = '/home/noamshahar/kmer_analysis/muscle_fastas'
for curr_csv in enrichment_files:

    enrichment_df = pd.read_csv(curr_csv, index_col=0)

    enrichment_mean = enrichment_df.enrichment_index.mean()
    enrichment_std = enrichment_df.enrichment_index.std()

    enrichment_std_thres = enrichment_mean + (enrichment_std * sd_thres)

    enrichment_df_filter = enrichment_df[enrichment_df.enrichment_index >= enrichment_std_thres]

    fasta_df = enrichment_df_filter[['seq']]

    curr_outpath = os.path.join(outpath, os.path.basename(curr_csv).replace('.csv', '_enriched_sd_' + str(sd_thres)) + '.fasta')

    fasta_utils.write_fasta(fasta_df, curr_outpath)

    muscle_fasta_file = os.path.basename(curr_outpath).replace('.fasta', '_muscle.fasta')
    muscle_fasta_outpath = os.path.join(outpath_muscle, muscle_fasta_file)

    muscle_cmd = [muscle_path, '-in', curr_outpath, '-out', muscle_fasta_outpath]

    subprocess.run(muscle_cmd)

    logo_cmd = ['weblogo', '-f', muscle_fasta_outpath, '-o',
                 muscle_fasta_outpath.replace('.fasta', '.png'),
                 '--fineprint', '""', '--size', 'large', '--errorbars', 'False',
                 '--scale-width', 'False', '-F', 'png', '--resolution', '300']
    
    subprocess.run(logo_cmd)