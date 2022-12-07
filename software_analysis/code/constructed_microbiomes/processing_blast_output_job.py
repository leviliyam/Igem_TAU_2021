# -*- coding: utf-8 -*-
import pandas as pd
import json
from Bio import SeqIO
from os import listdir
from os.path import isfile, join

nt_to_aa = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


def refseq_to_blast_name(refseq_name):
    blast_name = refseq_name.split('.')[0]
    return blast_name

def blast_to_refseq_name(blast_name, full_refseq_list):
    refseq_name = [i for i in full_refseq_list if refseq_to_blast_name(i) == blast_name][0]
    return refseq_name


def entry_to_data_files(entry:str, tls_files:list):
    ###find the data files for each master study
    fasta_files = [tls_dir + f for f in tls_files if 'fsa_nt' in f and entry in f]
    data_files = {fasta:fasta[:-7]+'_5_hits.csv' for fasta in fasta_files}
    return data_files

def tls_metadata(tls_metadata:pd.DataFrame(), tls_files:list):
    blast_results_dict = {}
    for entry_name, row in tls_metadata.iterrows():
        tls_dict = row.to_dict()
        del tls_dict['fasta']
        tls_dict['files'] = entry_to_data_files(entry_name, tls_files)
        blast_results_dict[entry_name] = tls_dict
    return blast_results_dict

def blast_df_to_dict(blast_df, full_refseq_list):
    'dict structure is that the key is the refseq id and the value is the highest percent of identity'
    'between the specific resfeq id (genome) and the tls scores (meaning the highest identity score between all tls '
    '(that correspond to a specific 16s denoted by its refseq id'

    n_seq = len(blast_df)/5
    blast_df.sort_values(by=['pident'], inplace=True)
    blast_df['pident'] = blast_df[blast_df['pident'] > 85]['pident']
    blast_df.dropna(inplace=True)
    blast_df.drop_duplicates(subset=['qseqid'], inplace=True, keep='last')
    blast_df['match_len'] = blast_df['qend']-blast_df['qstart']
    avg_match_len = blast_df['match_len'].mean()

    hit_names = [blast_to_refseq_name(q, full_refseq_list) for q in blast_df['qseqid'].to_list()]
    hit_evalue = blast_df['evalue']
    evalue_scores_dict = dict(zip(hit_names, hit_evalue))
    return evalue_scores_dict, n_seq, avg_match_len


def blastn_run(pident_scores_dict, genomes_df):

    cai_columns = list(nt_to_aa.keys())
    cai_df = genomes_df.drop([i for i in genomes_df.columns if i not in cai_columns], inplace=False, axis = 1)
    non_cai_df = genomes_df.drop(cai_columns, inplace=False, axis=1)
    match_data = {}
    for refseq, score in pident_scores_dict.items():
        cai_values = cai_df.loc[refseq, :].transpose().to_dict()
        other_values_dict = non_cai_df.loc[refseq, :].transpose().to_dict()
        match_data[refseq] = {
            'align_score': score,
            'cai': cai_values
        }
        match_data[refseq].update(other_values_dict)
    return match_data



def check_all_blast_res(genomes_df, tls_metadata:dict):
    blast_results_dict = {}
    for entry, entry_dict in tls_metadata.items():
        blast_results_list = []
        for fasta, blast in entry_dict['files'].items():
            blast_results_list.append(pd.read_csv(blast))
        if len(entry_dict['files'])>1:
            blast_df = pd.concat(blast_results_list)
        elif len(entry_dict['files'])>0:
            blast_df = blast_results_list[0]
        else:
            continue
        blast_df.columns = ['sseqid','qseqid',  'pident','length',
                        'mismatch', 'gapopen', 'qstart', 'qend',
                        'sstart', 'send', 'evalue', 'bitscore']
        full_refseq_list = genomes_df.index.to_list()
        evalue_scores_dict, n_seq, avg_match_len = blast_df_to_dict(blast_df, full_refseq_list)
        match_data = blastn_run(evalue_scores_dict, genomes_df)
        entry_dict['n_seq'] = n_seq
        entry_dict['avg_match_len'] = avg_match_len
        entry_dict['match_data'] = match_data
        print(entry_dict)
    return blast_results_dict


def save_data(blast_results_dict, out_fid):
    with open(out_fid + 'tls_genome_matches_new.json', "w") as outfile:
        json.dump(blast_results_dict, outfile)

if __name__ == "__main__":
    print('Start')

    genomes_df = pd.read_csv('../../data/processed_genomes/filtered/cai_and_16s_for_genomes_filtered.csv', index_col=0)
    genomes_df.drop(columns=['5s', '23s'], inplace=True)
    tls_metadata_df = pd.read_csv('../../data/processed_tls/tls_assembly_metadata.csv', index_col=0)
    tls_dir = '../../data/genbank_tls/'
    tls_files = [f for f in listdir(tls_dir) if isfile(join(tls_dir, f))]
    out_fid = '../../data/tls_genome_match/'
    blast_results_dict= tls_metadata(tls_metadata_df, tls_files)
    blast_results_dict = check_all_blast_res(genomes_df, blast_results_dict)
    save_data(blast_results_dict, out_fid)







