#!/powerapps/share/centos7/miniconda/miniconda3-4.7.12-environmentally/envs/Python3.9Plus/bin/python
import sys
sys.path.insert(1, '../')
from pairwise_align_function import align_seq_to_fasta
METAGENOME_SEQ = '../../../../data/tested_results/KDVY_example_metagenome/KDVY_first_entry_KDVY01000001.1.fasta'
GENOME_SLICE = '../../../../data/tested_results/genomes_16s_sliced/18_start_4500_end_4750.fasta'
OUT_DIR = '../../../../../data/tested_results/KDVY_example_metagenome/sliced_alignment_first_entry/'
align_seq_to_fasta(METAGENOME_SEQ, GENOME_SLICE, OUT_DIR)
