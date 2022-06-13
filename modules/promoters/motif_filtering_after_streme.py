import csv
import glob
import math
import typing
import random
from collections import defaultdict
from modules.promoters.globals_and_shared_methods import *
from modules.promoters.intersect_motifs_2_org_final import *
from pathlib import Path
import numpy as np

######################################################
###### Motif filtering for multi-organism model ######
######################################################


def is_optimized(fname, input_dict):
    """
    Checks if an intergenic STREME output file is from an optimized organism

    @param fname: name of a STREME output file

    @return: True if the file is from an optimized organism, False otherwise
    """
    fname = str(fname)
    name = fname.split(os.sep)[-2]
    org_name = '_'.join(name.split('_')[:2])
    org_name = org_name.replace("_", " ")
    if org_name in input_dict["organisms"] and input_dict["organisms"][org_name]["optimized"]:
        return True
    return False


def is_org_in_dict(fname, input_dict):
    fname = str(fname)
    name = fname.split(os.sep)[-2]
    org_name = '_'.join(name.split('_')[:2])
    org_name = org_name.replace("_", " ")
    return org_name in input_dict["organisms"]


def find_all_inter_files(base_directory=None):
    """
    Finds all xml STREME output files from intergenic runs

    @return: a list of file names
    """
    if base_directory is None:
        base_directory = Path(os.path.join(os.path.dirname(os.path.abspath(__file__)), start, "streme_outputs", "promoters_motifs"))
    inter_files = base_directory.glob(os.path.join("**", "*.xml"))
    return inter_files


def find_all_anti_motif_files(base_directory=None):
    """
    Finds all xml STREME output files from intergenic runs

    @return: a list of file names
    """
    if base_directory is None:
        base_directory = Path(os.path.join(os.path.dirname(os.path.abspath(__file__)), start, "streme_outputs", "intergenic_anti_motifs"))
    inter_files = base_directory.glob(os.path.join("**", "*.xml"))
    return inter_files


"""
Finds all xml STREME output files from selective runs

@return: a list of file names
"""
def find_all_selective_files():
    all_files = glob.glob(os.path.join(start, 'streme_outputs', "**", "*.xml"), recursive=True)
    selective_files = [f for f in all_files if 'inter' not in f]
    return selective_files


def unionize_motifs(input_dict) -> et.ElementTree:
    """
    Creates a new xml STREME file with all intergenic motifs for all the unwanted hosts

    @return: a pointer to an Element object containing the data in xml format
    """
    inter_files = [str(f) for f in find_all_inter_files() if is_optimized(f, input_dict)]
    base_file = inter_files[0]
    base_tree = et.parse(base_file)
    base_root = base_tree.getroot()
    base_element = base_root.find('.//motifs')

    for xml_file in inter_files[1:]:
        tree = et.parse(xml_file)
        root = tree.getroot()
        element = root.find('.//motifs')
        for motif in element.findall('motif'):
            base_element.append(motif)

    base_tree.write(os.path.join(start, 'unionized_motifs.xml'))
    return base_tree


"""
Finds all intergenic motifs with low correlation to selective ones

@param C_set: dictionary of the unionized intergenic motifs with pssms
@param file_list: a list of files to calculate their motifs' correlation to C_set
@threshold: used to decide which correlation scores are too low

@return: a list of motif indices to delete from the xml file
"""
def get_motifs_to_delete(C_set, file_list, threshold: float):
    to_delete = set()
    for file in file_list:
        if is_optimized(file):
            S_x = extract_pssm_from_xml(file)
            corr_df, pval_df = compare_pssm_sets(C_set, S_x)
            corr_df['max'] = corr_df.max(axis=1)
            low_idx = list(np.where(corr_df['max'] < threshold)[0])
            to_delete.update(low_idx) #mark all motifs with low correlation

    return to_delete


def generate_random_pssms_for_correlation_threshold_calculation():
    min_width = 6
    max_width = 20

    random_pssms_count = 100
    # iterations = math.ceil(random_pssms_count / (max_width - min_width + 1))
    random_pssms = {}
    # for i in range(iterations):
    for i in range(random_pssms_count):
        # for width in range(min_width, max_width + 1):
        width = random.randrange(min_width, max_width + 1)
        df = pd.DataFrame(index=['A', 'C', 'G', 'T'])
        for j in range(width):
            freqs = np.random.dirichlet(np.ones(4), size=1)
            df[j+1] = np.array(freqs[0], dtype=float)
        random_pssms[F"{i}-random"] = df
    return random_pssms


def calculate_experimental_threshold_per_motif_set(motif_files: typing.Sequence[str]) -> typing.Dict[str, float]:
    thresholds_per_motif_file = {}
    for file in motif_files:
        pssms = extract_pssm_from_xml(file)
        random_pssms = generate_random_pssms_for_correlation_threshold_calculation()
        corr_df = pd.DataFrame()
        for random_pssm_ind, random_pssm in random_pssms.items():
            for pssm_ind, pssm in pssms.items():
                corr, pval = compare_pssms(random_pssm, pssm)
                corr_df.loc[random_pssm_ind, pssm_ind] = corr
        thresholds_per_motif_file[file] = np.percentile(corr_df, 95)

    return thresholds_per_motif_file


def compare_with_motifs(candidate_pssms,
                        organisms_motifs_files: typing.Sequence[str],
                        threshold_per_motif_sets: typing.Dict[str, float]) -> typing.Dict[str, int]:
    correlated_organisms = defaultdict(int)

    organisms_count = len(organisms_motifs_files)
    for candidate_pssm_ind, candidate_pssm in candidate_pssms.items():
        for file in organisms_motifs_files:
            single_organism_pssms = extract_pssm_from_xml(file)
            corr_df = pd.DataFrame()
            for pssm_ind, pssm in single_organism_pssms.items():
                corr, pval = compare_pssms(candidate_pssm, pssm)
                corr_df.loc[candidate_pssm_ind, pssm_ind] = corr

            max_correlation = corr_df.max(axis=1)[0]
            if max_correlation >= threshold_per_motif_sets[file]:
                correlated_organisms[candidate_pssm_ind] += 1
        # Normalize count by the number of organisms
        correlated_organisms[candidate_pssm_ind] /= organisms_count

    return correlated_organisms


def multi_filter_motifs(tree, input_dict):
    """
    Creates a final STREME file with only motifs that are both selective and intergenic in all optimized organisms

    @param tree: an ElementTree object containing all intergenic motifs

    @return: the name of the new motif file created
    """
    # TODO modifications:
    #   2. Log experimental thresholds for every run
    #   3. Compare e-values of optimized vs. deoptimized organisms
    #   4. Compare e-values of optimized vs. deoptimized when looking at all promoters (not just HE)
    #   5. Repeat for subset of populations from mgnify DB - calcualte HE genes, calculate promoters, continue with analysis
    unified_motifs = os.path.join(start, 'unionized_motifs.xml')
    inter_files = [str(f) for f in find_all_inter_files() if is_optimized(f, input_dict)]
    anti_motif_files = [str(f) for f in find_all_anti_motif_files() if is_org_in_dict(f, input_dict) and not
                        is_optimized(f, input_dict)]
    candidates_pssms = extract_pssm_from_xml(unified_motifs)

    inter_files_thresholds = calculate_experimental_threshold_per_motif_set(inter_files)
    wanted_organisms_score = compare_with_motifs(candidates_pssms, inter_files, inter_files_thresholds)
    anti_motifs_thresholds = calculate_experimental_threshold_per_motif_set(anti_motif_files)
    unwanted_organisms_score = compare_with_motifs(candidates_pssms, anti_motif_files, anti_motifs_thresholds)

    Path(os.path.join(artifacts_path, "promoters")).mkdir(parents=True, exist_ok=True)
    with open(os.path.join(artifacts_path, "promoters", "experimental_p_value.csv"), "w") as p_value_file:
        csv_writer = csv.writer(p_value_file)
        csv_writer.writerow(["organism_motifs_file", "p_value"])
        for key, value in inter_files_thresholds.items():
            csv_writer.writerow([key, value])
        for key, value in anti_motifs_thresholds.items():
            csv_writer.writerow([key, value])

    final_motifs_score = {}
    tuning_parameter = input_dict["tuning_param"]
    for pssm_ind, pssm in candidates_pssms.items():
        final_motifs_score[pssm_ind] = tuning_parameter * wanted_organisms_score[pssm_ind] + (1 - tuning_parameter) * \
                                       unwanted_organisms_score[pssm_ind]

    final_percentile_cutoff = 75
    percentile_threshold = np.percentile(list(final_motifs_score.values()), final_percentile_cutoff)
    final_motif_set = [pssm_ind for pssm_ind, pssm_score in final_motifs_score.items() if
                       pssm_score >= percentile_threshold]

    for i, motif in enumerate(tree.findall('motif')):
        if i not in final_motif_set:
            tree.remove(motif)
    fname = 'final_motif_set.xml'
    fname = os.path.join(start, fname)
    tree.write(fname)
    
    return fname


def create_final_motif_xml(input_dict):
    """
    Calls the xml pipeline to create final motif set for MAST

    @param D1: threshold for intergenic correlation calculation
    @param D2: threshold for selective correlation calculation

    @return: the name of the new motif file created
    """

    tree = unionize_motifs(input_dict)
    motif_file = multi_filter_motifs(tree, input_dict)
    
    return motif_file

