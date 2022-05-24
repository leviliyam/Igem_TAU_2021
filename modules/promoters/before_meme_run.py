import sys
sys.path.append('..')
from modules.promoters.globals_and_shared_methods import *
from modules.shared_functions_and_vars import write_fasta

"""
Creates fasta files from input dictionary

@param data_dict: a dictionary containing all extracted data for each input organism

@return: name of the promoter file created
"""
def create_files_for_meme(data_dict):
    create_folder(opt_path)
    create_folder(deopt_path)
    orgs_dict = data_dict['organisms']
    for org in orgs_dict.keys():
        org_data = orgs_dict[org]
        org = org.replace(' ', '_')
        if org_data['optimized']:
            organism_dict['opt'].append(org)
            org_path = os.path.join(opt_path, org)
        else:
            organism_dict['deopt'].append(org)
            org_path = os.path.join(deopt_path, org)
        create_folder(org_path)
        file_path = os.path.join(org_path, org)
        #inter and 100 not necessary for control organisms???
        write_fasta(file_path + '_100_200', list(org_data['200bp_promoters'].values()), list(org_data['200bp_promoters'].keys()))
        write_fasta(file_path + '_inter',   list(org_data['intergenic'].values()),      list(org_data['intergenic'].keys()))
        write_fasta(file_path + '_33_200',  list(org_data['third_most_HE'].values()),   list(org_data['third_most_HE'].keys()))

    promoter_dict = data_dict['selected_prom']
    new_p_dict = dict()
    for p in promoter_dict.keys():
        new_p_key = p.replace(' ', '_')
        new_p_dict[new_p_key] = promoter_dict[p]

    fname = os.path.join(start, 'promoters')
    write_fasta(fname, list(new_p_dict.values()), list(new_p_dict.keys()))
    return fname + end


def create_unified_promoters_file(data_dict):
    selected_prom = {}
    for org, org_dict in data_dict['organisms'].items():
        if org_dict['optimized']:
            org_third_he_prom_dict = org_dict['third_most_HE']
            for prom_name, prom_Seq in org_third_he_prom_dict.items():
                selected_prom[F"{prom_name}_from_organism_{org}"] = prom_Seq

    fname = os.path.join(start, 'promoters')
    write_fasta(fname, list(selected_prom.values()), list(selected_prom.keys()))
    return fname + end
