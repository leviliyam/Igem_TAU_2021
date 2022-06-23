import typing
from modules.user_IO.input_functions import *
from modules.ORF.TAI import TAI
from modules.ORF.calculating_cai import general_geomean
from modules.logger_factory import LoggerFactory
import os
from statistics import mean, stdev
from pathlib import Path
from BCBio import GFF
# initialize the logger object
logger = LoggerFactory.create_logger("user_input")


class UserInputModule(object):
    @staticmethod
    def get_name() -> str:
        return "User Input"

    @classmethod
    def run_module(cls, user_inp_raw: typing.Dict) -> typing.Dict:
        logger.info('##########################')
        logger.info('# USER INPUT INFORMATION #')
        logger.info('##########################')
        return cls._parse_input(user_inp_raw)

    @classmethod
    def _parse_input(cls, usr_inp):
        '''
        :param usr_inp: in the following format
            {
        'ecoli': {'genome_path': '.gb' - mandatory
                    'genes_HE': '.csv'
                    'optimized': True - mandatory
                    }
         'bacillus': {'genome_path': '.gb'
                    'genes_HE': '.csv'
                    'optimized': False
                    }
            }
        :return: an extended dictionary with the following items:
        @selected_prom : final used list of promoters for MAST
        @sequence : the ORF to optimize
        @organism: contains the org names as keys, and for each organism- the key is the scientific organism name, and the value is _parse_single_input for the organism's
        input
        '''


        full_inp_dict = {}
        full_inp_dict['organisms'] = {}
        for key, val in usr_inp['organisms'].items():
            try:
                org_name, org_dict = cls._parse_single_input(val)
            except:
                raise ValueError(f'Error in input genome: {key}, re-check your input')
            if org_name in full_inp_dict['organisms'].keys():
                raise ValueError(f"Organism: {org_name}'s genome is inserted twice, re-check your input.")
            full_inp_dict['organisms'][org_name] = org_dict  # creating the sub dictionary for each organism- where the key is the scientific name and the value is the following dict:

        # #plots for intergenic sequence analysis: change pron length to 0 and th to 0 as well in the intergenic seuqence code before applying
        #     plt.hist([len(int_seq) for int_seq in list(org_dict['intergenic'].values())],
        #              label=org_name, bins=1000, alpha=0.3, density=True)
        # print(sum([len(int_seq) for int_seq in list(org_dict['intergenic'].values())])/len([len(int_seq) for int_seq in list(org_dict['intergenic'].values())]))
        # plt.title('Intergenic sequence length histogram for promoter analysis')
        # plt.xlabel("Intergenic sequence length")
        # plt.xlim(0,1500)
        # plt.legend()
        # plt.show()


        # add non org specific keys to dict
        orf_fasta_fid = usr_inp['sequence']
        if orf_fasta_fid is not None:
            try:
                orf_seq = str(SeqIO.read(orf_fasta_fid, 'fasta').seq)
            except:
                raise ValueError(
                    f'Error in protein .fasta file: {orf_fasta_fid}, make sure you inserted an undamaged .fasta file containing a single recored')
        prom_fasta_fid = usr_inp['selected_promoters']
        selected_prom = {}
        logger.info(f'\n\nSequence to be optimized given in the following file {orf_fasta_fid}')
        logger.info(f'containing this sequence: {orf_seq}')
        if prom_fasta_fid is not None:
            try:
                selected_prom = fasta_to_dict(prom_fasta_fid)
                logger.info(f'Promoter options ranked are given in the following file {prom_fasta_fid}, '
                            f'which contains {len(selected_prom)} promoters')
            except:
                raise ValueError(
                    f'Error in promoter fasta file: {prom_fasta_fid}, make sure you inserted an undamaged .fasta file')
        else:
            logger.info(
                f'External promoter options were not supplied. endogenous promoters will be used for optimization.'
                f'promoters from the 1/3 most highly expressed genes of all organisms are used- ')
            for org, org_dict in full_inp_dict['organisms'].items():
                if org_dict['optimized']:
                    org_third_he_prom_dict = org_dict['third_most_HE']
                    for prom_name, prom_Seq in org_third_he_prom_dict.items():
                        selected_prom[prom_name + ' from organism: ' + org] = prom_Seq
                    logger.info(f'{len(org_third_he_prom_dict)} promoters are selected from {org}')
            logger.info(f'Resulting in a total of {len(selected_prom)} used for promoter selection and optimization')

        full_inp_dict['selected_prom'] = selected_prom
        full_inp_dict['sequence'] = orf_seq
        full_inp_dict['tuning_param'] = usr_inp['tuning_param']
        return full_inp_dict

    @classmethod
    def run_for_mgnify(cls, user_inp_raw: typing.Dict) -> typing.Dict:
        logger.info('##########################')
        logger.info('# MGNIFY USER INPUT #')
        logger.info('##########################')
        full_inp_dict = {}
        full_inp_dict['organisms'] = {}

        for key, val in user_inp_raw['organisms'].items():
            try:
                org_name, org_dict = cls._parse_single_input_mgnify(val)
            except:
                raise ValueError(f'Error in input genome: {key}, re-check your input')
            full_inp_dict['organisms'][org_name] = org_dict

        # add non org specific keys to dict
        orf_fasta_fid = user_inp_raw['sequence']
        if orf_fasta_fid is not None:
            orf_seq = str(SeqIO.read(orf_fasta_fid, 'fasta').seq)

        selected_prom = {}
        logger.info(
            f'External promoter options were not supplied. endogenous promoters will be used for optimization.'
            f'promoters from the 1/3 most highly expressed genes of all organisms are used- ')

        for org, org_dict in full_inp_dict['organisms'].items():
            if org_dict['optimized']:
                org_third_he_prom_dict = org_dict['third_most_HE']
                for prom_name, prom_Seq in org_third_he_prom_dict.items():
                    selected_prom[prom_name + ' from organism: ' + org] = prom_Seq
                logger.info(f'{len(org_third_he_prom_dict)} promoters are selected from {org}')
        logger.info(f'Resulting in a total of {len(selected_prom)} used for promoter selection and optimization')

        full_inp_dict['selected_prom'] = selected_prom
        full_inp_dict['sequence'] = orf_seq
        full_inp_dict['tuning_param'] = user_inp_raw['tuning_param']
        return full_inp_dict


    @staticmethod
    def _parse_single_input_mgnify(val):
        genome = val['genome_path']
        gff_files = list(Path(genome).glob("*.gff"))
        if len(gff_files) != 1:
            print(F"Found {len(gff_files)} .gff files for {genome}. Abort!")
            exit(1)

        gff_file_path = gff_files[0]
        # for gff_file_path in gff_files:
        fasta_files = list(Path(genome).glob(F"{Path(genome).name}.fna"))
        if len(fasta_files) != 1:
            print(F"Wrong number of .fna files for {genome}..")
            exit(1)

        fasta_file = fasta_files[0]

        with open(fasta_file) as file:
            fasta_dict = {record.description: str(record.seq) for record in SeqIO.parse(file, "fasta")}
        fasta_keys = [x.split("-") for x in fasta_dict.keys()]

        gff_file = GFF.parse(gff_file_path)
        final_prom_dict = {}
        final_intergenic_dict = {}
        final_cds_dict = {}
        for record in gff_file:
            cds_seqs = []
            gene_names = []
            functions = []
            starts = []
            ends = []
            strands = []
            genome = None
            for x in fasta_keys:
                if record.id in x:
                    genome = fasta_dict.get("-".join(x))
                    break
            if genome is None:
                print(
                    F"Failed to find a matching genome sequence for: {record.id} in {file_path}. Skipping all record"
                    F" features.")
                continue
            for feature in record.features:
                if feature.type == "CDS":
                    if "gene" in feature.qualifiers.keys():
                        name = feature.qualifiers['gene'][0]
                    elif "locus_tag" in feature.qualifiers.keys():
                        name = feature.qualifiers["locus_tag"][0]
                    else:
                        continue
                    if feature.location is not None \
                            and name not in gene_names:
                        try:
                            cds = feature.extract(genome)
                        except:
                            print(F"Could not decode cds for feature {name} in {file_path}. Skipping.")
                            continue
                        function = " ".join(feature.qualifiers["product"])

                        if len(cds) % 3 != 0:
                            continue
                        gene_names.append(name)
                        cds_seqs.append(cds)
                        functions.append(function)
                        starts.append(feature.location.start)
                        ends.append(feature.location.end)
                        strands.append(feature.location.strand)
            entry_num = len(gene_names)
            name_and_function = [gene_names[i] + '|' + functions[i] for i in range(entry_num)]
            cds_dict = {name_and_function[i]: cds_seqs[i] for i in range(entry_num)}
            final_cds_dict.update(cds_dict)
            prom200_dict = extract_prom(starts, ends, strands, name_and_function, prom_length=200, genome=genome)
            final_prom_dict.update(prom200_dict)
            intergenic_dict = extract_intergenic(starts, ends, strands, prom_length=200, genome=genome, len_th=30)
            final_intergenic_dict.update(intergenic_dict)

        gene_names = list(final_cds_dict.keys())
        cai_weights = calculate_cai_weights_for_input(final_cds_dict, {}, None)
        cai_scores = general_geomean(sequence_lst=final_cds_dict.values(), weights=cai_weights)
        cai_scores_dict = {gene_names[i]: cai_scores[i] for i in range(len(gene_names))}

        highly_exp_promoters = extract_highly_expressed_promoters(cai_scores_dict, final_prom_dict, percent_used=1/3)

        org_dict = {
            '200bp_promoters': final_prom_dict,
            'third_most_HE': highly_exp_promoters,
            'intergenic': final_intergenic_dict,
            'cai_profile': cai_weights,  # {dna_codon:cai_score}
            'tai_profile': None,  # {dna_codon:tai_score}, if not found in tgcnDB it will be an empty dict.
            'cai_scores': cai_scores_dict,  # {'gene_name': score}
            'tai_scores': {},  # {'gene_name': score}, if not found in tgcnDB it will be an empty dict.
            'optimized': val['optimized'],
            'cai_avg': mean(cai_scores),
            'tai_avg': None,
            'cai_std': stdev(cai_scores),
            'tai_std': None
        }

        return genome, org_dict


    @staticmethod
    def _parse_single_input(val):
        '''
        create the relevant information for each organism
        :param val: the dictionary supplied for every organism:
        {'genome_path': '.gb'
                    'genes_HE': '.csv'
                    'optimized': False }
        :return:
        @org_name = scientific organism name
        @org_dict = {
            'tgcn': {anti codon:number of occurences}, for ORF model
            '200bp_promoters': {gene name and function: prom}, promoter model
            'third_most_HE': same as 200bp_promoters but only third most highly expressed promoters
            'gene_cds':  {gene name and function : cds}, for ORF model
            'intergenic':{position along the genome: intergenic sequence}, promoter model
            'expression_estimation_of_all_genes': {'gene name and function: estimated expression }when the expression csv is not given- the CAI is used as expression levels
            'cai_scores': {'gene_name': cai} ORF and promoter
            'cai_profile': {codon:cai score},  # ORF model
            'tai_scores': {'gene_name': tai} ORF and promoter
            'tai_profile': {codon:cai score},  # ORF model
            'optimized': bool- True if organism is optimized}
        '''
        gb_path = val['genome_path']
        exp_csv_fid = val['expression_csv']
        try:
            gb_file = SeqIO.read(gb_path, format='gb')
        except:
            raise ValueError(
                f'Error in genome GenBank file: {gb_path}, make sure you inserted an undamaged .gb file containing the full genome sequence and annotations')


        org_name = find_org_name(gb_file)
        logger.info(f'\nInformation about {org_name}:')
        if val['optimized']:
            logger.info('Organism is optimized')
        else:
            logger.info('Organism is deoptimized')

        prom200_dict, cds_dict, intergenic_dict, estimated_expression = extract_gene_data(gb_path, exp_csv_fid)
        logger.info(f'Number of genes: {len(cds_dict)}, number of intregenic regions: {len(intergenic_dict)}')
        gene_names = list(cds_dict.keys())
        cai_weights = calculate_cai_weights_for_input (cds_dict, estimated_expression, exp_csv_fid)
        cai_scores = general_geomean(sequence_lst= cds_dict.values(), weights= cai_weights)
        cai_scores_dict = {gene_names[i]:cai_scores[i] for i in range(len(gene_names))}
        #
        # try:
        #     tai_weights = TAI(tai_from_tgcnDB(org_name)).index
        #     tai_scores=general_geomean(sequence_lst= cds_dict.values(), weights= tai_weights)
        #     tai_scores_dict = {gene_names[i]:tai_scores[i] for i in range(len(gene_names))}
        #     tai_mean = mean(tai_scores)
        #     std_tai = stdev(tai_scores)
        # except:
        #     tai_weights = None
        #     tai_mean = None
        #     std_tai = None
        #     tai_scores_dict = {}
        tai_weights = None
        tai_mean = None
        std_tai = None
        tai_scores_dict = {}

        if len(estimated_expression):
            highly_exp_promoters = \
                extract_highly_expressed_promoters(estimated_expression, prom200_dict, percent_used=1/3)
        else:
            highly_exp_promoters = extract_highly_expressed_promoters(cai_scores_dict, prom200_dict, percent_used=1/3)


        org_dict = {
            '200bp_promoters': prom200_dict,
            'third_most_HE': highly_exp_promoters,
            'intergenic': intergenic_dict,
            'cai_profile': cai_weights,  # {dna_codon:cai_score}
            'tai_profile': tai_weights,  # {dna_codon:tai_score}, if not found in tgcnDB it will be an empty dict.
            'cai_scores': cai_scores_dict,  # {'gene_name': score}
            'tai_scores': tai_scores_dict,  # {'gene_name': score}, if not found in tgcnDB it will be an empty dict.
            'optimized': val['optimized'],
            'cai_avg': mean(cai_scores),
            'tai_avg':tai_mean,
            'cai_std':stdev(cai_scores),
            'tai_std': std_tai
            # 'gene_cds': cds_dict,  # cds dict {gene name and function : cds}, for ORF model
            # 'expression_estimation_of_all_genes': estimated_expression, # when the expression csv is not given- the CAI is used as expression levels
        }

        return org_name, org_dict
