import os
import shutil
import time
import traceback
from pathlib import Path
import typing

from modules.logger_factory import LoggerFactory

# Create clean artifacts directory
artifacts_directory = Path(os.path.join(str(Path(__file__).parent.resolve()), "artifacts"))
# if artifacts_directory.exists() and artifacts_directory.is_dir():
#     shutil.rmtree(artifacts_directory)
artifacts_directory.mkdir(parents=True, exist_ok=True)

from modules import Zscore_calculation, user_IO, RE, ORF, promoters
from modules import models

logger = LoggerFactory.create_logger("main")

current_directory = Path(__file__).parent.resolve()
# base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
# base_path = "/arabidopsis_microbiome"
base_path = "C:\\Users\\Kama\\Documents\\Moran\\biomedical-engineering\\microbiome-optimization\\arabidopsis_microbiome"

default_user_inp_raw = {
    'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
    'selected_promoters': None,
    'tuning_param': 0.75,
    'organisms': {
        'org1': {'genome_path': os.path.join(base_path, 'Agromyces allii.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org2': {'genome_path': os.path.join(base_path, 'Arthrobacter crystallopoietes.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org3': {'genome_path': os.path.join(base_path, 'Arthrobacter luteolus.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org4': {'genome_path': os.path.join(base_path, 'Arthrobacter pascens.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org5': {'genome_path': os.path.join(base_path, 'Arthrobacter subterraneus.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org6': {'genome_path': os.path.join(base_path, 'Arthrobacter tumbae.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org7': {'genome_path': os.path.join(base_path, 'Brevibacterium frigoritolerans.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org8': {'genome_path': os.path.join(base_path, 'Janibacter limosus.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org9': {'genome_path': os.path.join(base_path, 'Knoellia subterranea.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org10': {'genome_path': os.path.join(base_path, 'Mycolicibacterium smegmatis.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org11': {'genome_path': os.path.join(base_path, 'Nocardioides daejeonensis.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org12': {'genome_path': os.path.join(base_path, 'Nocardioides jensenii.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org13': {'genome_path': os.path.join(base_path, 'Nocardioides oleivorans.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org14': {'genome_path': os.path.join(base_path, 'Nocardioides sediminis.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org15': {'genome_path': os.path.join(base_path, 'Nocardioides terrigena.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org16': {'genome_path': os.path.join(base_path, 'Paenarthrobacter nitroguajacolicus.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org17': {'genome_path': os.path.join(base_path, 'Paenibacillus aceris.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org18': {'genome_path': os.path.join(base_path, 'Paenibacillus alginolyticus.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org19': {'genome_path': os.path.join(base_path, 'Paenibacillus oryzisoli.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org20': {'genome_path': os.path.join(base_path, 'Paenibacillus prosopidis.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org21': {'genome_path': os.path.join(base_path, 'Paenibacillus qinlingensis.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org22': {'genome_path': os.path.join(base_path, 'Pedococcus badiiscoriae.gb'),
                 'optimized': False,
                 'expression_csv': None,
                 },
        'org23': {'genome_path': os.path.join(base_path, 'Pedococcus bigeumensis.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org24': {'genome_path': os.path.join(base_path, 'Pedococcus dokdonensis.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org25': {'genome_path': os.path.join(base_path, 'Peribacillus muralis.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org26': {'genome_path': os.path.join(base_path, 'Peribacillus simplex.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org27': {'genome_path': os.path.join(base_path, 'Phycicoccus duodecadis.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org28': {'genome_path': os.path.join(base_path, 'Priestia flexa.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org29': {'genome_path': os.path.join(base_path, 'Pseudarthrobacter phenanthrenivorans.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org30': {'genome_path': os.path.join(base_path, 'Rhodanobacter denitrificans.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org31': {'genome_path': os.path.join(base_path, 'Rhodanobacter fulvus.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org32': {'genome_path': os.path.join(base_path, 'Terrabacter aerolatus.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org33': {'genome_path': os.path.join(base_path, 'Terrabacter tumescens.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
        'org34': {'genome_path': os.path.join(base_path, 'Yonghaparkia alkaliphila.gb'),
                  'optimized': False,
                  'expression_csv': None,
                  },
    }
}


def create_csv_for_organism(csv_directory_name: str,
                            csv_file_name: str,
                            organism_name: str,
                            mast_file_name: str,
                            user_input):
    import csv
    import xml.etree.ElementTree as et

    base_directory = os.path.join(artifacts_directory, "promoters_not_for_user")
    csv_dir_path = os.path.join(base_directory, csv_directory_name)
    Path(csv_dir_path).mkdir(parents=True, exist_ok=True)
    csv_file_path = os.path.join(csv_dir_path, csv_file_name)
    organism_info = user_input["organisms"][organism_name]
    with open(csv_file_path, "w") as csv_file:
        csv_writer = csv.writer(csv_file)
        header = ["gene_name", "gene_cai_score", "gene_promoter_evalue"]
        csv_writer.writerow(header)

        base_file = os.path.join(base_directory, mast_file_name)
        tree = et.parse(base_file)
        root = tree.getroot()

        cai_scores = []
        evalues = []
        for m in root.findall(".//sequence"):
            seq_name = m.get("name")
            if "hypothetical" in seq_name:
                logger.info(F"skipping hypothetical match {seq_name}")
                continue
            promoter_score = m.find("score")
            promoter_evalue = promoter_score.get("evalue")
            matching_genes = [k for k in organism_info['cai_scores'].keys() if k.startswith(seq_name)]
            if len(matching_genes) > 1:
                logger.error(F"Found more than one matching gene for promoter! {matching_genes}")
                exit(1)
            matching_gene = matching_genes[0]
            cai_score = organism_info["cai_scores"][matching_gene]

            # Filter only results of highly expressed genes
            percent_used = 1/3
            expression_list = list(organism_info["cai_scores"].values())
            expression_list.sort(reverse=True)
            expression_threshold = expression_list[round(len(expression_list) * percent_used)]

            # if cai_score < expression_threshold:
            #     logger.info(F"Skipping gene {matching_gene} because its cai score {cai_score} is lower than threshold "
            #                 F"{expression_threshold}")
            #     continue

            data = [matching_gene, cai_score, promoter_evalue]
            csv_writer.writerow(data)

            cai_scores.append(cai_score)
            evalues.append(promoter_evalue)

    return cai_scores, evalues


def create_csv_for_organism_selective(first_organism, second_organism, user_input, is_first: bool):
    suffix = "first" if is_first else "second"
    organism_name = first_organism if is_first else second_organism
    csv_directory_name = "selective_csv_files"
    csv_file_name = F"opt_{first_organism}_deopt_{second_organism}_{suffix}.csv"
    mast_file_name = F"mast {first_organism} {second_organism} {suffix}\\mast.xml".replace(" ", "_")
    return create_csv_for_organism(csv_directory_name=csv_directory_name,
                                   csv_file_name=csv_file_name,
                                   organism_name=organism_name,
                                   mast_file_name=mast_file_name,
                                   user_input=user_input)


def create_csv_for_organism_intergenic(organism_name: str, user_input):
    csv_file_name = F"{organism_name}.csv"
    csv_directory_name = "intergenic_csv_files"
    mast_file_name = F"intergenic_mast\\motif streme seq {organism_name} 100 200\\mast.xml".replace(" ", "_")
    return create_csv_for_organism(csv_directory_name=csv_directory_name,
                                   csv_file_name=csv_file_name,
                                   organism_name=organism_name,
                                   mast_file_name=mast_file_name,
                                   user_input=user_input)

from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator
from scipy.stats import spearmanr
import numpy as np


def analyze_intergenic(organism_name, cai_scores, evalues):
    if not cai_scores:
        logger.error(F"No cai values found for: {organism_name}")
        return

    percentile = 30

    cai_scores_array = np.array(cai_scores)
    top_threshold = np.percentile(cai_scores_array, percentile)
    bottom_threshold = np.percentile(cai_scores_array, 100-percentile)

    combined = [(cai_scores[i], evalues[i]) for i in range(len(cai_scores))]
    highly_expressed = [x for x in combined if x[0] >= top_threshold]
    lowly_expressed = [x for x in combined if x[0] < bottom_threshold]

    highly_coef, highly_p = spearmanr([x[0] for x in highly_expressed], [x[1] for x in highly_expressed])
    logger.info(F"Spearmans correlation coefficient for highly: {highly_coef}")
    lowly_coef, lowly_p = spearmanr([x[0] for x in lowly_expressed], [x[1] for x in lowly_expressed])
    logger.info(F"Spearmans correlation coefficient for lowly: {lowly_coef}")

    fig = pyplot.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.yaxis.set_major_locator(MultipleLocator(15))
    pyplot.scatter([x[0] for x in highly_expressed], [x[1] for x in highly_expressed],
                   label=F"highly expressed, rs: {highly_coef:.3f}, p-value: {highly_p:.3f}")
    pyplot.scatter([x[0] for x in lowly_expressed], [x[1] for x in lowly_expressed],
                   label=F"lowly expressed, rs: {lowly_coef:.3f}, p-value: {lowly_p:.3f}")
    pyplot.title(F"Intergenic motifs: {organism_name}")
    pyplot.xlabel("CAI score")
    pyplot.ylabel("e-value")
    pyplot.legend()

    base_directory = os.path.join(artifacts_directory, "promoters_not_for_user")
    plots_directory = os.path.join(base_directory, "plots")
    Path(plots_directory).mkdir(parents=True, exist_ok=True)
    plot_file_name = os.path.join(plots_directory, F"{organism_name}.png")
    pyplot.savefig(plot_file_name)
    pyplot.clf()


def analyze_selective(first_organism,
                      second_organism,
                      first_cai_scores,
                      first_evalues,
                      second_cai_scores,
                      second_evalues):
    if not first_cai_scores or not second_cai_scores:
        logger.error(F"No cai values found for: {first_organism} or {second_organism}")
        return

    percentile = 20

    first_cai_scores_array = np.array(first_cai_scores)
    first_top_threshold = np.percentile(first_cai_scores_array, percentile)
    first_bottom_threshold = np.percentile(first_cai_scores_array, 100 - percentile)

    second_cai_scores_array = np.array(second_cai_scores)
    second_top_threshold = np.percentile(second_cai_scores_array, percentile)
    second_bottom_threshold = np.percentile(second_cai_scores_array, 100 - percentile)

    combined_first = [(first_cai_scores[i], first_evalues[i]) for i in range(len(first_cai_scores))]
    combined_second = [(second_cai_scores[i], second_evalues[i]) for i in range(len(second_cai_scores))]
    highly_expressed_first = [x for x in combined_first if x[0] >= first_top_threshold]
    highly_expressed_second = [x for x in combined_second if x[0] >= second_top_threshold]


    coef_first, p_first = spearmanr([x[0] for x in highly_expressed_first], [x[1] for x in highly_expressed_first])
    coef_second, p_second = spearmanr([x[0] for x in highly_expressed_second], [x[1] for x in highly_expressed_second])
    logger.info(F"Spearmans correlation coefficient: {coef_first}, {coef_second}")

    fig = pyplot.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.yaxis.set_major_locator(MultipleLocator(15))
    pyplot.title(F"Optimized: {first_organism} Deoptimized: {second_organism}")
    pyplot.scatter([x[0] for x in highly_expressed_first], [x[1] for x in highly_expressed_first],
                   label=F"optimized rs: {coef_first:.3f}, p-value: {p_first:.3f}")
    pyplot.scatter([x[0] for x in highly_expressed_second], [x[1] for x in highly_expressed_second],
                   label=F"deoptimized rs: {coef_second:.3f}, p-value: {p_second:.3f}")
    pyplot.xlabel("CAI score")
    pyplot.ylabel("e-value")
    pyplot.legend()

    base_directory = os.path.join(artifacts_directory, "promoters_not_for_user")
    plots_directory = os.path.join(base_directory, "plots")
    Path(plots_directory).mkdir(parents=True, exist_ok=True)
    plot_file_name = os.path.join(plots_directory, F"opt_{first_organism}_deopt_{second_organism}.png")
    pyplot.savefig(plot_file_name)
    pyplot.clf()


def run_modules(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None,
                model_preferences_dict: typing.Optional[typing.Dict[str, str]] = None):
    user_inp_raw = user_input_dict or default_user_inp_raw

    model_preferences = models.ModelPreferences.init_from_dictionary(
        model_preferences_dict
    ) if model_preferences_dict is not None else models.ModelPreferences.init_from_config()

    try:
        import json

        # input_dict = user_IO.UserInputModule.run_module(user_inp_raw)   # keys: sequence, selected_prom, organisms

        # Store parsed input as json file
        # with open("parsed_input.json", "w") as user_input_file:
        #     json.dump(input_dict, user_input_file)

        # Read input from file
        with open("parsed_input.json", "r") as user_input_file:
            input_dict = json.load(user_input_file)

        # intergenic promoters
        # for organism_name in input_dict["organisms"].keys():
        #     cai_scores, evalues = create_csv_for_organism_intergenic(organism_name, input_dict)
        #     analyze_intergenic(organism_name, cai_scores, evalues)
        # logger.info("Intergenic end!")
        # exit(0)

        # selective promoters
        for organism_name in input_dict["organisms"].keys():
            for deopt_org_name in input_dict["organisms"].keys():
                if organism_name == deopt_org_name:
                    logger.info(F"Skip run for organism {organism_name} and deopt organism {deopt_org_name}")
                    continue

                first_cai_scores, first_evalues = create_csv_for_organism_selective(organism_name,
                                                                                    deopt_org_name,
                                                                                    input_dict,
                                                                                    is_first=True)

                second_cai_scores, second_evalues = create_csv_for_organism_selective(organism_name,
                                                                                      deopt_org_name,
                                                                                      input_dict,
                                                                                      is_first=False)
                analyze_selective(first_organism=organism_name,
                                  second_organism=deopt_org_name,
                                  first_cai_scores=first_cai_scores,
                                  first_evalues=first_evalues,
                                  second_cai_scores=second_cai_scores,
                                  second_evalues=second_evalues)

        logger.info("The end!")
        exit(0)

        ### unit 1 ############################################
        if model_preferences.restriction_enzymes or model_preferences.translation:
            final_cds, optimization_index, weakest_score = unit1(input_dict, model_preferences)
        else:
            final_cds = None
            optimization_index = None
            weakest_score = None
        #########################################################

        # ### unit 2 ############################################
        if model_preferences.transcription:
            p_name, native_prom, synth_promoter, evalue = promoters.promoterModule.run_module(input_dict)
        else:
            p_name = None
            native_prom = None
            synth_promoter = None
            evalue = None
        # #######################################################

        # TODO - get zip_directory from the user
        final_output, zip_file_path = user_IO.UserOutputModule.run_module(cds_sequence=final_cds,
                                                                          zscore=optimization_index,
                                                                          weakest_score=weakest_score,
                                                                          p_name=p_name,
                                                                          native_prom=native_prom,
                                                                          synth_promoter=synth_promoter,
                                                                          evalue=evalue,
                                                                          zip_directory=str(artifacts_directory))
        logger.info("Final output: %s, zip_file_path: %s", final_output, zip_file_path)
    except Exception as e:
        exception_str = traceback.format_exc()
        final_output = {
            'error_message': exception_str,
        }
        zip_file_path = None
        logger.error("Encountered unknown error when running modules. Error message: %s", exception_str)

    return final_output, zip_file_path


def unit1(input_dict, model_preferences: models.ModelPreferences):
    if model_preferences.translation:
        optimization_func = model_preferences.translation_function
        try:    # both CAI and tAI, select the one with the best optimization index tai optimization
            logger.info('\ntAI information:')
            cds_nt_final_tai = ORF.ORFModule.run_module(input_dict, 'tai', optimization_func)
            if model_preferences.restriction_enzymes:
                cds_nt_final_tai = RE.REModule.run_module(input_dict, cds_nt_final_tai)
            tai_mean_opt_index, tai_mean_deopt_index, tai_optimization_index, tai_weakest_score = \
                Zscore_calculation.ZscoreModule.run_module(cds_nt_final_tai, input_dict, optimization_type='tai')

            logger.info(f'Sequence:\n{cds_nt_final_tai}')
            logger.info(f'Optimized sequences score: {tai_mean_opt_index}, deoptimized sequence score: {tai_mean_deopt_index}')
            logger.info(f'Final optimization score: {tai_optimization_index}')

            # cai optimization
            logger.info('\nCAI information:')
            cds_nt_final_cai = ORF.ORFModule.run_module(input_dict, 'cai', optimization_func)
            if model_preferences.restriction_enzymes:
                cds_nt_final_cai = RE.REModule.run_module(input_dict, cds_nt_final_cai)  # todo: run both of them together to save time, or split creation of enzyme dict and the actual optimization (seems like a better solution)
            cai_mean_opt_index, cai_mean_deopt_index, cai_optimization_index , cai_weakest_score = \
                Zscore_calculation.ZscoreModule.run_module(cds_nt_final_cai, input_dict, optimization_type='cai')

            logger.info(f'Sequence:\n{cds_nt_final_cai}')
            logger.info(f'Optimized sequences score: {cai_mean_opt_index}, deoptimized sequence score: {cai_mean_deopt_index}')
            logger.info(f'Final optimization score: {cai_optimization_index}')

            if cai_optimization_index>tai_optimization_index:
                logger.info('CAI sequence was selected')
                final_cds = cds_nt_final_cai
                optimization_index = cai_optimization_index
                mean_opt_index = cai_mean_opt_index
                mean_deopt_index = cai_mean_deopt_index
                weakest_score = cai_weakest_score
            else:
                logger.info('tAI sequence was selected')
                final_cds = cds_nt_final_tai
                optimization_index = tai_optimization_index
                mean_opt_index = tai_mean_opt_index
                mean_deopt_index = tai_mean_deopt_index
                weakest_score = tai_weakest_score

        except:
            logger.info('\nCAI information:')
            final_cds = ORF.ORFModule.run_module(input_dict, 'cai', optimization_type=optimization_func)
            if model_preferences.restriction_enzymes:
                final_cds = RE.REModule.run_module(input_dict, final_cds)
            mean_opt_index, mean_deopt_index, optimization_index, weakest_score =\
                Zscore_calculation.ZscoreModule.run_module(final_cds, input_dict, 'cai')

    else:
        final_cds = RE.REModule.run_module(input_dict, input_dict['sequence'])
        mean_opt_index, mean_deopt_index, optimization_index, weakest_score = \
            Zscore_calculation.ZscoreModule.run_module(final_cds, input_dict, 'cai')

    logger.info(f'Sequence:\n{final_cds}')
    logger.info(f'Optimized sequences score: {mean_opt_index}, deoptimized sequence score: {mean_deopt_index}')
    logger.info(f'Weakest link score: {weakest_score}')
    logger.info(f'Final optimization score: {optimization_index}')
    return final_cds, optimization_index, weakest_score


if __name__ == "__main__":
    tic = time.time()
    run_modules()
    toc = time.time()
    modules_run_time = toc - tic
    print('Total modules run time: ', modules_run_time)
