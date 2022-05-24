import os
import statistics
import time
import traceback
from pathlib import Path
import typing
import csv
import json
import xml.etree.ElementTree as et

from matplotlib import pyplot
import numpy as np

from modules.logger_factory import LoggerFactory

# Create clean artifacts directory
artifacts_directory = Path(os.path.join(str(Path(__file__).parent.resolve()), "artifacts"))
artifacts_directory.mkdir(parents=True, exist_ok=True)

from modules import Zscore_calculation, user_IO, RE, ORF, promoters
from modules import models

logger = LoggerFactory.create_logger("main")

current_directory = Path(__file__).parent.resolve()
# base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
base_path = "/arabidopsis_microbiome"
# base_path = "C:\\Users\\Kama\\Documents\\Moran\\biomedical-engineering\\microbiome-optimization\\arabidopsis_microbiome"

default_user_inp_raw = {
    'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
    'selected_promoters': None,
    'tuning_param': 0.5,
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
            # if "hypothetical" in seq_name:
            #     logger.info(F"skipping hypothetical match {seq_name}")
            #     continue
            promoter_score = m.find("score")
            promoter_evalue = promoter_score.get("evalue")
            matching_genes = [k for k in organism_info['cai_scores'].keys() if k.startswith(seq_name)]
            if len(matching_genes) > 1:
                logger.error(F"Found more than one matching gene for promoter! {matching_genes}")
                exit(1)
            matching_gene = matching_genes[0]
            cai_score = organism_info["cai_scores"][matching_gene]

            data = [matching_gene, cai_score, promoter_evalue]
            csv_writer.writerow(data)

            cai_scores.append(cai_score)
            evalues.append(promoter_evalue)

    return cai_scores, evalues


def create_csv_for_organism_intergenic(organism_name: str, user_input):
    csv_file_name = F"{organism_name}.csv"
    csv_directory_name = "intergenic_csv_files"
    mast_file_name = F"intergenic_mast\\motif streme seq {organism_name} 100 200\\mast.xml".replace(" ", "_")
    return create_csv_for_organism(csv_directory_name=csv_directory_name,
                                   csv_file_name=csv_file_name,
                                   organism_name=organism_name,
                                   mast_file_name=mast_file_name,
                                   user_input=user_input)


def extract_intergenic_region_evalues(organism_name: str):
    base_directory = os.path.join(artifacts_directory, "promoters_not_for_user")
    mast_file_name = F"mast_inter_{organism_name}\\mast.xml".replace(" ", "_")
    base_file = os.path.join(base_directory, mast_file_name)
    tree = et.parse(base_file)
    root = tree.getroot()

    evalues = []
    for m in root.findall(".//sequence"):
        seq_name = m.get("name")
        # if "hypothetical" in seq_name:
        #     logger.info(F"skipping hypothetical match {seq_name}")
        #     continue
        promoter_score = m.find("score")
        promoter_evalue = promoter_score.get("evalue")

        evalues.append(promoter_evalue)

    return evalues

def analyze_intergenic(organism_name, cai_scores, evalues, inter_e_values):
    # if not cai_scores:
    #     logger.error(F"No cai values found for: {organism_name}")
    #     return

    base_directory = os.path.join(artifacts_directory, "promoters_not_for_user")
    plots_directory = os.path.join(base_directory, "plots")

    percentile = 20

    # cai_scores_array = np.array(cai_scores)
    # top_threshold = np.percentile(cai_scores_array, 100-percentile)
    # bottom_threshold = np.percentile(cai_scores_array, percentile)
    #
    # combined = [(cai_scores[i], float(evalues[i])) for i in range(len(cai_scores))]
    # combined_coef, combined_p = spearmanr([x[0] for x in combined], [x[1] for x in combined])
    # highly_expressed = [x for x in combined if x[0] >= top_threshold]
    # mean_highly_expression = statistics.mean([x[0] for x in highly_expressed])
    # mean_highly_e_value = statistics.mean([float(x[1]) for x in highly_expressed])
    # median_highly_e_value = statistics.median([float(x[1]) for x in highly_expressed])
    # std_highly_e_value = 0
    # if len(highly_expressed) > 1:
    #     std_highly_e_value = statistics.stdev([float(x[1]) for x in highly_expressed])
    #
    # lowly_expressed = [x for x in combined if x[0] < bottom_threshold]
    # mean_lowly_expression = statistics.mean([x[0] for x in lowly_expressed])
    # mean_lowly_e_value = statistics.mean([float(x[1]) for x in lowly_expressed])
    # median_lowly_e_value = statistics.median([float(x[1]) for x in lowly_expressed])
    # std_lowly_e_value = 0
    # if len(lowly_expressed) > 1:
    #     std_lowly_e_value = statistics.stdev([float(x[1]) for x in lowly_expressed])

    # --------------------------
    # Calcualte e-value histogram for promoters and intergenic
    # --------------------------
    float_evalues = [float(x) for x in evalues]
    if not evalues or not inter_e_values:
        logger.info(F"missing values for {organism_name}")
        return
    bottom_evalues_threshold = np.percentile(float_evalues, percentile)
    # float_evalues = [x for x in float_evalues if x < bottom_evalues_threshold]
    mean_evalues = statistics.mean(float_evalues)
    median_evalues = statistics.median(float_evalues)
    float_inter_e_values = [float(x) for x in inter_e_values]
    bottom_inter_evalues_threshold = np.percentile(float_inter_e_values, percentile)
    # float_inter_e_values = [x for x in float_inter_e_values if x < bottom_inter_evalues_threshold]
    mean_inter_evalues = statistics.mean(float_inter_e_values)
    median_inter_evalues = statistics.median(float_inter_e_values)

    from scipy.stats import mannwhitneyu
    w, p = mannwhitneyu(x=float_evalues, y=float_inter_e_values)

    ax = pyplot.subplot(121)
    pyplot.title(F"Promoters\n{organism_name}")
    pyplot.xlabel("e-value")
    pyplot.ylabel("# sequences")
    pyplot.hist(float_evalues, bins=20, color="tab:blue", edgecolor='k')
    pyplot.axvline(mean_evalues, color='k', linestyle='dashed', linewidth=2, label=F"mean: {mean_evalues:.3f}")
    pyplot.axvline(median_evalues, color='k', linestyle='dashdot', linewidth=1, label=F"median: {median_evalues:.3f}")
    pyplot.legend()
    pyplot.subplot(122, sharex=ax)
    pyplot.title(F"Intergenic regions\n{organism_name}")
    pyplot.xlabel("e-value")
    pyplot.ylabel("# sequences")
    pyplot.hist(float_inter_e_values, color="tab:red", bins=20, edgecolor='k')
    pyplot.axvline(mean_inter_evalues, color='k', linestyle='dashed', linewidth=2, label=F"mean: {mean_inter_evalues:.3f}")
    pyplot.axvline(median_inter_evalues, color='k', linestyle='dashdot', linewidth=1, label=F"median: {median_inter_evalues:.3f}")
    pyplot.legend()
    pyplot.tight_layout()
    directory = os.path.join(plots_directory, "e-value histogram")
    Path(directory).mkdir(parents=True, exist_ok=True)
    logger.info(F"{organism_name} P value: {p}")
    plot_file_name = os.path.join(directory, F"{organism_name}_p_value_{p}.png")
    pyplot.savefig(plot_file_name)
    pyplot.clf()

    return

    # --------------------------
    # Calcualte cai score histogram
    # --------------------------
    # pyplot.title(F"CAI score histogram: {organism_name}")
    # fig = pyplot.figure()
    # pyplot.hist([x[0] for x in highly_expressed], label=F"highly_expressed", color="b", bins=20)
    # pyplot.hist([x[0] for x in lowly_expressed], label=F"lowly_expressed", color="r", bins=20)
    # pyplot.hist([x[0] for x in combined if x not in highly_expressed and x not in lowly_expressed],
    #             label=F"highly_expressed", color="g", bins=20)
    # pyplot.legend()
    # pyplot.xlabel("cai score")
    # directory = os.path.join(plots_directory, "cai histogram")
    # Path(directory).mkdir(parents=True, exist_ok=True)
    # plot_file_name = os.path.join(directory, F"{organism_name}.png")
    # pyplot.savefig(plot_file_name)
    # pyplot.clf()

    # --------------------------
    # Mean/Median e-value for each group
    # --------------------------
    pyplot.title(F"Mean/Median e-values: {organism_name}")
    x = [mean_highly_expression, mean_lowly_expression]
    y = [mean_highly_e_value, mean_lowly_e_value]
    e = [std_highly_e_value, std_lowly_e_value]
    pyplot.errorbar(x, y, e, linestyle='dashed', marker='^', ecolor="red", capsize=10, label="mean")
    x = [mean_highly_expression+1, mean_lowly_expression+1]
    y1 = [median_highly_e_value, median_lowly_e_value]
    text = ["highly\nexpressed", "lowly\nexpressed"]
    pyplot.errorbar(x, y1, e, linestyle='dashed', marker='^', ecolor="red", capsize=10, label="median")
    pyplot.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    for i in range(len(x)):
        pyplot.annotate(text[i], (x[i]-1, y[i]))
        pyplot.annotate(text[i], (x[i], y1[i]))
    pyplot.ylabel("Mean/Median e-value")
    pyplot.legend()
    directory = os.path.join(plots_directory, "intergenic_median")
    Path(directory).mkdir(parents=True, exist_ok=True)
    plot_file_name = os.path.join(directory, F"{organism_name}.png")
    pyplot.savefig(plot_file_name)
    pyplot.clf()

    # --------------------------
    # Plot cai_score as function of e-value
    # --------------------------
    pyplot.title(F"CAI score by e-value: {organism_name}")
    e_value_array = np.array([float(x) for x in evalues])
    top_threshold_e = np.percentile(e_value_array, 100 - percentile)
    bottom_threshold_e = np.percentile(e_value_array, percentile)
    high_e = [x for x in combined if x[1] >= top_threshold_e]
    low_e = [x for x in combined if x[1] < bottom_threshold_e]
    # fig = pyplot.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # ax.yaxis.set_major_locator(MultipleLocator(15))
    pyplot.scatter([x[1] for x in high_e], [x[0] for x in high_e],
                   label=F"high e-value")
    pyplot.scatter([x[1] for x in low_e], [x[0] for x in low_e],
                   label=F"low e-value")
    pyplot.title(F"Intergenic motifs: {organism_name}")
    pyplot.xlabel("e-value")
    pyplot.ylabel("CAI score")
    pyplot.legend()
    directory = os.path.join(plots_directory, "intergenic_cai_score_by_evalue")
    Path(directory).mkdir(parents=True, exist_ok=True)
    plot_file_name = os.path.join(directory, F"{organism_name}.png")
    pyplot.savefig(plot_file_name)
    pyplot.clf()

    # Plot mean/medican of cai score per e-value
    pyplot.title(F"Mean/Median CAI score: {organism_name}")
    mean_high_cai = statistics.mean([float(x[0]) for x in high_e])
    median_high_cai = statistics.median([float(x[0]) for x in high_e])
    mean_low_cai = statistics.mean([float(x[0]) for x in low_e])
    median_low_cai = statistics.median([float(x[0]) for x in low_e])
    std_high_e_value = 0
    if len(high_e) > 1:
        std_high_e_value = statistics.stdev([float(x[0]) for x in high_e])
    std_low_e_value = 0
    if len(low_e) > 1:
        std_low_e_value = statistics.stdev([float(x[0]) for x in low_e])
    x = [1, 2]
    y = [mean_high_cai, mean_low_cai]
    e = [std_high_e_value, std_low_e_value]
    pyplot.errorbar(x, y, e, linestyle='dashed', marker='^', ecolor="red", capsize=10, label="mean")
    x1 = [3, 4]
    y1 = [median_high_cai, median_low_cai]
    text = ["high\ne-value", "low\ne-value"]
    pyplot.errorbar(x1, y1, e, linestyle='dashed', marker='^', ecolor="red", capsize=10, label="median")
    pyplot.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False,  # ticks along the top edge are off
        labelbottom=False)  # labels along the bottom edge are off
    for i in range(len(x)):
        pyplot.annotate(text[i], (x[i], y[i]))
        pyplot.annotate(text[i], (x1[i], y1[i]))
    pyplot.ylabel("Mean/Median cai_score")
    pyplot.legend()
    directory = os.path.join(plots_directory, "intergenic_e_value_mean")
    Path(directory).mkdir(parents=True, exist_ok=True)
    plot_file_name = os.path.join(directory, F"{organism_name}.png")
    pyplot.savefig(plot_file_name)
    pyplot.clf()

    # --------------------------
    # Histograms
    # --------------------------
    pyplot.title(F"e-value distribution: {organism_name}")
    ax = pyplot.subplot(311)
    pyplot.ylabel("# promoters")
    pyplot.hist([x[1] for x in combined], label=F"all", color="b", bins=20)
    pyplot.legend()
    ax2 = pyplot.subplot(312, sharex=ax)
    pyplot.hist([x[1] for x in highly_expressed], label=F"highly_expressed", color="g", bins=20)
    pyplot.legend()
    pyplot.subplot(313, sharex=ax, sharey=ax2)
    pyplot.hist([x[1] for x in lowly_expressed], label=F"lowly_expressed", color="r", bins=20)
    pyplot.legend()
    pyplot.xlabel("e-value")
    directory = os.path.join(plots_directory, "intergenic_histogram")
    Path(directory).mkdir(parents=True, exist_ok=True)
    plot_file_name = os.path.join(directory, F"{organism_name}.png")
    pyplot.savefig(plot_file_name)
    pyplot.clf()

    # --------------------------
    # Scatter plots
    # --------------------------
    # fig = pyplot.figure()
    # ax = fig.add_subplot(1, 1, 1)
    # ax.yaxis.set_major_locator(MultipleLocator(15))
    # pyplot.scatter([x[0] for x in highly_expressed], [x[1] for x in highly_expressed],
    #                label=F"highly expressed, rs: {highly_coef:.3f}, p-value: {highly_p:.3f}")
    # pyplot.scatter([x[0] for x in lowly_expressed], [x[1] for x in lowly_expressed],
    #                label=F"lowly expressed, rs: {lowly_coef:.3f}, p-value: {lowly_p:.3f}")
    # pyplot.title(F"Intergenic motifs: {organism_name}")
    # pyplot.xlabel("CAI score")
    # pyplot.ylabel("e-value")
    # pyplot.legend()

    # Path(plots_directory).mkdir(parents=True, exist_ok=True)
    # plot_file_name = os.path.join(plots_directory, F"{organism_name}.png")
    # pyplot.savefig(plot_file_name)
    # pyplot.clf()


def run_modules(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None,
                model_preferences_dict: typing.Optional[typing.Dict[str, str]] = None):
    # user_inp_raw = user_input_dict or default_user_inp_raw
    #
    # for organism_name in user_inp_raw["organisms"].keys():
    #     user_inp_raw["organisms"][organism_name]["optimized"] = True

    model_preferences = models.ModelPreferences.init_from_dictionary(
        model_preferences_dict
    ) if model_preferences_dict is not None else models.ModelPreferences.init_from_config()

    try:
        # input_dict = user_IO.UserInputModule.run_module(user_inp_raw)   # keys: sequence, selected_prom, organisms

        # # Store parsed input as json file
        # with open("parsed_input.json", "w") as user_input_file:
        #     json.dump(input_dict, user_input_file)

        # exit(0)

        # Read input from file
        with open("parsed_input.json", "r") as user_input_file:
             input_dict = json.load(user_input_file)

        organisms_names = list(input_dict["organisms"].keys())
        organisms_count = len(organisms_names)
        wanted_count = 2
        unwanted_count = 2

        import random
        selected_orgs_indices = random.sample(range(organisms_count), wanted_count + unwanted_count)
        selected_orgs = {}
        for i in selected_orgs_indices:
            org = organisms_names[i]
            selected_orgs[org] = input_dict["organisms"][org]
            if i >= wanted_count:
                selected_orgs[org]["optimized"] = False

        input_dict["organisms"] = selected_orgs
        promoters.promoterModule.run_module(input_dict)

        exit(0)

        # intergenic promoters
        for organism_name in input_dict["organisms"].keys():
            cai_scores, evalues = create_csv_for_organism_intergenic(organism_name, input_dict)
            inter_e_values = extract_intergenic_region_evalues(organism_name)
            analyze_intergenic(organism_name, cai_scores, evalues, inter_e_values)
        logger.info("Intergenic end!")
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
