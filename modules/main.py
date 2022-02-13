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
base_path = "/arabidopsis_microbiome"
# base_path = "C:\\Users\\Kama\\Documents\\Moran\\biomedical-engineering\\microbiome-optimization\\arabidopsis_microbiome"

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


def run_modules(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None,
                model_preferences_dict: typing.Optional[typing.Dict[str, str]] = None):
    user_inp_raw = user_input_dict or default_user_inp_raw

    model_preferences = models.ModelPreferences.init_from_dictionary(
        model_preferences_dict
    ) if model_preferences_dict is not None else models.ModelPreferences.init_from_config()

    from modules.promoters.motif_filtering_after_streme import find_all_selective_files
    from modules.promoters.run_meme import run_mast
    from modules.promoters.globals_and_shared_methods import deopt_path

    for index, selective_motif_file in enumerate(find_all_selective_files()):
        splitted_motif_file_name = str(Path(selective_motif_file).parent.name).split("_")
        # mast for first org
        first_organism_name = "_".join(splitted_motif_file_name[:2])
        second_organism_name = "_".join(splitted_motif_file_name[-4:-2])
        logger.info(F"Run mast for first organism: {first_organism_name}")
        first_org_promoter_dir = os.path.join(deopt_path, first_organism_name)
        first_org_promoter_file_path = os.path.join(first_org_promoter_dir, F"{first_organism_name}_100_200.fasta")
        run_mast(selective_motif_file, first_org_promoter_file_path, F"{first_organism_name}_{second_organism_name}_first")
        # mast for second org
        logger.info(F"Run mast for second organism: {second_organism_name}")
        second_org_promoter_dir = os.path.join(deopt_path, second_organism_name)
        second_org_promoter_file_path = os.path.join(second_org_promoter_dir, F"{second_organism_name}_100_200.fasta")
        run_mast(selective_motif_file, second_org_promoter_file_path, F"{first_organism_name}_{second_organism_name}_second")

    logger.info("The end")
    exit(0)

    try:
        # TODO - convert the input dictionary to a model we can pass to the other models
        input_dict = user_IO.UserInputModule.run_module(user_inp_raw)   # keys: sequence, selected_prom, organisms

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
