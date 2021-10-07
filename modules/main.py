import os
import time
from pathlib import Path

from modules.logger_factory import LoggerFactory
from modules import Zscore_calculation, user_IO, RE, ORF, promoters

tic = time.time()
logger = LoggerFactory.create_logger("main")

current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")
user_inp_raw = {
    'sequence': os.path.join(base_path, 'mCherry_original.fasta'),
    'selected_promoters': None,
    'tuning_param':0.75,
    'organisms': {
                    # 'opt1': {'genome_path': os.path.join(base_path, 'Escherichia coli.gb'),
                    #          'optimized': True,
                    #          'expression_csv': os.path.join(base_path, 'ecoli_mrna_level.csv')},
                    #
                    # 'deopt1': {'genome_path': os.path.join(base_path, 'Bacillus subtilis.gb'),
                    #            'optimized': False,
                    #            'expression_csv': os.path.join(base_path, 'bacillus_mrna_level.csv')},

                    'deopt2': {'genome_path': os.path.join(base_path, 'Sulfolobus acidocaldarius.gb'),
                              'optimized': False,
                              'expression_csv': None},

                    'opt2': {'genome_path': os.path.join(base_path, 'Mycobacterium tuberculosis.gb'),
                             'optimized': True,
                             'expression_csv': None},
                    #
#                     'opt3': {'genome_path': os.path.join(base_path, 'Pantoea ananatis.gb'),
#                              'optimized': True,
#                              'expression_csv': None},

#                     'opt4': {'genome_path': os.path.join(base_path, 'Azospirillum brasilense.gb'),
#                              'optimized': True,
#                              'expression_csv': None}
            }
    }

model_preferences = {'RE': True, #todo: test restcition enzymes
                     'translation': True,
                     'transcription': False,
                     'translation_function': 'zscore_hill_climbing_average'#, 'single_codon_global', 'single_codon_local’, 'zscore_hill_climbing_average', 'zscore_hill_climbing_weakest_link'
}

def run_modules(user_inp_raw, model_preferences = model_preferences):
    try:
        input_dict = user_IO.UserInputModule.run_module(user_inp_raw) #keys: sequence, selected_prom, organisms

        ### unit 1 ############################################
        if model_preferences['RE'] or model_preferences['translation']:
            final_cds, optimization_index, weakest_score= unit1(input_dict, model_preferences)
        else:
            final_cds= None
            optimization_index= None
            weakest_score= None
        #########################################################

        # ### unit 2 ############################################
        if model_preferences['transcription']:
             p_name, native_prom, synth_promoter, evalue = promoters.promoterModule.run_module(input_dict)
        else:
            p_name = None
            native_prom = None
            synth_promoter = None
            evalue = None
        # #######################################################

        # TODO - get zip_directory from the user
        zip_directory_path = os.path.join(str(Path(__file__).parent.resolve()), "artifacts")
        Path(zip_directory_path).mkdir(parents=True, exist_ok=True)
        final_output = user_IO.UserOutputModule.run_module(cds_sequence=final_cds,
                                                           zscore=optimization_index,
                                                           weakest_score = weakest_score,
                                                           p_name=p_name,
                                                           native_prom = native_prom,
                                                           synth_promoter = synth_promoter,
                                                           evalue = evalue,
                                                           zip_directory=zip_directory_path)
        logger.info("Final output: %s", final_output)
    except Exception as e:
        final_output = {
            'error_message': str(e),
        }

    return final_output


def unit1(input_dict, model_preferences ):

    if model_preferences['translation']:
        optimization_func = model_preferences['translation_function']
        try: #both CAI and tAI, select the one with the best optimization index
            #tai optimization
            logger.info('\ntAI information:')
            cds_nt_final_tai = ORF.ORFModule.run_module(input_dict, 'tai', optimization_func)
            if model_preferences['RE']:
                cds_nt_final_tai = RE.REModule.run_module(input_dict, cds_nt_final_tai)
            tai_mean_opt_index, tai_mean_deopt_index, tai_optimization_index, weakest_score = \
                Zscore_calculation.ZscoreModule.run_module(cds_nt_final_tai, input_dict, optimization_type='tai')

            logger.info(f'Sequence:\n{cds_nt_final_tai}')
            logger.info(f'Optimized sequences score: {tai_mean_opt_index}, deoptimized sequence score: {tai_mean_deopt_index}')
            logger.info(f'Final optimization score: {tai_optimization_index}')

            #cai optimization
            logger.info('\nCAI information:')
            cds_nt_final_cai = ORF.ORFModule.run_module(input_dict, 'cai', optimization_func)
            if model_preferences['RE']:
                cds_nt_final_cai = RE.REModule.run_module(input_dict, cds_nt_final_cai)  # todo: run both of them together to save time, or split creation of enzyme dict and the actual optimization (seems like a better solution)
            cai_mean_opt_index, cai_mean_deopt_index, cai_optimization_index , weakest_score=\
                Zscore_calculation.ZscoreModule.run_module(cds_nt_final_cai, input_dict, optimization_type='cai')

            logger.info(f'Sequence:\n{cds_nt_final_cai}')
            logger.info(f'Optimized sequences score: {cai_mean_opt_index}, deoptimized sequence score: {cai_mean_deopt_index}')
            logger.info(f'Final optimization score: {cai_optimization_index}')

            if cai_optimization_index>tai_optimization_index:
                logger.info('CAI sequence was selected')
                final_cds = cds_nt_final_cai
                optimization_index = cai_optimization_index
            else:
                logger.info('tAI sequence was selected')
                final_cds = cds_nt_final_tai
                optimization_index = tai_optimization_index


        except:
            logger.info('\nCAI information:')
            final_cds = ORF.ORFModule.run_module(input_dict, 'cai', optimization_type=optimization_func)
            if model_preferences['RE']:
                final_cds = RE.REModule.run_module(input_dict, final_cds)
            mean_opt_index, mean_deopt_index, optimization_index, weakest_score =\
                Zscore_calculation.ZscoreModule.run_module(final_cds, input_dict, 'cai')


    else:
        final_cds = RE.REModule.run_module(input_dict, input_dict['sequence'])
        mean_opt_index, mean_deopt_index, optimization_index ,weakest_score = \
            Zscore_calculation.ZscoreModule.run_module(final_cds, input_dict, 'cai')

    logger.info(f'Sequence:\n{final_cds}')
    logger.info(f'Optimized sequences score: {mean_opt_index}, deoptimized sequence score: {mean_deopt_index}')
    logger.info(f'Weakest link score: {weakest_score}')
    logger.info(f'Final optimization score: {optimization_index}')
    return final_cds, optimization_index, weakest_score


if __name__ == "__main__":
    run_modules(user_inp_raw)
toc = time.time()

print('time: ', toc-tic)


