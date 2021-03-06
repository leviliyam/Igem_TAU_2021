import os
import shutil
import time
import traceback
from pathlib import Path
import typing

from logger_factory.logger_factory import LoggerFactory
from modules.testing_for_modules import generate_testing_data


# Create clean artifacts directory
artifacts_directory = Path(os.path.join(str(Path(__file__).parent.resolve()), "artifacts"))
if artifacts_directory.exists() and artifacts_directory.is_dir():
    shutil.rmtree(artifacts_directory)
artifacts_directory.mkdir(parents=True, exist_ok=True)
from modules import user_IO, ORF, sequence_family
from modules.stats.evaluation import ZscoreModule
from modules import models

current_directory = Path(__file__).parent.resolve()
base_path = os.path.join(Path(current_directory).parent.resolve(), "example_data")

logger = LoggerFactory.get_logger()


def run_modules(user_input_dict: typing.Optional[typing.Dict[str, typing.Any]] = None):
    user_inp_raw = user_input_dict or default_user_inp_raw

    try:
        before_parsing_input = time.time()
        user_input = user_IO.UserInputModule.run_module(user_inp_raw)
        after_parsing_input = time.time()

        logger.info(F"Total input processing time: {after_parsing_input-before_parsing_input}")

        # ####################################### family of sequences #####################################
        # in this part, user input is split into different inputs according to the sequence family theory
        clustered_user_inputs = sequence_family.SequenceFamilyModule.run_module(user_input)
        final_cds = None
        optimization_index = None
        weakest_score = None
        for input_cluster in clustered_user_inputs:
            # TODO - what do we want to display for each run? We should store the results differently
            final_cds, optimization_index, weakest_score = unit1(input_cluster)
        ##################################################################################################
        zip_directory = user_input.zip_directory or str(artifacts_directory)
        final_output = user_IO.UserOutputModule.run_module(cds_sequence=final_cds,
                                                           zscore=optimization_index,
                                                           weakest_score=weakest_score,
                                                           zip_directory=zip_directory)
        logger.info("Final output: %s", final_output)
    except Exception as e:
        exception_str = traceback.format_exc()
        final_output = {
            "error_message": exception_str,
        }
        logger.error("Encountered unknown error when running modules. Error message: %s", exception_str)

    return final_output


def unit1(user_input: models.UserInput):
    optimization_method = user_input.optimization_method
    try:
        # both CAI and tAI, select the one with the best optimization index tai optimization
        logger.info('tAI information:')
        cds_nt_final_tai = ORF.ORFModule.run_module(user_input, 'tai', optimization_method)

        tai_mean_opt_index, tai_mean_deopt_index, tai_optimization_index, tai_weakest_score = \
            ZscoreModule.run_module(cds_nt_final_tai, user_input, optimization_type='tai')

        logger.info(f'Sequence:\n{cds_nt_final_tai}')
        logger.info(f'Optimized sequences score: {tai_mean_opt_index}, '
                    f'deoptimized sequence score: {tai_mean_deopt_index}')
        logger.info(f'Final optimization score: {tai_optimization_index}')

        # cai optimization
        logger.info('CAI information:')
        cds_nt_final_cai = ORF.ORFModule.run_module(user_input, 'cai', optimization_method)

        cai_mean_opt_index, cai_mean_deopt_index, cai_optimization_index, cai_weakest_score = \
            ZscoreModule.run_module(cds_nt_final_cai, user_input, optimization_type='cai')

        logger.info(f'Sequence:\n{cds_nt_final_cai}')
        logger.info(f'Optimized sequences score: {cai_mean_opt_index}, '
                    f'deoptimized sequence score: {cai_mean_deopt_index}')
        logger.info(f'Final optimization score: {cai_optimization_index}')

        if cai_optimization_index > tai_optimization_index:
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
        logger.info('CAI information:')
        final_cds = ORF.ORFModule.run_module(user_input, 'cai', optimization_method=optimization_method)
        mean_opt_index, mean_deopt_index, optimization_index, weakest_score =\
            ZscoreModule.run_module(final_cds, user_input, 'cai')

    logger.info(f'Sequence:\n{final_cds}')
    logger.info(f'Optimized sequences score: {mean_opt_index}, deoptimized sequence score: {mean_deopt_index}')
    logger.info(f'Weakest link score: {weakest_score}')
    logger.info(f'Final optimization score: {optimization_index}')
    return final_cds, optimization_index, weakest_score


if __name__ == "__main__":
    tic = time.time()
    default_user_inp_raw = generate_testing_data(n_organisms=4,
                                                 percent_optimized=0.7,
                                                 clusters_count=1,
                                                 tuning_param=0.5)
    run_modules()
    toc = time.time()
    modules_run_time = toc - tic
    logger.info("Total modules run time: ", modules_run_time)
