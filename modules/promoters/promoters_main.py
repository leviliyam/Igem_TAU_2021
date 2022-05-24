from modules.logger_factory import LoggerFactory
from modules.promoters.before_meme_run import *
from modules.promoters.run_meme import *
from modules.promoters.motif_filtering_after_streme import *
from modules.promoters.promoter_choosing_after_mast import *

#todo:
# run_mast should return the f_name as an output
# as to all IO: write everything into a directory called "promoters_not_for_user", except the mast html file which
# should be copied into the "logs" directory
# divide code into different sub-modules
# make run_module return the best sequence and it's e value.
# re-rank promoters

logger = LoggerFactory.create_logger("promoters")


class promoterModule(object):

    @staticmethod
    def run_module(full_input_dict):
        logger.info('##########################')
        logger.info('# PROMOTERS INFORMATION #')
        logger.info('##########################')

        # FIXME - extracting motifs
        # promoter_file_path = create_files_for_meme(full_input_dict)
        # run_streme()

        # FIXME - only for analysis (mast runs)
        # # Calculate mast for intergenic motifs
        # for intergenic_motif_file in find_all_inter_files():
        #     organism_name = "_".join(str(Path(intergenic_motif_file).parent.name).split("_")[:2])
        #     org_promoter_dir = os.path.join(deopt_path, organism_name)
        #     org_promoter_file_path = os.path.join(org_promoter_dir, F"{organism_name}_100_200.fasta")
        #     run_mast(intergenic_motif_file, org_promoter_file_path)
        #
        # # Calculate mast for selective motifs
        # for index, selective_motif_file in enumerate(find_all_selective_files()):
        #     splitted_motif_file_name = str(Path(selective_motif_file).parent.name).split("_")
        #     # mast for first org
        #     first_organism_name = "_".join(splitted_motif_file_name[:2])
        #     logger.info(F"Run mast for first organism: {first_organism_name}")
        #     first_org_promoter_dir = os.path.join(deopt_path, first_organism_name)
        #     first_org_promoter_file_path = os.path.join(first_org_promoter_dir, F"{first_organism_name}_100_200.fasta")
        #     run_mast(selective_motif_file, first_org_promoter_file_path)
        #     # mast for second org
        #     second_organism_name = "_".join(splitted_motif_file_name[-4:-2])
        #     logger.info(F"Run mast for first organism: {second_organism_name}")
        #     second_org_promoter_dir = os.path.join(deopt_path, second_organism_name)
        #     second_org_promoter_file_path = os.path.join(second_org_promoter_dir, F"{second_organism_name}_100_200.fasta")
        #     run_mast(selective_motif_file, second_org_promoter_file_path, "1")
        #
        # exit(0)

        motif_file_path = create_final_motif_xml(full_input_dict)
        logger.info("motif file path: %s", motif_file_path)

        promoter_file_path = create_unified_promoters_file(full_input_dict)
        mast_output_folder = run_mast(motif_file_path, promoter_file_path)
        logger.info("mast output folder: %s", mast_output_folder)
        # return modify_promoter(promoter_file_path, mast_output_folder)


