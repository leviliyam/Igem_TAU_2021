# -*- coding: utf-8 -*-
import subprocess

destination_dir = '../../data/'

def run_cmd(cmd, verbose=False, *args, **kwargs):
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True
    )
    std_out, std_err = process.communicate()
    if verbose:
        print(std_out.strip(), std_err)


def download_files():
    with open("cds_files.txt", "r") as cds_file:
        file_names = cds_file.readlines()
        for file_name in file_names:
            cds_destination_file = destination_dir+'ncbi_genome_cds'
            command = 'wget' +  file_name.strip() + ' ‐P ' + cds_destination_file
            run_cmd(command)
            print(command)
            command = 'gzip -d' +  cds_destination_file
            run_cmd(command)

    with open("rna_files.txt", "r") as rna_files:
        file_names = rna_files.readlines()
        for file_name in file_names:
            rna_destination_file = destination_dir+'ncbi_genome_rna'
            command = 'wget ' +  file_name.strip() + ' ‐P ' + rna_destination_file
            run_cmd(command)
            print(command)
            command = 'gzip -d' +  rna_destination_file
            run_cmd(command)



if __name__ == "__main__":
    print("Start")
    download_files()
    print("The end")










# # -*- coding: utf-8 -*-
# import codecs
# import subprocess
#
# destination_dir = '../../data/genbank_tls'
# index_file = codecs.open("tls_index.html", 'r')
#
#
# def run_cmd(cmd, verbose=False, *args, **kwargs):
#     process = subprocess.Popen(
#         cmd,
#         stdout=subprocess.PIPE,
#         stderr=subprocess.PIPE,
#         shell=True
#     )
#     std_out, std_err = process.communicate()
#     if verbose:
#         print(std_out.strip(), std_err)
#
#
# def download_files():
#     Lines = index_file.readlines()
#
#
#     count = 0
#     for line in Lines:
#         count += 1
#         line = line.strip()
#         if "fsa_nt.gz" in line or ".mstr.gbff.gz" in line:
#             f_name = line.split('"')[1]
#             run_cmd("wget ‐P ncbi_tls https://ftp.ncbi.nlm.nih.gov/genbank/tls/K/" +f_name +" --no-check-certificate")
#             # print("wget ‐P " + destination_dir + " https://ftp.ncbi.nlm.nih.gov/genbank/tls/K/" +f_name +" --no-check-certificate")
#             destination_file = destination_dir + '/' + f_name
#             run_cmd("gzip -d ncbi_tls/" + f_name)
#
#
# if __name__ == "__main__":
#     print("Start")
#     download_files()
#     print("The end")
#
#
#
# # wget -P ../../data/genbank_tls
# # https://ftp.ncbi.nlm.nih.gov/genbank/tls/K/tls.KAAS.mstr.gbff.gz --no-check-certificate
# # wget ‐P.. /../ data / genbank_tls
# # https: // ftp.ncbi.nlm.nih.gov / genbank / tls / K / tls.KAAW.mstr.gbff.gz - -no - check - certificate &
