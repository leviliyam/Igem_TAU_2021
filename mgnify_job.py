import csv
import json
import os
import requests
import shutil

from pathlib import Path

from Bio import SeqIO
from BCBio import GFF

supported_formats = ("FASTA", "GFF")


def download_genome_files(genome, download_url):
    download_folder = os.path.join("mgnify_files", genome)
    Path(download_folder).mkdir(parents=True, exist_ok=True)
    downloads = requests.get(download_url)
    downloads_json = downloads.json()
    files_list = downloads_json["data"]
    for file in files_list:
        file_format = file["attributes"]["file-format"]["name"]
        if file_format not in supported_formats:
            continue
        file_url = file["links"]["self"]
        file_response = requests.get(file_url, allow_redirects=True)
        file_name = file["id"]
        with open(os.path.join(download_folder, file_name), "wb") as fd:
            fd.write(file_response.content)


def download_files():
    url = "https://www.ebi.ac.uk/metagenomics/api/v1/genomes"
    page_count = 0
    while url is not None:
        print(F"Starting page no. {page_count}")
        genomes_response = requests.get(url)
        genomes_json = genomes_response.json()

        genomes_list = genomes_json["data"]
        for genome in genomes_list:
            print(F"Iterating id {genome['id']}")
            genome_type = genome["attributes"]["type"]
            if genome_type != "MAG":
                print(F"Skip id {genome['id']}. Type: {genome_type}")
                continue
            download_link = genome["relationships"]["downloads"]["links"]["related"]
            download_genome_files(genome=genome["id"], download_url=download_link)

        page_count += 1

        url = genomes_json["links"]["next"]


def group_catalogue_genomes(catalogue_genomes_url, catalogue_folder):
    while catalogue_genomes_url is not None:
        catalogue_genomes_response = requests.get(catalogue_genomes_url)
        catalogue_genomes = catalogue_genomes_response.json()
        genomes_list = catalogue_genomes["data"]
        for genome in genomes_list:
            genome_id = genome["id"]
            genome_folder = os.path.join("mgnify_files", genome_id)
            try:
                shutil.move(genome_folder, catalogue_folder)
            except Exception as e:
                print(e)
                print(F"Skipping genome {genome_id}")
            new_genome_folder = os.path.join(catalogue_folder, genome_id)
            Path(new_genome_folder).mkdir(parents=True, exist_ok=True)
            genome_summary_path = os.path.join(new_genome_folder, F"{genome_id}_summary.json")
            with open(genome_summary_path, "w") as genome_summary_file:
                json.dump(genome, genome_summary_file)
        catalogue_genomes_url = catalogue_genomes["links"]["next"]


def group_by_catalogues():
    url = "https://www.ebi.ac.uk/metagenomics/api/v1/genome-catalogues"
    catalogues_response = requests.get(url)
    catalogues_json = catalogues_response.json()
    for catalogue in catalogues_json["data"]:
        catalogue_name = catalogue["id"]
        catalogue_folder = os.path.join("mgnify_files", catalogue_name)
        Path(catalogue_folder).mkdir(parents=True, exist_ok=True)
        genomes_url = catalogue["relationships"]["genomes"]["links"]["related"]
        group_catalogue_genomes(genomes_url, catalogue_folder)


def filter_files():
    catalogues_list = Path("mgnify_files").glob("*")
    for catalogue_path in catalogues_list:
        selected_catalogue_files = []
        completeness_values = []
        contamination_values = []
        genomes_list = Path(catalogue_path).glob("*")
        for genome in genomes_list:
            summary_files = Path(genome).glob("*.json")
            for summary_file in summary_files:
                with open(str(summary_file)) as file:
                    genome_summary = json.load(file)
                    completeness_values.append(genome_summary["attributes"]["completeness"])
                    contamination_values.append(genome_summary["attributes"]["contamination"])
                    # if genome_summary["attributes"]["completeness"] < 99:
                    #     print(F"Skipping {str(genome)}, completeness is: {genome_summary['attributes']['completeness']}")
                    #     continue
                    # if genome_summary["attributes"]["contamination"] > 0.5:
                    #     print(F"Skipping {str(genome)}, contamination is: {genome_summary['attributes']['contamination']}")
                    #     continue
                    # selected_catalogue_files.append(str(genome))
        selected_files_path = os.path.join(str(catalogue_path), "selected_genome_files.txt")
        histogram_path = os.path.join(str(catalogue_path), "histogram.csv")
        with open(histogram_path, "w") as histogram_file:
            writer = csv.writer(histogram_file)
            writer.writerow(["completeness", "contamination"])
            for i in range(len(completeness_values)):
                writer.writerow([completeness_values[i], contamination_values[i]])
        # with open(selected_files_path, "w") as selected_genome_file:
        #     selected_genome_file.writelines(selected_catalogue_files)


# write ideas for the promoter model
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse


def extract_prom(cds_start, cds_stop, cds_strand, cds_names, prom_length, genome):
    prom_dict = {}
    genome_fwd = genome
    genome_rev = genome
    for idx, vals in enumerate(zip(cds_start, cds_stop, cds_strand)):
        start, end, strand = vals
        if strand == 1:
            genome_fwd = genome[:start] + '-' * (end - start) + genome[end:]
        else:
            genome_rev = genome[:start] + '-' * (end - start) + genome[end:]
    for idx, vals in enumerate(zip(cds_start, cds_strand, cds_names)):
        start, strand, name = vals
        if strand == 1:
            prom = genome_fwd[start - prom_length:start]
        else:
            prom = reverse_complement(genome_rev[start:start + prom_length])
        # TODO - should we also exclude overlaps with cds?
        if prom != '-' * prom_length and len(prom.replace('-', '')) > 0:
            prom_dict[name] = prom.replace('-', '')
    return prom_dict


def extract_intergenic(cds_start, cds_stop, cds_strand, prom_length, genome, len_th=30):
    """
    extract intergenic sequences for streme intergenic motifs
    :param cds_start:
    :param cds_stop:
    :param cds_strand: -1 if the sequence is on the reverse strand, 1 for forward
    :param prom_length: length of the sequences considered as promoters
    :param genome: string of the whole 5'->3' genome on the fwd strand
    :param len_th:shortest sequence to add to the file
    :return: list of intergenic sequences
    """
    for idx, vals in enumerate(zip(cds_start, cds_stop, cds_strand)):
        start, end, strand = vals
        if strand == 1:
            genome = genome[:start - prom_length] + '-' * (end - start + prom_length) + genome[end:]
        else:  # strand ==-1
            genome = genome[:start] + '-' * (end - start + prom_length) + genome[end + prom_length:]
    if not genome:
        intergenic_list = []
    else:
        intergenic_list = genome.split('-')
    return {k: i for k, i in enumerate(intergenic_list) if len(i) > len_th}


def write_fasta(file_path, list_seq, list_name):
    with open(file_path, "w") as ofile:
        for i in range(len(list_seq)):
            ofile.write(">" + str(list_name[i]) + "\n" + str(list_seq[i]) + "\n")


def read_fasta(file_path):
    with open(file_path) as file:
        return {record.description: str(record.seq) for record in SeqIO.parse(file, "fasta")}


def extract_promoters_data_for_file(file_path, sequece_fasta_file, directory_path="."):
    fasta_dict = read_fasta(sequece_fasta_file)
    fasta_keys = [x.split("-") for x in fasta_dict.keys()]

    gff_file = GFF.parse(file_path)
    final_prom_dict = {}
    final_intergenic_dict = {}
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
            print(F"Failed to find a matching genome sequence for: {record.id} in {file_path}. Skipping all record"
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
        prom200_dict = extract_prom(starts, ends, strands, name_and_function, prom_length=200, genome=genome)
        final_prom_dict.update(prom200_dict)
        intergenic_dict = extract_intergenic(starts, ends, strands, prom_length=200, genome=genome, len_th=30)
        final_intergenic_dict.update(intergenic_dict)

    prom_file_path = os.path.join(directory_path, "promoters.fasta")
    write_fasta(prom_file_path, list(final_prom_dict.values()), list(final_prom_dict.keys()))
    intergenic_file_path = os.path.join(directory_path, "intergenic.fasta")
    write_fasta(intergenic_file_path, list(final_intergenic_dict.values()), list(final_intergenic_dict.keys()))


def extract_promoters_data():
    catalogues_list = Path("mgnify_files").glob("*")
    for catalogue_path in catalogues_list:
        if "human-oral-v1-0" not in str(catalogue_path):
            continue
        genomes_list = Path(catalogue_path).glob("*")
        for genome in genomes_list:
            gff_files = Path(genome).glob("*.gff")
            for gff_file in gff_files:
                fasta_files = list(Path(genome).glob(F"{Path(genome).name}.fna"))
                if len(fasta_files) != 1:
                    print(F"Skipping {genome}. Wrong number of .fna files.")
                extract_promoters_data_for_file(str(gff_file), str(fasta_files[0]))


if __name__ == "__main__":
    print("Start")
    # extract_promoters_data()
    extract_promoters_data_for_file(
         file_path=r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\mgnify_files\MGYG000290000\MGYG000290000.gff",
         # file_path=r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\mgnify_files\MGYG000295586\MGYG000295586.gff",
         sequece_fasta_file=r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\mgnify_files\MGYG000290000\MGYG000290000.fna",
         # r"C:\Users\Kama\Documents\Moran\biomedical-engineering\microbiome-optimization\mgnify_files\MGYG000295586\MGYG000295586.fna"
    )
    print("The end")

