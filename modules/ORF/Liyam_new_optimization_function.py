from modules import models
from modules.stats.optimization import OptimizationModule
from modules.shared_functions_and_vars import nt_to_aa


def change_all_codons_of_aa(seq: str, selected_codon: str):
    '''
    change all synonymous codons in the str seq to the selected codon, return str as well
    '''
    split_seq = [seq[i:i+3].upper() for i in range(0, len(seq), 3)]
    new_split_seq = []
    for codon in split_seq:
        if nt_to_aa[codon] == nt_to_aa[selected_codon]:
            new_split_seq.append(selected_codon)
        else:
            new_split_seq.append(codon)
    return ''.join(new_split_seq)

import csv
import time

# in each round - check all single synonymous codon changes and calculate optimization score - take the best one
def hill_climbing_optimize_by_zscore(seq: str,
                                     user_input: models.UserInput,
                                     cai_or_tai: str,
                                     max_iter: int,
                                     optimization_method: models.OptimizationMethod):
    """
    hill climbing function for performing codon optimization
    in each iteration - for each codon, change all synonymous codons to a specific one and test the zscore of the new
    sequence after each iteration, select the sequence with the best zscore - if it was not changed since the last
    iteration, break. The maximum number of iterations allowed is "max_iter"
    @seq: str, tested seq
    @inp_dict: input dict after usr_inp code
    @opt_type: 'cai' or 'tai'
    @max_iter: maximal number of iterations to perform
    return: seq, which is the optimized sequence
    """
    seq_options = {}
    score = OptimizationModule.run_module(
        final_seq=seq,
        user_input=user_input,
        cai_or_tai=cai_or_tai,
        optimization_method=optimization_method,
    )
    seq_options[seq] = score

    # Single codon replacement
    with open('codon_change_by_iteration.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["iteration", "codon", "amino acid", "score"])
        for run in range(max_iter):
            # TODO - we change in each iteration a single codon. We may consider changing at most X codons at a time to
            #  reduce risk of falling to a local maxima.
            tested_seq_to_codon = {}
            for codon in nt_to_aa.keys():
                tested_seq = change_all_codons_of_aa(seq, codon)
                tested_seq_to_codon[tested_seq] = codon

                score = OptimizationModule.run_module(
                    final_seq=tested_seq,
                    user_input=user_input,
                    cai_or_tai=cai_or_tai,
                    optimization_method=optimization_method,
                )
                seq_options[tested_seq] = score
                # TODO - score, what codon what chosen and of which AA - changes as function of iteration

            new_seq = max(seq_options, key=seq_options.get)

            changed_codon = tested_seq_to_codon[new_seq]
            aa = nt_to_aa[changed_codon]
            writer.writerow([run+1, changed_codon, aa, seq_options[new_seq]])

            if new_seq == seq:
                break
            else:
                seq = new_seq

        writer.writerow([seq])
    return seq


# in each round - check all single synonymous codon changes and calculate optimization score - take the best one
def hill_climbing_optimize_aa_bulk_by_zscore(seq: str,
                                             user_input: models.UserInput,
                                             cai_or_tai: str,
                                             max_iter: int,
                                             optimization_method: models.OptimizationMethod):
    """
    hill climbing function for performing codon optimization
    in each iteration - for each codon, change all synonymous codons to a specific one and test the zscore of the new
    sequence after each iteration, select the sequence with the best zscore - if it was not changed since the last
    iteration, break. The maximum number of iterations allowed is "max_iter"
    @seq: str, tested seq
    @inp_dict: input dict after usr_inp code
    @opt_type: 'cai' or 'tai'
    @max_iter: maximal number of iterations to perform
    return: seq, which is the optimized sequence
    """
    optimization_method = models.OptimizationMethod.hill_climbing_average
    seq_options = {}
    score = OptimizationModule.run_module(
        final_seq=seq,
        user_input=user_input,
        cai_or_tai=cai_or_tai,
        optimization_method=optimization_method,
    )
    seq_options[seq] = score

    with open("aa_score_change_by_iteration.csv", "w") as f, open("aa_change_by_iteration.csv", "w") as aa_file:
        aa_score_writer = csv.writer(f)
        aa_score_writer.writerow(["iteration", "score"])

        from modules.shared_functions_and_vars import synonymous_codons
        import collections

        iterations_summary = collections.defaultdict(list)
        i = 0
        for run in range(max_iter):
            def find_best_aa_synonymous_codon(codons_list, seq_to_change):
                aa_seq_options = {}
                for aa_codon in codons_list:
                    option_seq = change_all_codons_of_aa(seq_to_change, aa_codon)
                    aa_seq_options[aa_codon] = OptimizationModule.run_module(
                        final_seq=option_seq,
                        user_input=user_input,
                        cai_or_tai=cai_or_tai,
                        optimization_method=optimization_method,
                    )
                aa_new_seq = max(aa_seq_options, key=aa_seq_options.get)
                return aa_new_seq

            aa_to_selected_codon = {}
            # Find best synonymous_codon per aa
            for aa in synonymous_codons.keys():
                selected_aa_codon = find_best_aa_synonymous_codon(codons_list=synonymous_codons[aa], seq_to_change=seq)
                aa_to_selected_codon[aa] = selected_aa_codon
                iterations_summary[aa].append(selected_aa_codon)

            # create new seq by replacing all synonymous codons
            new_seq = seq
            for aa in aa_to_selected_codon:
                new_seq = change_all_codons_of_aa(new_seq, aa_to_selected_codon[aa])

            # Calculate score after all replacements
            score = OptimizationModule.run_module(
                final_seq=new_seq,
                user_input=user_input,
                cai_or_tai=cai_or_tai,
                optimization_method=optimization_method,
            )
            seq_options[new_seq] = score

            aa_score_writer.writerow([run+1, seq_options[new_seq]])

            if new_seq == seq:
                break

            else:
                seq = new_seq
                i += 1

        aa_score_writer.writerow([seq])
        aa_iterations = csv.writer(aa_file)
        headers = ["aa"]
        headers.extend(list(range(i)))
        aa_iterations.writerow(headers)
        for aa in iterations_summary:
            row = [aa]
            row.extend(iterations_summary[aa])
            aa_iterations.writerow(row)

    return seq
