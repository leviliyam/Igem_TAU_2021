from Bio.Seq import Seq

synonymous_codons = {
    "C": ["TGT", "TGC"],
    "D": ["GAT", "GAC"],
    "S": ["TCT", "TCG", "TCA", "TCC", "AGC", "AGT"],
    "Q": ["CAA", "CAG"],
    "M": ["ATG"],
    "N": ["AAC", "AAT"],
    "P": ["CCT", "CCG", "CCA", "CCC"],
    "K": ["AAG", "AAA"],
    "*": ["TAG", "TGA", "TAA"],
    "T": ["ACC", "ACA", "ACG", "ACT"],
    "F": ["TTT", "TTC"],
    "A": ["GCA", "GCC", "GCG", "GCT"],
    "G": ["GGT", "GGG", "GGA", "GGC"],
    "I": ["ATC", "ATA", "ATT"],
    "L": ["TTA", "TTG", "CTC", "CTT", "CTG", "CTA"],
    "H": ["CAT", "CAC"],
    "R": ["CGA", "CGC", "CGG", "CGT", "AGG", "AGA"],
    "W": ["TGG"],
    "V": ["GTA", "GTC", "GTG", "GTT"],
    "E": ["GAG", "GAA"],
    "Y": ["TAT", "TAC"],
}
# --------------------------------------------------------------

def optimize_sequence(target_gene, high_expression_organisms, low_expression_organisms, n_initiation_codons=12):
    """

    :param target_gene: Gene object, which is to be optimized
    :param high_expression_organisms: list of Organism objects. The organisms where we want to express the target gene in
    :param low_expression_organisms: list of Organism objects. The organisms where we go not want the expression
    :param n_initiation_codons: number of codons of the sequence which need to be optimized due to the initiation rules
    :return: Optimised gene sequence according to the organisms' features

    The function calculates the difference between the features of each codon. Each feature has its own weight (ratio)
    """

    optimized_sequence = Seq('')

    optimal_codons = find_optimal_codons(high_expression_organisms, low_expression_organisms) # optimal codons->dict(AA:codon)

    target_protein = target_gene.protein_seq

    # optimize the initiation

    optimized_sequence += optimize_initiation(target_gene[:n_initiation_codons*3])

    for aa in target_protein[n_initiation_codons:]:
        optimized_sequence += optimal_codons[aa]

    return optimized_sequence

# --------------------------------------------------------------
def calc_diff(expression_organisms, low_expression_organisms, codons):
    """

    :param expression_organisms: list of expression organisms
    :param low_expression_organisms: list of no_expression organisms
    :param codons: list of codons to calculate the difference for
    :return: dict (codon:score).

    The function finds the difference between the features for each codon in codons list.
    Each feature has ratio, which is actually a weight of the feature in the optimization.



    The function works only for one expression organism and one no-expression organism
    """

    expression_organisms_features = expression_organisms[0].features
    no_expression_organisms_features = low_expression_organisms[0].features

    diff = {}
    for codon in codons:
        diff[codon] = 0
        for i in range(len(expression_organisms_features)):
            if no_expression_organisms_features[i].weights[codon] == 0:
                diff[codon] += 1000
                continue

            diff[codon] += expression_organisms_features[i].ratio * \
                           (expression_organisms_features[i].weights[codon] / no_expression_organisms_features[i].weights[codon])

        # we need to turn the values upside down, because we are looking for minimal value in find_optimal_codons_function

        diff[codon] = 1/diff[codon]

    return diff

# --------------------------------------------------------------
def loss_function(high_expression_organisms, low_expression_organisms, codons, local_maximum):
    """

    :param high_expression_organisms: list of expression organisms
    :param low_expression_organisms: list of no_expression organisms
    :param codons: list of codons to calculate the loss for
    :return: loss - dict (codon:score)

    The function iterates through each feature in each organism, and sums up loss for each codon
    """

    loss = {}
    if local_maximum:
        for high_expression_organism in high_expression_organisms:
            loss = iterate_through_feature([high_expression_organism], codons, loss, high_expression=True)

        for low_expression_organism in low_expression_organisms:
            loss = iterate_through_feature([low_expression_organism], codons, loss, high_expression=False)
    else:
        loss = iterate_through_feature(high_expression_organisms, codons, loss, high_expression=True)
        loss = iterate_through_feature(low_expression_organisms, codons, loss, high_expression=True)

    return loss


# --------------------------------------------------------------
def iterate_through_feature(organisms, codons, loss, high_expression):
    """

    :param organisms: List of organism objects for which the sequence is optimized
    :param codons: list of codons to choose from
    :param loss: loss(dict) taken from a previous iteration
    :param high_expression: Whether current organism is being optimized for
    high expression or low expression
    :return: updated loss dictionary

    The funciton calculates loss for each codon for the organism,
    and adds it to the loss from the previously calculated organisms.
    """


    for feature_name in [feature.index_name for feature in organisms[0].features]:

        max_value = find_max_value_per_feature(organisms, feature_name, codons)

        for organism in organisms:

            feature = [feature for feature in organism.features if feature.index_name == feature_name]
            f = feature[0]

            for codon in codons:
                loss[codon] = 0
                if high_expression:
                    loss[codon] += f.ratio * ((f.weights[codon] / max_value - 1) ** 2)
                else:
                    loss[codon] += f.ratio * ((f.weights[codon] / max_value) ** 2)

    return loss

# --------------------------------------------------------------
def find_max_value_per_feature(organisms, feature_name, codons):

    values = []
    for organism in organisms:
        for feature in organism.features:
            if feature.index_name == feature_name:
                values.extend([feature.weights[codon] for codon in codons])
    max_value = max(values)

    if max_value == 0:
        max_value = 0.000001

    return  max_value
# --------------------------------------------------------------
def find_optimal_codons(high_expression_organisms, low_expression_organisms, evaluation_function=loss_function):
    """

    :param high_expression_organisms: list of Organism objects. The organisms where we want to express the target gene in
    :param low_expression_organisms: list of Organism objects. The organisms where we do not want the expression in
    :param evaluation_function: function which evaluates loss
    :return: Dictionary in the format Amino Acid: Optimal codon.


    """

    optimal_codons = {}

    for aa, codons in synonymous_codons.items():
        loss = evaluation_function(high_expression_organisms, low_expression_organisms, codons, local_maximum=True)
        optimal_codons[aa] = min(loss, key=loss.get)

    return optimal_codons



# --------------------------------------------------------------
def optimize_initiation(seq):
    """

    :param seq: Seq object (only the initiation part)
    :return: optimized initiation sequence

    For now, the function returns the same sequence as in the input, tbd.
    """
    return seq
