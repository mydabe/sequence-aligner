"""
COMPSCI 260 Python library

A basic library module that contains some useful features.
"""
import re
import sys


def get_fasta_dict(filename):
    """Given a fasta input file, return a dictionary with the name of each
       sequence as a key, and the actual sequence as the value.  Thus, if you
       called dict = get_fasta_dict(filename), you could see the sequence named
       "sequencename" (given that it exists with that name in the fasta file)
       by refering to dict["sequencename"].  To see all of the sequence names, 
       you could use dict.keys()."""

    filename.rstrip()
    fasta_dict = {}

    f = open(filename, "r")

    try:
        curr_seq = ''
        curr_seq_key = ''

        for line in f.readlines():
            # get rid of leading and trailing line
            line = line.strip()

            if line == '' or line.isspace() or line[0] == '#':
                continue
            elif line[0] == '>':
                if not curr_seq_key == '' and not curr_seq == '':
                    fasta_dict[curr_seq_key] = curr_seq
                curr_seq = ''
                curr_seq_key = line[1:]
            else:
                curr_seq += line


            # update the latest line
        fasta_dict[curr_seq_key] = curr_seq
        f.close()

    except:
        print("Error reading fasta file: make sure the file format is correct")
        raise

    return fasta_dict


def complement_base(base, seq_type='DNA'):
    """Returns the complement of a base.
    Args:
        base (str): one nucleotide to complement.
        seq_type (str): the type of nucleotide sequence 'RNA' or 'DNA'.

    Returns:
        The complement of the input base as a string.
    """

    # Check the seq_type
    assert seq_type in ['DNA', 'RNA'], "seq_type must be DNA or RNA"

    # Lookup table for the input bases
    complements_DNA = {"A": "T", "C": "G", "G": "C", "T": "A",
                       "a": "t", "c": "g", "g": "c", "t": "a"}
    complements_RNA = {"A": "U", "C": "G", "G": "C", "U": "A",
                       "a": "u", "c": "g", "g": "c", "u": "a"}
    try:
        if seq_type == 'DNA':
            assert base in 'ACGTacgt', "Make sure sequence contains exact "\
                "DNA nucleotide characters only (i.e. A, C, G, or T)"
            return complements_DNA[base]
        elif seq_type == 'RNA':
            assert base in 'ACGUacgu', "Make sure sequence contains exact "\
                "RNA nucleotide characters only (i.e. A, C, G, or U)"
            return complements_RNA[base]

    except KeyError:
        print ("Sequence type", seq_type, "does not agree with provided base",
               base)


def reverse_complement(sequence, seq_type='DNA'):
    """Reverse complement an input DNA or RNA nucleotide sequence.

    Args:
        sequence (str): nucleotide sequence to reverse complement.
        seq_type (str): the type of nucleotide sequence 'RNA' or 'DNA'.

    Returns:
        A string of the reverse complement of the input sequence.
    """

    # Complement the sequence
    complement = ''
    for base in sequence:
        complement += complement_base(base, seq_type)

    # Reverse the complement sequence
    return complement[::-1]


codon_to_aa_dict = {}
codon_to_aa_dict[r"GC."] = "A"                          # Alanine
codon_to_aa_dict[r"[TU]G[TUC]"] = "C"                   # Cysteine
codon_to_aa_dict[r"GA[TUC]"] = "D"                      # Aspartic Acid
codon_to_aa_dict[r"GA[AG]"] = "E"                       # Glutamic Acid
codon_to_aa_dict[r"[TU][TU][TUC]"] = "F"                # Phenylalanine
codon_to_aa_dict[r"GG."] = "G"                          # Glycine
codon_to_aa_dict[r"CA[TUC]"] = "H"                      # Histidine
codon_to_aa_dict[r"A[TU][TUCA]"] = "I"                  # Isoleucine
codon_to_aa_dict[r"AA[AG]"] = "K"                       # Lysine
codon_to_aa_dict[r"[TU][TU][AG]|C[TU]."] = "L"          # Leucine
codon_to_aa_dict[r"A[TU]G"] = "M"                       # Methionine
codon_to_aa_dict[r"AA[TUC]"] = "N"                      # Asparagine
codon_to_aa_dict[r"CC."] = "P"                          # Proline
codon_to_aa_dict[r"CA[AG]"] = "Q"                       # Glutamine
codon_to_aa_dict[r"CG.|AG[AG]"] = "R"                   # Arginine
codon_to_aa_dict[r"[TU]C.|AG[TUC]"] = "S"               # Serine
codon_to_aa_dict[r"AC."] = "T"                          # Threonine
codon_to_aa_dict[r"G[TU]."] = "V"                       # Valine
codon_to_aa_dict[r"[TU]GG"] = "W"                       # Tryptophan
codon_to_aa_dict[r"[TU]A[TUC]"] = "Y"                   # Tyrosine
codon_to_aa_dict[r"[TU]A[AG]|[TU]GA"] = "*"             # Stop


def codon_2_amino_acid(codon):
    """Uses regular expressions to translate one codon to the respective amino acid
       (works whether the codon is DNA or RNA, and whether it is provided in lower
       or uppercase (since we uppercase before translating))."""

    for codon_regex in codon_to_aa_dict:
        if re.search(codon_regex, codon.upper()) is not None:
            return codon_to_aa_dict[codon_regex]

    print("Unresolvable codon: " + codon)
    sys.exit('Error!')  # if failed, then output "Error!" to console


def translate(nucleic_acid_seq):
    """A subroutine to translate an input nucleic acid sequence into a peptide."""

    protein = ""

    for i in range(0, len(nucleic_acid_seq)-2, 3):
        protein += codon_2_amino_acid(nucleic_acid_seq[i:i+3])

    return protein


def max_over_indices(array, indices):
    """Return the maximum value over the indices (list) given. If indices
       are out of range, a None object will be returned."""

    if any(t < 0 for t in indices):
        return None
    if max(indices) >= len(array):
        return None

    temp = []
    for i in indices:
        temp.append(array[i])

    return max(temp)


def print_abbrv(string):
    """Abbreviate long strings for printing to the console: certain versions/
       configurations of Eclipse have a bug where long strings will be printed
       "invisibly" for some mysterious reason.  That is, a variable may
       actually contain a long string, but when you ask Eclipse to print it to
       the screen, it doesn't seem to appear.  Hence, this function can help
       you ensure your variable really contains what you think it contains."""

    l = len(string)
    if l < 100:
        print(string)
    else:
        ellipsis = "...%d more characters..." % (l-70)
        print(string[:50] + ellipsis + string[-20:])
