from aligner_helpers import *
from compsci260lib import *


def run_global_aligner():
    """
    Given two sequences of either DNA or amino acids, initialize the
    appropriate substitution matrix, run the global aligner and report the
    optimal alignment score.
    """
    dict = get_fasta_dict("students.fasta")
    print(global_aligner(dict["Student_1"], dict["Student_7"], ))
    seq1 = dict["Student_1"]
    seq2 = dict["Student_7"]
    match = 2  # Using the scoring function indicated in Problem 1
    mismatch = -1
    gap_penalty = 1  # Gap penalty specified in Problem 1

    seq_type = validate_sequences(seq1, seq2)
    if seq_type == 1:
        # Both the sequences are DNA sequences so use the scores for match and
        # mismatch
        subst_matrix = create_subst_matrix_dna(match, mismatch)
    elif seq_type == 2:
        # Both the sequences are protein sequences so read in the BLOSUM62
        # substitution matrix
        subst_matrix = create_subst_matrix_aa("BLOSUM62.txt")
    else:
        sys.exit("Input sequences are of different types: not both DNA or both protein")

    # Obtain a dictionary of scores for aligning a pair of characters
    subst_dict = create_subst_matrix_dict(subst_matrix)

    optimal_score = solve_global_aligner(seq1, seq2, subst_dict, gap_penalty)

    print(f"Optimal Score: {optimal_score}")


def solve_global_aligner(seq1, seq2, subst_dict, gap_penalty):
    """The overall procedure for collecting the inputs, running the aligner,
    and displaying the table and return the optimal alignment score

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): penalty for a column containing a gap char (g: use a
            positive value because this value will be subtracted)

    Returns:
        (int) the optimal alignment score
    """

    # Initialize the DP table's data structure
    # as a list of lists of ints
    dp_table = [[0] * (len(seq2)+1) for _ in range(len(seq1)+1)]

    # Compute the score of the optimal global alignment
    max_value = global_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table)

    # Display the dp table
    display_dp_table(seq1, seq2, dp_table)

    return max_value


def global_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table):
    """A dynamic programming algorithm that takes two sequences and returns the
    score of the optimal alignment.

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dict): substitution matrix stored as a dictionary, with
            keys that reference the two characters being aligned, and values
            being the corresponding score.  See the create_subst_matrix_dict()
            function to know how this works.

        gap_penalty (int): linear gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

        dp_table (list of list of ints): dynamic programming table, in the
            structure of dp_table[i][j]

    Returns:
        (int): the optimal alignment score
    """

    # the dp table has len(seq1) + 1 rows and len(seq2) + 1 columns
    I = len(dp_table)      # so I is 1 more than m
    J = len(dp_table[0])   # so J is 1 more than n

    # Initialize the dp table with solutions to base cases using linear gap
    # penalty
    gap = -gap_penalty  # Initialize a variable to subtract the gap penalty
    row_penalty = gap   # Use these variables to initialize the base cases of the DP table
    col_penalty = gap
    for k in range(I):
        for l in range(J):
            if k == l == 0:  # Skip this cell; display_dp_table handles this case.
                continue
            elif k == 0:  # Fill all columns with the appropriate penalty
                dp_table[k][l] = col_penalty
                col_penalty = col_penalty + gap
            elif l == 0:  # Fill all rows with the appropriate penalty
                dp_table[k][l] = row_penalty
                row_penalty = row_penalty + gap
            else:  # Initialize other cells to None for later processing.
                dp_table[k][l] = None
    # Compute the scores for the rest of the matrix,
    # i.e. all the elements in dp_table[i][j] for i > 0 and j > 0.
    #

    for row in range(I):
        for col in range(J):
            if dp_table[row][col] is None:  # Update all "None" cells to their appropriate score using the
                # update function discussed in class
                pair = seq1[row - 1] + seq2[col - 1]  # Find the appropriate pair or "column"
                score = subst_dict[pair]  # Find its score in the substitution dictionary
                T1 = score + dp_table[row - 1][col - 1]  # Calculate type I, II, and III, scores
                T2 = dp_table[row - 1][col] - gap_penalty
                T3 = dp_table[row][col - 1] - gap_penalty
                dp_table[row][col] = max(T1, T2, T3)  # The max of these three scores is the optimal choice
    # The optimal score is found at the lower right corner of the dp table:
    return dp_table[I-1][J-1]



if __name__ == "__main__":
    run_global_aligner()
