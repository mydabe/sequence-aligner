from aligner_helpers import *
from helperlibrary import *
from global_aligner import *


def run_global_aligner_plus():
    """Generate the optimal global alignments between:
    For each alignment, report the optimal alignment score,
    the top-most alignment, and the bottom-most alignment.
    """
    # Retrieve appropriate sequences to be aligned (.fasta files).
    dict1 = get_fasta_dict("")
    dict2 = get_fasta_dict("")
    list1 = list(dict1.values())
    list2 = list(dict2.values())

    seq1 = list1[0]
    seq2 = list2[0]
    match = 2  # Using an arbitrary scoring function
    mismatch = -1
    gap_penalty = 8  # Using an arbitrary gap penalty.

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
    #  Call solve_global_aligner_plus to get the details of the top and bottom alignments,
    #  along with their score
    optimal_details = solve_global_aligner_plus(seq1, seq2, subst_dict, gap_penalty)
    #  Print the optimal score to the console
    print(f"Optimal Score: {optimal_details[0]}")
    #  If both the top and bottom alignments are the same, print only one of them
    if len(optimal_details) == 2:
        print("Optimal Alignment: ")
        print_alignment(optimal_details[1][0], optimal_details[1][1])
    #  Else print both of them
    else:
        print("Top Alignment: ")
        print_alignment(optimal_details[1][0], optimal_details[1][1])
        print("Bottom Alignment: ")
        print_alignment(optimal_details[2][0], optimal_details[2][1])



def solve_global_aligner_plus(seq1, seq2, subst_dict, gap_penalty):
    """The overall procedure for collecting the inputs, running the aligner,
    filling in the DP table, and returning the final value and alignments.

    Args:
        seq1 (str): first sequence to be aligned
        seq2 (str): second sequence to be aligned
        subst_dict (dictionary string -> int): dictionary representation of the
            substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character); this
            value should be positive because we will subtract it

    Returns a tuple of:
        (the optimal alignment score as an int,
         the top-most alignment achieving this score as a tuple of strings,
         the bottom-most alignment achieving this score as a tuple of strings)
    """

    # Initialize the DP table's data structure
    # as a list of lists of ints
    dp_table = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]

    # Compute the score and optimal alignments
    solved = global_aligner_plus(seq1, seq2, subst_dict, gap_penalty, dp_table)

    # Return
    return solved
def global_aligner_plus(seq1, seq2, subst_dict, gap_penalty, dp_table):
    # the dp table has len(seq1) + 1 rows and len(seq2) + 1 columns
    I = len(dp_table)      # so I is 1 more than m
    J = len(dp_table[0])   # so J is 1 more than n

    # Initialize the dp table with solutions to base cases using linear gap
    # penalty
    gap = -gap_penalty  # Initialize a variable to subtract the gap penalty
    row_penalty = gap   # Use these variables to initialize the base cases of the DP table
    col_penalty = gap

    choices = []
    for k in range(I):
        for l in range(J):
            if k == l == 0:  # Add a None object so that each value of the dp_table is a tuple.
                dp_table[k][l] = 0, None
            elif k == 0:  # Fill all columns with the appropriate penalty
                dp_table[k][l] = col_penalty, ["T3"]
                col_penalty = col_penalty + gap
            elif l == 0:  # Fill all rows with the appropriate penalty
                dp_table[k][l] = row_penalty, ["T2"]
                row_penalty = row_penalty + gap
            else:  # Initialize other cells to None for later processing.
                dp_table[k][l] = None
    # Compute the scores for the rest of the matrix,
    # i.e. all the elements in dp_table[i][j] for i > 0 and j > 0.
    for row in range(I):
        for col in range(J):
            if dp_table[row][col] is None:  # Update all "None" cells to their appropriate score using the
                # update function discussed in class
                pair = seq1[row - 1] + seq2[col - 1]  # Find the appropriate pair or "column"
                score = subst_dict[pair]  # Find its score in the substitution dictionary
                T1 = score + dp_table[row - 1][col - 1][0]  # Calculate type I, II, and III, scores
                T2 = dp_table[row - 1][col][0] - gap_penalty
                T3 = dp_table[row][col - 1][0] - gap_penalty
                m = max(T1, T2, T3)
                # The max of these three scores is the optimal choice. Also record the column type of
                # the cell to compute top and bottom alignments.
                types = generate_types(m, T1, T2, T3)
                # Each cell should contain the optimal score of the i, j column and the possible types that
                # yield that score.
                dp_table[row][col] = m, types
    # Recover the top and bottom alignment sequences from the table
    top_aligned_1, top_aligned_2 = compute_alignment(dp_table, seq1, seq2, I, J)
    bottom_aligned_1, bottom_aligned_2 = compute_alignment(dp_table, seq1, seq2, I, J, bottom=True)
    # Return the requested information:
    if top_aligned_1 == bottom_aligned_1 and top_aligned_2 == bottom_aligned_2:
        return dp_table[I - 1][J - 1][0], (top_aligned_1, top_aligned_2)
    return dp_table[I-1][J-1][0], (top_aligned_1, top_aligned_2), (bottom_aligned_1, bottom_aligned_2)


def generate_types(m, T1, T2, T3): # Generate a list of column types that constitute the optimal choice
    types = []
    if m == T2:
        types.append("T2")
    if m == T1:
        types.append("T1")
    if m == T3:
        types.append("T3")
    return types


def compute_alignment(dp_table, orig_seq1, orig_seq2, I, J, bottom=False):
    # Compute the top alignment by iterating through
    #  the dp_table from the cell in the lower right corner.
    x = I - 1
    y = J - 1
    top_seq_1 = ""
    top_seq_2 = ""
    while x > 0:
        if x >= 1 and y == 0:  # Avoid infinite loop
            break
        while y > 0:
            if y >= 1 and x == 0:  # Avoid infinite loop
                break
            # Access the types list, the second value of the tuple stored in each cell of dp_table
            types = dp_table[x][y][1]
            # Access the last element of the list; this will always be the topmost
            # available choice given that generate_types adds elements from bottom-most to top-most order
            val = types[len(types) - 1]
            # If the bottom argument is true, compute the bottom alignment by selecting the bottom-most option
            # available
            if bottom:
                val = types[0]
            # Add the appropriate character to the appropriate sequence
            if val == "T3":
                top_seq_2 = orig_seq2[y - 1] + top_seq_2
                top_seq_1 = "-" + top_seq_1
                y -= 1
            elif val == "T1":
                top_seq_2 = orig_seq2[y - 1] + top_seq_2
                top_seq_1 = orig_seq1[x - 1] + top_seq_1
                y -= 1
                x -= 1
            else:
                top_seq_1 = orig_seq1[x - 1] + top_seq_1
                top_seq_2 = "-" + top_seq_2
                x -= 1
    # Add any remaining characters of either sequence
    while y != 0:
        top_seq_2 = orig_seq2[y - 1] + top_seq_2
        top_seq_1 = "-" + top_seq_1
        y -= 1

    while x != 0:
        top_seq_1 = orig_seq1[x - 1] + top_seq_1
        top_seq_2 = "-" + top_seq_2
        x -= 1
    # Return the two sequences as a tuple
    return top_seq_1, top_seq_2
if __name__ == "__main__":
   run_global_aligner_plus()
