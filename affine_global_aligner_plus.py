from helperlibrary import *
from aligner_helpers import *


def run_ag_aligner_plus():
    """Align atpa_Hs.fasta and atpaEc.fasta and report the optimal
    alignment score, the top-most alignment, and the bottom-most
    alignment.
    """
    # Import the appropriate sequences to be aligned
    dict = get_fasta_dict("students.fasta")
    seq1 = dict["Student_3"]
    seq2 = dict["Student_7"]
    match = 2  # Using the scoring function indicated in Problem 1. If protein sequence, ignore this
    mismatch = -1
    gap_penalty = 1  # Gap penalty specified in Problem 3
    affine_penalty = 11

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
    # Retrieve the details of the optimal top and bottom affine alignments along with their score
    optimal_details = solve_ag_aligner_plus(seq1, seq2, subst_dict, gap_penalty, affine_penalty)
    # Print the score to the console
    print(f"Optimal Score w/ Affine Gap Penalty: {optimal_details[0]}")
    # If the bottom and top alignment are the same (one unique alignment), print only one
    if len(optimal_details) == 2:
        print("Unique Affine Alignment: ")
        print_alignment(optimal_details[1][0], optimal_details[1][1])
    # Print both alignments
    else:
        print("Top Affine Alignment: ")
        print_alignment(optimal_details[1][0], optimal_details[1][1])
        print("Bottom Affine Alignment: ")
        print_alignment(optimal_details[2][0], optimal_details[2][1])

def solve_ag_aligner_plus(seq1, seq2, subst_dict, gap_penalty, affine_penalty):
    """The procedure for collecting the inputs, running the aligner,
    and returning the score and optimal alignments.

    Args:
        seq1 (str): first sequence to match
        seq2 (str): second sequence to match
        subst_dict (dictionary string -> int): dictionary
            representation of the substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character);
            this value should be positive because we will subtract it
        affine_penalty (int): affine penalty; as a positive integer

    Returns a tuple of:
        (the optimal alignment score as an int,
         the top-most alignment achieving this score as a tuple of
         strings, the bottom-most alignment achieving this score as a
         tuple of strings)

        Example output:
            (6, ("AT-AGG", "ATCCGG"), ("ATA-GG", "ATCCGG"))
    """

    # Initialize the DP table's data structure
    # as a list of lists of ints
    dp_table = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]

    # Compute the details of the alignments using the global_ag_aligner_plus function
    alignment_details = global_ag_aligner_plus(seq1, seq2, subst_dict, gap_penalty, affine_penalty, dp_table)
    #  Return a tuple as described above
    return alignment_details

def global_ag_aligner_plus(seq1, seq2, subst_dict, gap_penalty, affine_penalty, t1_table):
    I = len(t1_table)  # so I is 1 more than m
    J = len(t1_table[0])  # so J is 1 more than n
    gap = -gap_penalty  # Initialize a variable to subtract the gap penalty
    affine = -affine_penalty  # Initialize a variable to subtract the affine penalty
    # Create two additional tables for types two and three
    t2_table = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    t3_table = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]
    # Initialize base case cells
    for k in range(I):
        for l in range(J):
            if k == l == 0:  # Initialize the top left corner of each table as zero
                t1_table[k][l] = 0
                t2_table[k][l] = 0
                t3_table[k][l] = 0
            elif k == 0:  # Fill all columns with the appropriate penalty
                # These cells should be infinity since T3 columns cannot appear in tables 1 or 2
                t1_table[k][l] = -float("inf")
                t2_table[k][l] = -float("inf")
                if l == 1:  # If the current cell signifies the beginning of a T3 block...
                    t3_table[k][l] = affine + gap
                else:
                    t3_table[k][l] = t3_table[k][l - 1] + gap
            elif l == 0:  # Fill all rows with the appropriate penalty
                # These cells should be infinity since T2 columns cannot appear in tables 1 or 3
                t1_table[k][l] = -float("inf")
                t3_table[k][l] = -float("inf")
                if k == 1:  # If the current cell signifies the beginning of a T2 block...
                    t3_table[k][l] = affine + gap
                else:
                    t2_table[k][l] = t2_table[k - 1][l] + gap
            else:  # Initialize other cells to None for later processing.
                t1_table[k][l] = None
                t2_table[k][l] = None
                t3_table[k][l] = None

    # Compute the scores for the rest of the matrix,
    # i.e. all the elements in dp_table[i][j] for i > 0 and j > 0.
    for row in range(I):
        for col in range(J):
            if t1_table[row][col] is None:  # Update all "None" cells to their appropriate score using the
                # update function discussed in class
                pair = seq1[row - 1] + seq2[col - 1]  # Find the appropriate pair or "column"
                score = subst_dict[pair]  # Find its score in the substitution dictionary
                # Access the necessary dependencies to fill table 1
                t1_val = t1_table[row - 1][col - 1]
                t2_val = t2_table[row - 1][col - 1]
                t3_val = t3_table[row - 1][col - 1]
                # Update the score of table 1 (match/mismatch table)
                m = max(t2_val, t3_val, t1_val)
                t1_table[row][col] = score + m
                # Access the necessary dependencies to fill table 2
                t1_val = t1_table[row - 1][col] + affine
                t2_val = t2_table[row - 1][col]
                t3_val = t3_table[row - 1][col] + affine
                # Update the score of table 2 (type 2 matrix)
                t2_table[row][col] = gap + max(t1_val, t2_val, t3_val)
                # Access the necessary dependencies to fill table 3
                t1_val = t1_table[row][col - 1] + affine
                t2_val = t2_table[row][col - 1] + affine
                t3_val = t3_table[row][col - 1]
                # Update the score of table 3 (type 3 matrix)
                t3_table[row][col] = gap + max(t1_val, t2_val, t3_val)
    # Compute the top alignment (see comments on compute_affine_top_alignment)
    aligned_sequences_top = compute_affine_alignment(t1_table, t2_table, t3_table, seq1, seq2, I, J)
    # Compute the bottom alignment
    aligned_sequences_bottom = compute_affine_alignment(t1_table, t2_table, t3_table, seq1, seq2, I, J, bottom=True)
    # Check if the alignments are the same
    top_seq_1 = aligned_sequences_top[0]
    top_seq_2 = aligned_sequences_top[1]
    bottom_seq_1 = aligned_sequences_bottom[0]
    bottom_seq_2 = aligned_sequences_bottom[1]
    if top_seq_1 == bottom_seq_1 and top_seq_2 == bottom_seq_2: # Only add unique alignments to the return tuple
        return max(t1_table[I - 1][J - 1], t2_table[I - 1][J - 1], t3_table[I - 1][J - 1]), (top_seq_1, top_seq_2)
    # Return the optimal score and the sequences of both the bottom and top alignments.
    return (max(t1_table[I - 1][J - 1], t2_table[I - 1][J - 1], t3_table[I - 1][J - 1]), (top_seq_1, top_seq_2),
            (bottom_seq_1, bottom_seq_2))


def generate_types(m, T1, T2, T3): # Generate a list of column types that constitute the optimal choice
    types = []
    if m == T2:
        types.append("T2")
    if m == T1:
        types.append("T1")
    if m == T3:
        types.append("T3")
    return types


def compute_affine_alignment(t1_table, t2_table, t3_table, orig_seq1, orig_seq2, I, J, bottom=False):
    x = I - 1
    y = J - 1
    top_seq_1 = ""
    top_seq_2 = ""
    # Access the values of the last cell in each table
    t1_val = t1_table[I - 1][J - 1]
    t2_val = t2_table[I - 1][J - 1]
    t3_val = t3_table[I - 1][J - 1]
    # The maximum of these values should begin the traceback
    m = max(t1_val, t2_val, t3_val)
    # Generate the types that give the optimal score m; each type now represents which table to search in next.
    traceback = generate_types(m, t1_val, t2_val, t3_val)

    while x > 0:
        if x >= 1 and y == 0:  # Avoid infinite loop
            break
        while y > 0:
            if y >= 1 and x == 0:  # Avoid infinite loop
                break
            # Access the last element of the list; this will always be the topmost
            # available choice given that generate_types adds elements from bottom-most to top-most order
            val = traceback[len(traceback) - 1]
            # If the bottom argument is true, compute the bottom alignment by selecting the bottom-most option
            # available
            if bottom:
                val = traceback[0]
            # Add the appropriate character to the appropriate sequence
            if val == "T3":
                top_seq_2 = orig_seq2[y - 1] + top_seq_2
                top_seq_1 = "-" + top_seq_1
                # Given that the previous column was T3, find the optimal type of the x, y - 1 column
                t1_val = t1_table[x][y - 1]
                t2_val = t2_table[x][y - 1]
                t3_val = t3_table[x][y - 1]
                m = max(t1_val, t2_val, t3_val)
                # Generate a new traceback list
                traceback = generate_types(m, t1_val, t2_val, t3_val)
                y -= 1
            elif val == "T1":
                top_seq_2 = orig_seq2[y - 1] + top_seq_2
                top_seq_1 = orig_seq1[x - 1] + top_seq_1
                # Given that the previous column was T1, find the optimal type of the x - 1, y - 1 column
                t1_val = t1_table[x - 1][y - 1]
                t2_val = t2_table[x - 1][y - 1]
                t3_val = t3_table[x - 1][y - 1]
                m = max(t1_val, t2_val, t3_val)
                # Generate a new traceback list
                traceback = generate_types(m, t1_val, t2_val, t3_val)
                y -= 1
                x -= 1
            else:
                top_seq_1 = orig_seq1[x - 1] + top_seq_1
                top_seq_2 = "-" + top_seq_2
                # Given that the previous column was T2, find the optimal type of the x - 1, y column
                t1_val = t1_table[x - 1][y]
                t2_val = t2_table[x - 1][y]
                t3_val = t3_table[x - 1][y]
                m = max(t1_val, t2_val, t3_val)
                # Generate a new traceback list
                traceback = generate_types(m, t1_val, t2_val, t3_val)
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
    run_ag_aligner_plus()
