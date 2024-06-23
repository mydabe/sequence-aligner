from compsci260lib import *
from aligner_helpers import *
from random import randint


def run_local_aligner_plus():
    """Locally align O18381.fasta and P63015.fasta and report the
    optimal score and optimal local alignment information.
    """
    # Import the appropriate sequences to be locally aligned
    dict1 = get_fasta_dict("P63015.fasta")
    dict2 = get_fasta_dict("O18381.fasta")
    list1 = list(dict1.values())
    list2 = list(dict2.values())

    seq1 = list1[0]
    seq2 = list2[0]
    match = 2  # Using the scoring function indicated in Problem 1, ignore if aligning protein sequences
    mismatch = -1
    gap_penalty = 8  # Gap penalty specified in Problem 4

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

    optimal_score = solve_local_aligner_plus(seq1, seq2, subst_dict, gap_penalty)

    print(f"Local Alignment Optimal Score: {optimal_score[0]}")
    # Indices of optimal score are based on the elements in the tuple returned on line 37
    print_alignment(optimal_score[3][0], optimal_score[3][1])
    print(f"Local Alignment Details: {optimal_score}")

def solve_local_aligner_plus(seq1, seq2, subst_dict, gap_penalty):
    """The procedure for collecting the inputs, running the aligner,
    and returning the score and optimal alignment(s).

    Note that for each returned local alignment, starting positions
    also need to be returned. These are the positions of the first
    character in each aligned sequence relative to the original
    sequence.

    Args:
        seq1 (str): first sequence to match
        seq2 (str): second sequence to match
        subst_dict (dictionary string -> int): dictionary
            representation of the substitution matrix
        gap_penalty (int): gap penalty (penalty per gap character);
            this value should be positive because we will subtract it

    A max score may be in multiple locations, so return the optimal
    score, the locations of all the maxima, and any one optimal
    alignment as a tuple.

    Returns a tuple of:
        (the optimal alignment score as an int,

         the locations of the maxima in the dp table as a list of
         tuples. these positions will include the offset of the
         initialized penalty row and column, so that location (i,j)
         refers to the i-prefix of X and the j-prefix of Y, just as in
         lecture,

         tuple for an optimal alignment)

        The alignment will be in the form:

              (tuple of indices of the characters of the first aligned
               sequence used in the alignment),

              (tuple of indices of the characters of the second aligned
               sequence used in the alignment),

              the first aligned sequence as a string,

              the second aligned sequence as a string)

        As an example with the sequences:

            Sequence 1: TAG
            Sequence 2: TAAGATAAG

        A possible return may be:

            (11, # the optimal score

             # the two maximal locations in the dp table
             [(3, 4), (3, 9)],

             # one possible alignment:
             ((1, 3), # the nt positions mapping TA-G from TAG
              (1, 4), # the nt positions mapping TAAG from TAAGATAAG

              "TA-G", "TAAG") # the sequences of the alignment
            )

        Corresponding to two maxima with an optimal
        alignment score of 11.
    """
    # Initialize the DP table's data structure
    # as a list of lists of ints
    dp_table = [[0] * (len(seq2) + 1) for _ in range(len(seq1) + 1)]

    # Compute the score of the optimal local alignment, the amount of optimal cells and their (i, j) locations,
    # and the nucleotide positions where the alignments begin and end
    max_value, optimal_cells, count, new_seqs, beginnings, ends = local_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table)
    return max_value, optimal_cells, count, new_seqs, [ends, beginnings]

def local_aligner(seq1, seq2, subst_dict, gap_penalty, dp_table):
    # the dp table has len(seq1) + 1 rows and len(seq2) + 1 columns
    I = len(dp_table)  # so I is 1 more than m
    J = len(dp_table[0])  # so J is 1 more than n

    # Initialize the dp table with solutions to base cases using linear gap
    # penalty
    gap = -gap_penalty  # Initialize a variable to subtract the gap penalty
    row_penalty = gap  # Use these variables to initialize the base cases of the DP table
    col_penalty = gap
    max_seen = -float("inf")  # Variable to track the maximum cell seen ("running max")
    optimal_cell_count = 0  # Count variable to track how many optimal cells there are

    optimal_cells = []  # Store the locations of the optimal cells as a list of tuples
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
                # If the maximal score is negative, the current cell is a "STOP" cell.
                if m < 0:
                    dp_table[row][col] = 0, "STOP"
                    continue
                elif m > max_seen:  # If the max is a new global maximum, update max_seen and reset
                    # optimal_cell_count and optimal_cells
                    optimal_cell_count = 1
                    optimal_cells = [(row, col)]
                    max_seen = m
                elif m == max_seen:  # If the current max is another occurrence of the global max, add its
                    # location to the optimal cells list and increment optimal_cell_count
                    optimal_cell_count += 1
                    optimal_cells.append((row, col))
                # Generate types
                types = generate_types(m, T1, T2, T3)
                dp_table[row][col] = m, types
    # Recover the top and bottom alignment sequences from the table using a random optimal cell from optimal_cells
    # list.
    aligned_seq_1, aligned_seq_2, x_start, x_end, y_start, y_end = (
        compute_alignment(dp_table, optimal_cells[randint(0, len(optimal_cells) - 1)], seq1, seq2))
    # Print the requested information:
    optimal_cell = optimal_cells[0]
    return (dp_table[optimal_cell[0]][optimal_cell[1]][0], optimal_cells, optimal_cell_count, (aligned_seq_1, aligned_seq_2),
            (x_start, x_end), (y_start, y_end))

def generate_types(m, T1, T2, T3):  # Generate a list of column types that constitute the optimal choice
    types = ["0"]
    if m == T2:
        types.append("T2")
    if m == T1:
        types.append("T1")
    if m == T3:
        types.append("T3")
    return types


def compute_alignment(dp_table, optimal_location, orig_seq1, orig_seq2):
    # Compute the top alignment by iterating through
    #  the dp_table from an optimal_location
    x = optimal_location[0]
    y = optimal_location[1]
    end_x = x
    end_y = y
    start_x = 0
    start_y = 0
    top_seq_1 = ""
    top_seq_2 = ""
    while x > 0:
        if x >= 1 and y == 0:  # Avoid infinite loop
            break
        if dp_table[x][y][0] == 0:
            break
        while y > 0:
            if y >= 1 and x == 0:  # Avoid infinite loop
                break
            # If a STOP cell is reached, record the starting index of the alignment for both sequences
            if dp_table[x][y][0] == 0:
                start_x = x + 1
                start_y = y + 1
                break
            # Otherwise, continue traversing through the table
            types = dp_table[x][y][1]
            # Access the last element of the list; this will always be the topmost
            # available choice given that generate_types adds elements from bottom-most to top-most order
            # Add the appropriate character to the appropriate sequence
            val = types[len(types) - 1]
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
    # Return the alignment and their locations within the original sequences as a tuple
    return top_seq_1, top_seq_2, start_x, end_x, start_y, end_y

if __name__ == "__main__":
    run_local_aligner_plus()
