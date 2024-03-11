from pathlib import Path

# Constants
GAP_PENALTY = -2
MATCH = 1
MISMATCH = -1

# Local alignment
INITIAL_GAP_PENALTY = 0

# ----------------------------------
SEQ1 = "AAAGCTCCGATCTCG" 
SEQ2 = "TAAAGCAATTTTGGTTTTTTTCCGA"

# Found two optimal alignments: 
# Score 5, indexes: 6,5     | CGAAA
# Score 5, indexes: 25,10   | AGCCT

# ----------------------------------
SEQ1 = "".join(Path("src/resources/DNA_A_California_2009_pandemicH1N1_segment7.txt").read_text().split("\n")[1:])
SEQ2 = "".join(Path("src/resources/DNA_A_Brisbane_2007_H1N1_M2_CDS.txt").read_text().split("\n")[1:])

# Found one optimal alignment:
# Score 206, indexes: 294,994

ALL_ALIGNMENTS = []

def create_table(SEQ1: str, SEQ2: str) -> list[list[int]]:
    """
    Creates a dynamic programming table for sequence alignment.

    Args:
        SEQ1 (str): The first sequence.
        SEQ2 (str): The second sequence.

    Returns:
        list: A 2D list representing the dynamic programming table.
    """

    # empty table
    table = [[None for _ in range(len(SEQ1) + 1)] for _ in range(len(SEQ2) + 1)]

    # always 0
    table[0][0] = 0

    # setup minus for first row
    for i in range(len(table[0])):
        table[0][i] = INITIAL_GAP_PENALTY * i

    # setup minus for all rows
    for i in range(len(table)):
        table[i][0] = INITIAL_GAP_PENALTY * i

    # loop through table
    for row_index in range(1, len(table)):
        for col_index in range(1, len(table[0])):

            seq1_base = SEQ1[col_index - 1]
            seq2_base = SEQ2[row_index - 1]

            # check above
            ABOVE_SCORE = table[row_index - 1][col_index] + GAP_PENALTY

            # check left
            LEFT_SCORE = table[row_index][col_index - 1] + GAP_PENALTY

            # check diag
            DIAG_SCORE = table[row_index - 1][col_index - 1]
            # if dna base match
            if seq1_base == seq2_base:
                DIAG_SCORE += MATCH
            else:
                DIAG_SCORE += MISMATCH

            table[row_index][col_index] = max(ABOVE_SCORE, LEFT_SCORE, DIAG_SCORE,0)

    # for i in table:
    #     print(i)

    return table


def alignment_score(align1: str, align2: str) -> int:
    """
    Calculates the score of an alignment.

    Args:
        align1 (str): The first aligned sequence.
        align2 (str): The second aligned sequence.

    Returns:
        int: The score of the alignment.
        int: amount of matches
        int: amount of gaps
    """
    score = 0
    matches = 0
    gaps = 0
    for a, b in zip(align1, align2):
        if a == "-" or b == "-":
            score += GAP_PENALTY
            gaps += 1
        elif a == b:
            score += MATCH
            matches += 1
        else:
            score += MISMATCH
    return score, matches, gaps


def prettyprint_alignments(a1: str, a2: str) -> None:
    """
    Pretty prints the alignments of two sequences.

    Args:
        a1 (str): The first aligned sequence.
        a2 (str): The second aligned sequence.

    Returns:
        None
    """
    blocks = 40
    build_middle_str = ""

    # full block size, remainder size
    full_blocks = len(a1) // blocks
    remainder = len(a1) % blocks

    for block_index in range(full_blocks + 1):

        start = block_index * blocks

        if block_index == full_blocks:  # last one not full
            end = start + remainder
        else:  # full block
            end = start + blocks

        for char_a, char_b in zip(a1[start:end], a2[start:end]):
            # create middle alignment segment
            if char_a == char_b:
                build_middle_str += "|"
            else:
                build_middle_str += " "

        first_line = len(f"SEQ1: {start} ")
        print(f"SEQ1: {start} {a1[start:end]} {end}")
        print(f"{' '*first_line}{build_middle_str}")
        print(f"SEQ2: {start} {a2[start:end]} {end}")
        print("\n")
        build_middle_str = ""


def prettyprint_score(a1: str, a2: str) -> None:
    score,matches,gaps = alignment_score(a1,a2)
    mismatches = len(a1) - matches

    match_percent =     round((matches/len(a1)) * 100, 2)
    mismatch_percent =  round((mismatches/len(a1)) * 100, 2)
    gaps_percent =      round((gaps/len(a1)) * 100, 2)

    print(f"Sequence identity: {matches}/{len(a1)} ({match_percent}%) Mismatches: {mismatches}/{len(a1)} ({mismatch_percent}%) Gaps {gaps}/{len(a1)} ({gaps_percent}%) \n")


def get_max_score_indices(table: list[list[int]]) -> list[tuple[int,int]]:
    """
    Returns the indices of the maximum score(s) in the given table.
    May return multiple

    Args:
        table (list): A 2D list representing the table of scores.

    """
    max_score = 0
    max_indices = []

    #find highest
    for i in range(len(table)):
        for j in range(len(table[0])):
            if table[i][j] > max_score:
                max_score = table[i][j]

    #find all indexes that contain max score
    for i in range(len(table)):
        for j in range(len(table[i])):

            if table[i][j] == max_score:
                max_indices.append( (i, j) )

    return max_indices, max_score


def find_best_alignment(i: int, j: int, a1: str, a2: str, table: list[list], seq1: str, seq2: str) -> int:
    """
    Recursively finds the best alignment score between two sequences.

    Args:
        i (int): The current row index in the table.
        j (int): The current column index in the table.
        A1 (str): The partial alignment of sequence 1.
        A2 (str): The partial alignment of sequence 2.
        table (list[list]): The scoring table.

    Returns:    
        int: The best alignment score.

    """
    global ALL_ALIGNMENTS

    if (i > 0) and (j > 0):

        if table[i][j] == 0:
            ALL_ALIGNMENTS.append( (a1[::-1],a2[::-1]) ) 
            return
            
        match_score = MATCH if seq1[j - 1] == seq2[i - 1] else MISMATCH

        # check diag
        if i > 0 and j > 0 and table[i][j] == table[i - 1][j - 1] + match_score:
            return find_best_alignment(
                i - 1, j - 1, a1 + seq1[j - 1], a2 + seq2[i - 1], table, seq1, seq2
            )

        # check above
        elif i > 0 and table[i][j] == table[i - 1][j] + GAP_PENALTY:
            return find_best_alignment(i - 1, j, a1 + "-", a2 + seq2[i - 1], table, seq1, seq2)

        # check left
        elif j > 0 and table[i][j] == table[i][j - 1] + GAP_PENALTY:
            return find_best_alignment(i, j - 1, a1 + seq1[j - 1], a2 + "-", table, seq1, seq2)

    else:
        ALL_ALIGNMENTS.append( (a1[::-1],a2[::-1]) )
        return

if __name__ == "__main__":

    table = create_table(SEQ1, SEQ2)

    # find incides of table with highest score
    max_score_indices, max_score_found = get_max_score_indices(table)

    print("POSSIBLE OPTIMAL ALIGNMENTS: \n")
    for max_score_index in max_score_indices:

        print("------------------- ALIGNMENT FOUND: -------------------")

        i,j = max_score_index

        print(f"Indexes with max score ({max_score_found}) {i},{j}) \n")

        find_best_alignment(
            i,j, "", "", table, SEQ1, SEQ2
        )

    #print good alignments
    for alignment in ALL_ALIGNMENTS:
        a1,a2 = alignment
        score,_,_ = alignment_score(a1[::-1], a2[::-1])
        prettyprint_alignments(a1[::-1], a1[::-1])
        prettyprint_score(a1[::-1], a2[::-1])

    print("------------------------------")
    print(f"\nAmount of alignments found with score {max_score_found} = {len(ALL_ALIGNMENTS)}")



# todo fix first line in files