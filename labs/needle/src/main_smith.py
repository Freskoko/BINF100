"""
# BINF100 LAB 1

## 1. (7 points) Test your program with the following sample sequences and see if it can find all matches

### a. Example sequences

    AAAGCTCCGATCTCG and TAAAGCAATTTTGGTTTTTTTCCGA

**Found two optimal alignments:**

Score 5, indexes: 6,5     | CGAAA

Score 5, indexes: 25,10   | AGCCT

### b. Real life sequences

    The segment 7 sequence for the 2009 H1N1 pandemic influenza virus
    (DNA_A_California_2009_pandemicH1N1_segment7.txt) and the coding
    region of the M2 gene from the Brisbane seasonal strain
    (DNA_A_Brisbane_2007_H1N1_M2_CDS.txt), both available from mittUiB

**Found one optimal alignment:**

Score 206, indexes: 294,994

### Other questions

Part 2 and 3 can be seen by running the code.
"""

from pathlib import Path

# Constants
GAP_PENALTY = -2
MATCH = 1
MISMATCH = -1

# Local alignment
INITIAL_GAP_PENALTY = 0

# Global alignments
ALL_ALIGNMENTS = []

# REPLACE THIS WITH YOUR FILE PATH TO THE FILES
YOUR_FILE_SEQ_1 = "src/resources/DNA_A_California_2009_pandemicH1N1_segment7.txt"
YOUR_FILE_SEQ_2 = "src/resources/DNA_A_Brisbane_2007_H1N1_M2_CDS.txt"

# ----------------------------------
# SEQ1 = "AAAGCTCCGATCTCG"
# SEQ2 = "TAAAGCAATTTTGGTTTTTTTCCGA"

# ----------------------------------
SEQ1 = "".join(Path(YOUR_FILE_SEQ_1).read_text().split("\n")[1:])
SEQ2 = "".join(Path(YOUR_FILE_SEQ_2).read_text().split("\n")[1:])


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

    # setup first row
    for i in range(len(table[0])):
        table[0][i] = INITIAL_GAP_PENALTY * i

    # setup first col
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

            table[row_index][col_index] = max(ABOVE_SCORE, LEFT_SCORE, DIAG_SCORE, 0)

    # for i in table:
    #     print(i)

    return table


def alignment_stats(align1: str, align2: str) -> tuple[int, int, int]:
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


def prettyprint_alignments(a1: str, a2: str, start_seq1: int, start_seq2: int) -> None:
    """
    Pretty prints the alignments of two sequences.

    Args:
        a1 (str): The first aligned sequence.
        a2 (str): The second aligned sequence.
        start_seq1 (int): Starting index of a1 in original sequence
        start_seq2 (int): Starting index of a2 in original sequence

    Returns:
        None
    """
    blocks = 20
    build_middle_str = ""

    full_blocks = len(a1) // blocks
    remainder = len(a1) % blocks

    for block_index in range(full_blocks + 1):
        start = block_index * blocks

        if block_index == full_blocks:  # last one not full
            end = start + remainder
            end_len = remainder
        else:  # full block
            end = start + blocks
            end_len = blocks

        for char_a, char_b in zip(a1[start:end], a2[start:end]):
            if char_a == char_b:
                build_middle_str += "|"
            else:
                build_middle_str += " "

        seq1_prefix = f"SEQ1: {start_seq1 + start}"
        seq2_prefix = f"SEQ2: {start_seq2 + start}"

        offset = max(len(seq1_prefix), len(seq2_prefix))

        print(
            seq1_prefix.ljust(offset)
            + f" {a1[start:end]} {start_seq1 + start + end_len}"
        )
        print("".ljust(offset) + " " + build_middle_str)
        print(
            seq2_prefix.ljust(offset)
            + f" {a2[start:end]} {start_seq2 + start + end_len}"
        )
        print("\n")
        build_middle_str = ""


def get_max_score_indices(table: list[list[int]]) -> tuple[list[tuple[int, int]], int]:
    """
    Returns the indices of the maximum score(s) in the given table.
    May return multiple

    Args:
        table (list): A list of tuples containting indices which have the highest score in the table.

    """
    max_score = 0
    max_indices = []

    # find highest
    for i in range(len(table)):
        for j in range(len(table[0])):
            if table[i][j] > max_score:
                max_score = table[i][j]

    # find all indexes that contain max score
    for i in range(len(table)):
        for j in range(len(table[i])):
            if table[i][j] == max_score:
                max_indices.append((i, j))

    return max_indices, max_score


def prettyprint_score(a1: str, a2: str) -> None:
    score, matches, gaps = alignment_stats(a1, a2)
    mismatches = len(a1) - matches

    match_percent = round((matches / len(a1)) * 100, 2)
    mismatch_percent = round((mismatches / len(a1)) * 100, 2)
    gaps_percent = round((gaps / len(a1)) * 100, 2)

    print(
        f"Sequence identity: {matches}/{len(a1)} ({match_percent}%) Mismatches: {mismatches}/{len(a1)} ({mismatch_percent}%) Gaps {gaps}/{len(a1)} ({gaps_percent}%) \n"
    )


def find_best_alignment(
    i: int, j: int, a1: str, a2: str, table: list[list], seq1: str, seq2: str
) -> None:
    """
    Recursively finds the best alignment score between two sequences.
    Updates the global ALL_ALIGNMENTS list with the alignments found.

    Args:
        i (int): The current row index in the table.
        j (int): The current column index in the table.
        A1 (str): The partial alignment of sequence 1.
        A2 (str): The partial alignment of sequence 2.
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        table (list[list]): The scoring table.

    Returns:
        None

    """
    global ALL_ALIGNMENTS

    if (i > 0) and (j > 0):
        # found 0 -> STOP
        if table[i][j] == 0:
            ALL_ALIGNMENTS.append((a1[::-1], a2[::-1]))
            return

        match_score = MATCH if seq1[j - 1] == seq2[i - 1] else MISMATCH

        # check diag
        if i > 0 and j > 0 and table[i][j] == table[i - 1][j - 1] + match_score:
            return find_best_alignment(
                i - 1, j - 1, a1 + seq1[j - 1], a2 + seq2[i - 1], table, seq1, seq2
            )

        # check above
        elif i > 0 and table[i][j] == table[i - 1][j] + GAP_PENALTY:
            return find_best_alignment(
                i - 1, j, a1 + "-", a2 + seq2[i - 1], table, seq1, seq2
            )

        # check left
        elif j > 0 and table[i][j] == table[i][j - 1] + GAP_PENALTY:
            return find_best_alignment(
                i, j - 1, a1 + seq1[j - 1], a2 + "-", table, seq1, seq2
            )

    else:
        # got to top of table -> STOP
        ALL_ALIGNMENTS.append((a1[::-1], a2[::-1]))
        return


if __name__ == "__main__":
    table = create_table(SEQ1, SEQ2)
    max_score_indices, max_score_found = get_max_score_indices(table)
    print("\nOPTIMAL ALIGNMENTS: \n")
    for max_score_index in max_score_indices:
        i, j = max_score_index
        print(f"Indexes with max score ({max_score_found}) {i},{j}) \n")
        find_best_alignment(i, j, "", "", table, SEQ1, SEQ2)

    for alignment, startseq in zip(ALL_ALIGNMENTS, max_score_indices):
        print("------------------- ALIGNMENT FOUND: ------------------- \n")

        a1, a2 = alignment
        seq2_start, seq1_start = startseq

        prettyprint_alignments(a1[::-1], a2[::-1], seq1_start, seq2_start)
        prettyprint_score(a1[::-1], a2[::-1])
    print("------------------------------")
    print(
        f"\nAmount of alignments found with score {max_score_found} = {len(ALL_ALIGNMENTS)} \n"
    )
