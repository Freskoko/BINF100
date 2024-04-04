"""
# BINF100 LAB 3: Evaluation of Alignment Scores

a) (2 points) Compute the optimal local alignment score for the sequence pair by using
dynamic programming.

### SEQ1A and SEQ 2: (score 6)

Indexes with max score (6) 12,6) 

------------------- ALIGNMENT FOUND: ------------------- 

SEQ1: 6  TTTTAA 12
         ||||||
SEQ2: 12 TTTTAA 18


Sequence identity: 6/6 (100.0%) Mismatches: 0/6 (0.0%) Gaps 0/6 (0.0%) 

------------------------------

Amount of alignments found with score 6 = 1 


### SEQ1B and SEQ 2: (score 8)


Indexes with max score (8) 12,8) 

------------------- ALIGNMENT FOUND: ------------------- 

SEQ1: 8  TTTTAACG 16
         ||||||||
SEQ2: 12 TTTTAACG 20


Sequence identity: 8/8 (100.0%) Mismatches: 0/8 (0.0%) Gaps 0/8 (0.0%) 

------------------------------

Amount of alignments found with score 8 = 1 


b) (2 points) Generate 1000 random sequences by shuffling the characters from seq2.

See `generate_n_random_sequences()`


c) (2 points) Compute the optimal scores for aligning each of the 1000 random sequences with
seq1a or seq1b, respectively. Hint: You only need to compute the optimal scores, not the
alignments

See code below `if __name__ == "__main__"`


d) (2 points) Estimate the p-value of the original score from a) based on the statistics of scores
from c). Hint: Use the right statistics, i.e. the one for local alignment scores

See code below `if __name__ == "__main__"`


e) (2 points) Based on a significance level alpha=0.01, give answers to questions 1. and 2. from above


1. Are seq1a and seq2 homologous?

SEQ1A:

P VALUE OF 0.06, ALPHA = 0.01
P VALUE OF 0.06 OBSERVED, SEQUENCES ARE NOT HOMOLOGOUS

P values observed are of course affected by the random sequences, 
but seemed to stay around 0.055

No, these sequences are not homologous

2. Are seq1b and seq2 homologous?

SEQ1B:

P VALUE OF 0.004, ALPHA = 0.01
P VALUE OF 0.004 OBSERVED, SEQUENCES ARE HOMOLOGOUS

P values observed are of course affected by the random sequences, 
but seemed to stay around 0.004

Yes, these sequences are homologous

"""

import random

# I just put in a seed so that its reproducable
random.seed("BINF100 IS FUN :)")

GAP_PENALTY = -2
MATCH = 1
MISMATCH = -1

INITIAL_GAP_PENALTY = 0

ALL_ALIGNMENTS = []

def create_table(SEQ1: str, SEQ2: str) -> list[list[int]]:
    """
    Creates a dynamic programming table for sequence alignment.

    Args:
        SEQ1 (str): The first sequence.
        SEQ2 (str): The second sequence.

    Returns:
        list: A 2D list of ints representing the dynamic programming table.
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



def full_waterman_process(seq1_inp,seq2_inp) -> int:
    """
    Carries out the full smith-waterman alogrithm between
    two sequences. 

    Returns the score found, and the seq2 
    """
    table = create_table(seq1_inp, seq2_inp)
    max_score_indices, max_score_found = get_max_score_indices(table)
    print("\nOPTIMAL ALIGNMENTS: \n")
    for max_score_index in max_score_indices:
        i, j = max_score_index
        print(f"Indexes with max score ({max_score_found}) {i},{j}) \n")
        find_best_alignment(i, j, "", "", table, seq1_inp, seq2_inp)

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

    return max_score_found


# b) (2 points) Generate 1000 random sequences by shuffling the characters from seq2.
def generate_n_random_sequences(n: int,string_inp: str) -> list[str]:
    """
    Generates n randomly shuffled variations of string_inp
    
    Returns a list of length n, containing the randomly shuffled strings

    Used this source:
    https://stackoverflow.com/questions/2668312/shuffle-string-in-python
    """
    rand_sequences = []
    for _ in range(n):
        rand_seq = 2

        rand_seq= list(string_inp)
        random.shuffle(rand_seq)
        rand_seq = ''.join(rand_seq)

        rand_sequences.append(rand_seq)

    return rand_sequences

if __name__ == "__main__":
    ALPHA = 0.01

    SEQ1A="AATTTT"
    SEQ1B="GCAATTTT"
    SEQ2="TAAAGCAATTTTGGTTTTTTTCCGA"

    SEQ1 = SEQ1B
    SEQ2 = SEQ2
    
    #a) (2 points) Compute the optimal local alignment score for the sequence pair by using
    # dynamic programming. 
    max_score = full_waterman_process(SEQ1,SEQ2)

    #b) (2 points) Generate 1000 random sequences by shuffling the characters from seq2.
    random_seq2_sequences = generate_n_random_sequences(1000,SEQ2)

    #c) (2 points) Compute the optimal scores for aligning each of the 1000 random sequences with
    # seq1a or seq1b, respectively. 
    score_dict = {}
    for rand_seq in random_seq2_sequences:
        # i know i dont need to do the full alignment, but this works 
        score = full_waterman_process(SEQ1,rand_seq)
        score_dict[rand_seq] = score

    #d) (2 points) Estimate the p-value of the original score from a) based on the statistics of scores
    #from c).

    # count how many reach max score
    # HMM UNSURE ABOUT THIS PART !!! 
    matching_max_count = len([k for k in score_dict if score_dict[k] >= max_score])
    print(f"AMOUNT MATCHING OR ABOVE MAX COUNT : {matching_max_count}")
    random_matching_max = matching_max_count / 1000

    #e) (2 points) Based on a significance level alpha=0.01, give answers to questions 1. and 2. from above
    print(f"P VALUE OF {random_matching_max}, ALPHA = {ALPHA}")
    print(f"P VALUE OF {random_matching_max} OBSERVED, SEQUENCES ARE {'NOT ' if random_matching_max>ALPHA else ''}HOMOLOGOUS")
    
    # SEQ1A:

    # P VALUE OF 0.06, ALPHA = 0.01
    # P VALUE OF 0.06 OBSERVED, SEQUENCES ARE NOT HOMOLOGOUS

    # SEQ1B:

    # P VALUE OF 0.004, ALPHA = 0.01
    # P VALUE OF 0.004 OBSERVED, SEQUENCES ARE HOMOLOGOUS





    


    


