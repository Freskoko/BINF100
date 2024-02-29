from pathlib import Path


GAP_PENALTY = -2
MATCH = 1
MISMATCH = -1
INITIAL_GAP_PENALTY = 0

SEQ1 = Path("resources/DNA_A_California_2009_pandemicH1N1_segment7.txt").read_text()
SEQ2 = Path("resources/DNA_A_Brisbane_2007_H1N1_M2_CDS.txt").read_text()

def create_table(SEQ1: str, SEQ2: str):
    table = [[0 for _ in range(len(SEQ1)+1)] for _ in range(len(SEQ2)+1)]
    max_value = 0
    max_position = None

    for row_index in range(1, len(table)):
        for col_index in range(1, len(table[0])):
            seq1_base = SEQ1[col_index-1]
            seq2_base = SEQ2[row_index-1]

            if seq1_base == seq2_base:
                DIAG_SCORE = table[row_index-1][col_index-1] + MATCH
            else:
                DIAG_SCORE = table[row_index-1][col_index-1] + MISMATCH

            ABOVE_SCORE = table[row_index-1][col_index] + GAP_PENALTY
            LEFT_SCORE = table[row_index][col_index-1] + GAP_PENALTY

            cell_score = max(0, ABOVE_SCORE, LEFT_SCORE, DIAG_SCORE)

            if cell_score > max_value:
                max_value = cell_score
                max_position = (row_index, col_index)

            table[row_index][col_index] = cell_score

    return table, max_position


def find_backtrack(i: int, j: int, A1: str, A2: str, k: int, table: list[list]):
    #if reached top right STOP
    if (i > 0) or (j>0):

        print(f"recursion : {k} | A1 = {A1}, A2 = {A2} | i,j = {i},{j}")

        #check above
        if ((i>0) and table[i][j] == (table[i-1][j]+GAP_PENALTY)) or i == 0:
            A1mod = A1 + "-"
            A2mod = A2 + SEQ2[i-1]
            find_backtrack(i-1,j,A1mod,A2mod, k+1, table)

        #check left
        if ((j>0) and table[i][j] == (table[i][j-1]+GAP_PENALTY)):
            A1mod = A1 + SEQ1[j-1]
            A2mod = A2 + "-"
            find_backtrack(i,j-1,A1mod,A2mod, k+1, table)


        #check diag
        if (i>0) and (j>0) and ( # see if not on top right
            (table[i][j] == table[i-1][j-1]+MATCH and (SEQ2[i-1] == SEQ1[j-1])) # check for match in seq and scores
            or (table[i][j] == table[i-1][j-1]+MISMATCH and (SEQ2[i-1] != SEQ1[j-1]))): # check for diag if mismatch and score mismatch 
            
            A1mod = A1 + SEQ1[j-1]
            A2mod = A2 + SEQ2[i-1]
            find_backtrack(i-1,j-1,A1mod,A2mod, k+1, table)

    else:
        print(f"found alignments {A1[::-1]:8} | {A2[::-1]:8}")


table, start_position = create_table(SEQ1, SEQ2)
find_backtrack(start_position[0], start_position[1], "", "", 0, table)

# a) AAAGCTCCGATCTCG and TAAAGCAATTTTGGTTTTTTTCCGA
# found alignments AAAGC    | AAAGC 

# b) The segment 7 sequence for the 2009 H1N1 pandemic influenza virus
# (DNA_A_California_2009_pandemicH1N1_segment7.txt) and the coding
# region of the M2 gene from the Brisbane seasonal strain
# (DNA_A_Brisbane_2007_H1N1_M2_CDS.txt), both available from mittUiB

# found alignments GCCTACCAGAAGCGAATGGGAGTGCAGATGCAGCGATTCAAGTG
# ATCCTCTCGTCATTGCAGCAAATATCATTGGGATCTTGCACCTGATATTGTGGATTACTGATCGTCTTTT
# TTTCAAATGTATTTATCGTCGCTTTAAATACGGTTTGAAAAGAGGGCCTTCTACGGAAGGAGTGCCTGAG
# TCCATGAGGGAAGAATATCAACAGGAACAGCAGAGTGCTGTGGATGTTGACGATGGTCATTTTGTCAACA
# TAGAGCTAGAGTAA  
# |
#  GCCTATCAGAAACGAATGGGGGTGCAGATGCAACGATTCAAGTG
# ATCCTCTTGTTGTTGCCGCAAGTATAATTGGGATTGTGCACTTGATATTGTGGATTATTGATCGCCTTTT
# TTCCAAAAGCATTTATCGTATCTTTAAACACGGTTTAAAAAGAGGGCCTTCTACGGAAGGAGTACCAGAG
# TCTATGAGGGAAGAATATCGAGAGGAACAGCAGAATGCTGTGGATGCTGACGATGATCATTTTGTCAGCA
# TAGAGCTAGAGTAA
