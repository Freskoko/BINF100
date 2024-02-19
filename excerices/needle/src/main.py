GAP_PENALTY = -1
MATCH = 1
MISMATCH = 0

def create_table(seq1: str, seq2: str):

    #empty table
    table = [   [None for _ in range(len(seq1)+1)]  for _ in range(len(seq2)+1)]

    #always 0
    table[0][0] = 0

    #setup minus for first row
    for i in range(len(table[0])):
        table[0][i] = -1*i

    #setup minus for all rows
    for i in range(len(table)):
        table[i][0] = -1*i

    # ok but dont loop first row and col!!! #TODO
        

    for row_index in range(1,len(table)):
        for col_index in range(1, len(table[0])):

            ABOVE_SCORE = -999
            LEFT_SCORE = -999
            DIAG_SCORE = -999

            seq1_base = seq1[col_index-1]
            seq2_base = seq2[row_index-1]

            print(seq1_base,seq2_base)

            #check above
            if table[row_index-1][col_index] is not None:
                ABOVE_SCORE = table[row_index-1][col_index] + GAP_PENALTY

            #check left
            if table[row_index][col_index-1] is not None:
                LEFT_SCORE = table[row_index][col_index-1] + GAP_PENALTY

            #check diag
            if table[row_index-1][col_index-1] is not None:
                DIAG_SCORE = table[row_index-1][col_index-1]
                if seq1_base == seq2_base:
                    DIAG_SCORE += MATCH
                else:
                    DIAG_SCORE += MISMATCH

            table[row_index][col_index] = max(ABOVE_SCORE,LEFT_SCORE,DIAG_SCORE)


    # print(table)
    for i in table:
        print(i)


def find_backtrack():
    ...

create_table("CACGTAT","CGCA")


