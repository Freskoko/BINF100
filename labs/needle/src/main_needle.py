GAP_PENALTY = -1
MATCH = 1
MISMATCH = 0

# local alignment
INITAL_GAP_PENALTY = 0

SEQ1 = "CACGTAT"
SEQ2 = "CGCA"


def create_table(SEQ1: str, SEQ2: str):
    # empty table
    table = [[None for _ in range(len(SEQ1) + 1)] for _ in range(len(SEQ2) + 1)]

    # always 0
    table[0][0] = 0

    # setup minus for first row
    for i in range(len(table[0])):
        table[0][i] = INITAL_GAP_PENALTY * i

    # setup minus for all rows
    for i in range(len(table)):
        table[i][0] = GAP_PENALTY * i

    # loop through table
    for row_index in range(1, len(table)):
        for col_index in range(1, len(table[0])):
            seq1_base = SEQ1[col_index - 1]
            seq2_base = SEQ2[row_index - 1]

            # check above
            if table[row_index - 1][col_index] is not None:
                ABOVE_SCORE = table[row_index - 1][col_index] + GAP_PENALTY

            # check left
            if table[row_index][col_index - 1] is not None:
                LEFT_SCORE = table[row_index][col_index - 1] + GAP_PENALTY

            # check diag
            if table[row_index - 1][col_index - 1] is not None:
                DIAG_SCORE = table[row_index - 1][col_index - 1]
                # if dna base match
                if seq1_base == seq2_base:
                    DIAG_SCORE += MATCH
                else:
                    DIAG_SCORE += MISMATCH

            table[row_index][col_index] = max(ABOVE_SCORE, LEFT_SCORE, DIAG_SCORE)

    # print(table)
    for i in table:
        print(i)

    return table


def find_backtrack(i: int, j: int, A1: str, A2: str, k: int, table: list[list]):
    # start low right corner, find highest score, keep going.

    # if reached top right STOP
    if (i > 0) and (j > 0):
        print(f"recursion : {k} | A1 = {A1}, A2 = {A2} | i,j = {i},{j}")

        # check above
        if ((i > 0) and table[i][j] == (table[i - 1][j] + GAP_PENALTY)) or i == 0:
            A1mod = A1 + "-"
            A2mod = A2 + SEQ2[i - 1]
            find_backtrack(i - 1, j, A1mod, A2mod, k + 1, table)

        # check left
        if (j > 0) and table[i][j] == (table[i][j - 1] + GAP_PENALTY):
            A1mod = A1 + SEQ1[j - 1]
            A2mod = A2 + "-"
            find_backtrack(i, j - 1, A1mod, A2mod, k + 1, table)

        # check diag
        if (
            (i > 0)
            and (j > 0)
            and (  # see if not on top right
                (
                    table[i][j] == table[i - 1][j - 1] + MATCH
                    and (SEQ2[i - 1] == SEQ1[j - 1])
                )  # check for match in seq and scores
                or (
                    table[i][j] == table[i - 1][j - 1] + MISMATCH
                    and (SEQ2[i - 1] != SEQ1[j - 1])
                )
            )
        ):  # check for diag if mismatch and score mismatch
            A1mod = A1 + SEQ1[j - 1]
            A2mod = A2 + SEQ2[i - 1]
            find_backtrack(i - 1, j - 1, A1mod, A2mod, k + 1, table)

    else:
        print(f"found alignments {A1[::-1]:8} | {A2[::-1]:8}")
        return (A1[::-1],)


table = create_table(SEQ1, SEQ2)

A1 = ""
A2 = ""

# start at highest vwlu in last row!!
# TODO
right_corner_col = len(table) - 1
right_corner_row = len(table[0]) - 1

largest_last_row_value = max(table[-1])
largest_last_row_value_col_index = table[-1].index(largest_last_row_value)


find_backtrack(
    largest_last_row_value, largest_last_row_value_col_index, A1, A2, 0, table
)
