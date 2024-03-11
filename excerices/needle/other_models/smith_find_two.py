from pathlib import Path

# THIS ONLY FINDS TWO

GAP_PENALTY = -2
MATCH = 1
MISMATCH = -1
INITIAL_GAP_PENALTY = 0


def smith_waterman(seq1, seq2):
    m, n = len(seq1), len(seq2)
    score = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_pos = None

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = (-1, -1)
            if seq1[i - 1] == seq2[j - 1]:
                match = (i - 1, j - 1)
            score[i][j] = max(
                0,
                score[i - 1][j] + GAP_PENALTY,
                score[i][j - 1] + GAP_PENALTY,
                score[match[0]][match[1]]
                + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH),
            )
            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    align1, align2 = "", ""
    i, j = max_pos
    while i > 0 and j > 0:
        current_score = score[i][j]
        diagonal_score = score[i - 1][j - 1]
        up_score = score[i][j - 1]
        left_score = score[i - 1][j]

        if current_score == 0:
            break

        if current_score == diagonal_score + (
            MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
        ):
            align1 += seq1[i - 1]
            align2 += seq2[j - 1]
            i -= 1
            j -= 1
        elif current_score == up_score + GAP_PENALTY:
            align1 += "-"
            align2 += seq2[j - 1]
            j -= 1
        elif current_score == left_score + GAP_PENALTY:
            align1 += seq1[i - 1]
            align2 += "-"
            i -= 1

    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2


# seq1 = 'AAAGCTCCGATCTCG'
# seq2 = 'TAAAGCAATTTTGGTTTTTTTCCGA'

seq1 = Path("resources/DNA_A_California_2009_pandemicH1N1_segment7.txt").read_text()
seq2 = Path("resources/DNA_A_Brisbane_2007_H1N1_M2_CDS.txt").read_text()

align1, align2 = smith_waterman(seq1, seq2)
print(align1)
print("\n")
print(align2)
