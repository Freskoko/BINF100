from pathlib import Path

# THIS FINDS A LOT OF THEM, WRONG I THINK

GAP_PENALTY = -2
MATCH = 1
MISMATCH = -1
INITIAL_GAP_PENALTY = 0

align_num = 0


def smith_waterman(seq1: str, seq2: str) -> tuple[list[list[int]], int]:
    """
    Performs the Smith-Waterman algorithm to find the optimal local alignment between two sequences.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        tuple: A tuple containing the score matrix and the maximum score.
    """
    m, n = len(seq1), len(seq2)
    score = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0

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
    return score, max_score


def traceback(
    i: int,
    j: int,
    seq1: str,
    seq2: str,
    score: list[list[int]],
    align1: str = "",
    align2: str = "",
) -> None:
    """
    Recursively traces back the alignment path and prints the aligned sequences.

    Args:
        i (int): The current row index in the score matrix.
        j (int): The current column index in the score matrix.
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.
        score (list[list[int]]): The score matrix.
        align1 (str, optional): The aligned sequence 1. Defaults to an empty string.
        align2 (str, optional): The aligned sequence 2. Defaults to an empty string.

    Returns:
        None
    """
    if score[i][j] == 0:
        global align_num
        align_num += 1
        print(f"ALIGNMENT FOUND || {align_num}")
        print(f"SEQ1: {align1[::-1]}")
        print(f"SEQ2: {align2[::-1]}")
        return

    if (
        i > 0
        and j > 0
        and score[i][j]
        == score[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH)
    ):
        traceback(
            i - 1, j - 1, seq1, seq2, score, align1 + seq1[i - 1], align2 + seq2[j - 1]
        )

    if i > 0 and score[i][j] == score[i - 1][j] + GAP_PENALTY:
        traceback(i - 1, j, seq1, seq2, score, align1 + seq1[i - 1], align2 + "-")

    if j > 0 and score[i][j] == score[i][j - 1] + GAP_PENALTY:
        traceback(i, j - 1, seq1, seq2, score, align1 + "-", align2 + seq2[j - 1])


def alignment_score(align1: str, align2: str) -> int:
    """
    Calculates the score of an alignment.

    Args:
        align1 (str): The first aligned sequence.
        align2 (str): The second aligned sequence.

    Returns:
        int: The score of the alignment.
    """
    score = 0
    for a, b in zip(align1, align2):
        if a == "-" or b == "-":
            score += GAP_PENALTY
        elif a == b:
            score += MATCH
        else:
            score += MISMATCH
    return score


def get_all_alignments(seq1: str, seq2: str) -> None:
    """
    Calculate all possible alignments between two sequences using the Smith-Waterman algorithm.

    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        None
    """
    score, max_score = smith_waterman(seq1, seq2)

    for i in range(len(score)):
        for j in range(len(score[i])):
            if score[i][j] == max_score:
                traceback(i, j, seq1, seq2, score)


if __name__ == "__main__":
    seq1 = "AAAGCTCCGATCTCG"
    seq2 = "TAAAGCAATTTTGGTTTTTTTCCGA"

    seq1 = Path("resources/DNA_A_California_2009_pandemicH1N1_segment7.txt").read_text()
    seq2 = Path("resources/DNA_A_Brisbane_2007_H1N1_M2_CDS.txt").read_text()
    get_all_alignments(seq1, seq2)

    # todo we need to score the alignment from the first smith watermna
    # that will be the get score, so then we can look for other possible tracebacks
    # which also achieve the same score.
