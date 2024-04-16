# Original code written by Susanna Roeblitz
# Reused by Hermann Holstad Walaunet for BINF100 2023
"""Mandatory assignment 3"""

import numpy as np

# Sequences:
seq1a = "AATTTT"
seq1b = "GCAATTTT"
seq2 = "TAAAGCAATTTTGGTTTTTTTCCGA"


# Initialization
N = len(seq1a)
M = len(seq2)

# initialize (M+1) by (N+1) matrix
H_A = np.zeros((M + 1, N + 1))

# scoring scheme
gap = -2
mismatch = -1
match = 1


# Step 1: Build Alignment matrix
def build_alignment_matrix(seq1, seq2):

    # scoring scheme
    gap = -2
    mismatch = -1
    match = 1

    N = len(seq1)
    M = len(seq2)
    H = np.zeros((M + 1, N + 1))

    for j in range(1, N + 1, 1):
        # initialize 1st row with zeros (free initial gaps in string 2)
        H[0][j] = 0
    for i in range(1, M + 1, 1):
        # initialize 1st column with zeros (free initial gaps in string 1)
        H[i][0] = 0
    for i in range(1, M + 1, 1):
        for j in range(1, N + 1, 1):
            if seq1[j - 1] == seq2[i - 1]:
                score1 = H[i - 1][j - 1] + match
            else:
                score1 = H[i - 1][j - 1] + mismatch
            score2 = H[i][j - 1] + gap
            score3 = H[i - 1][j] + gap
            # take max of candidate scores and zero:
            H[i][j] = max(0, score1, score2, score3)

    return H


# Function to scramble a sequence into a random new sequence
def scramble_sequence(seq: str) -> str:
    scrambled_seq = "".join(np.random.permutation(list(seq)))
    return scrambled_seq


# Function to calculate the Gumbel distribution p-value
def gumbel_distribution(y: float, y_values: list) -> float:
    y_bar = np.mean(y_values)
    s = np.std(y_values, ddof=1)

    # Constant, 90th percentile of the Gumbel distribution
    GUMBEL_SCALE_PARAMETER = 1.282
    # Constant, approximate value of the Euler-Macheroni constant, denoted gamma
    VARIANCE_CONSTANT = 0.577

    # Calculate lambda and mu
    lambda_ = GUMBEL_SCALE_PARAMETER / s
    mu = y_bar - VARIANCE_CONSTANT / lambda_
    p_value = 1 - np.exp(-np.exp(-lambda_ * (y - mu)))

    return p_value


# a)
print("a)")

H_A = build_alignment_matrix(seq1a, seq2)
H_B = build_alignment_matrix(seq1b, seq2)
print("optimal local alignment score for sequences 1a and 2:", np.amax(H_A))
print("optimal local alignment score for sequences 1b and 2:", np.amax(H_B), "\n")

# b)
random_alignments = [scramble_sequence(seq2) for _ in range(1000)]
print("b)\n1000 random alignments have been generated\n")

# c)
optimal_scores_a = [
    np.amax(build_alignment_matrix(seq1a, alignment)) for alignment in random_alignments
]
optimal_scores_b = [
    np.amax(build_alignment_matrix(seq1b, alignment)) for alignment in random_alignments
]
print(
    "c)\ncomputed optimal scores for 1000 randomly scrambled sequences from seq2 with sequence 1a and 1b separately"
)

# d)
print("\nd) ")

original_score_a = np.amax(H_A)
original_score_b = np.amax(H_B)
alpha = 0.01
p_value_a = gumbel_distribution(original_score_a, optimal_scores_a)
p_value_b = gumbel_distribution(original_score_b, optimal_scores_b)
print("Estimated p-values of 1a and 1b using Gumbel distribution")

# e)
print("\ne) ")
print("alpha:", alpha)

print(
    "Q1: p-value for the local alignment scores of sequences 1a and 2:",
    p_value_a.round(5),
)
if p_value_a < alpha:
    print(
        "Results are statistically significant, therefore The null hypothesis is rejected."
    )
    print("Conclusion: The sequences 1a and 2 are homologous.")
else:
    print(
        "Results are not statistically significant, The null hypothesis is not rejected."
    )
    print("Conclusion: The sequences are not homologous.")


print(
    "\nQ2: p-value for the local alignment scores of sequences 1b and 2:",
    p_value_b.round(5),
)
if p_value_b < alpha:
    print(
        "Results are statistically significant, therefore The null hypothesis is rejected."
    )
    print("Conclusion: The sequences 1b and 2 are homologous.")
else:
    print(
        "Results are not statistically significant, The null hypothesis is not rejected."
    )
    print("Conclusion: The sequences are not homologous.")