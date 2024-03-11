from io import TextIOWrapper

from utils.codon_dict import CODES

# group 7


def file_to_reversed_string(f: TextIOWrapper) -> str:
    """
    Open file as list, throw away first line, join to string and reverse string
    """
    # linjer = f.readlines()[1:]
    # linjer_til_str = "".join(linjer)
    # reverse_linjer = linjer_til_str[::-1]

    return "".join(f.readlines()[1:])[::-1]


def translate(sequence: str) -> str:
    aa_seq = ""
    for i in range(
        0, len(sequence) - 3, 3
    ):  # fra 0 til len(seq)-3, hopp over 3 hver gang
        aa_seq += CODES[sequence[i : i + 3]]

    return aa_seq


def compare_sequences(outfile: TextIOWrapper, seq1: str, seq2: str):
    ms_match = 0
    for index, aas in enumerate(zip(seq1, seq2)):
        if aas[0] != aas[1]:
            ms_match += 1

            outstr = (
                f"Mismatch found at index {index}: seq1: {aas[0]}, seq2: {aas[1]} \n"
            )
            print(outstr)
            outfile.write(outstr)

    if ms_match > 0:
        print(f"mismatches (if any) found : {ms_match}")


def main(fname):
    with open("src/inputs/wild_type.txt", "r") as person1, open(
        f"src/inputs/{fname}.txt", "r"
    ) as person2, open(f"src/outputs/outfile_{fname}.txt", "w") as outfile:
        print(f"DATA FOR FILE : {fname}")
        dna_1 = file_to_reversed_string(person1)
        dna1_sequence = dna_1.replace("T", "U")  # transcription
        dna1_sequence = translate(dna1_sequence.upper())

        dna_2 = file_to_reversed_string(person2)
        dna2_sequence = dna_2.replace("T", "U")  # transcription
        dna2_sequence = translate(dna2_sequence).upper()

        compare_sequences(outfile, dna_1, dna_2)
        print(f"END DATA FOR FILE : {fname}")

    return


if __name__ == "__main__":
    main("mutant_1")
    main("mutant_2")
