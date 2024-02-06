from io import TextIOWrapper
from codon_dict import CODES    


def file_to_reversed_string(f: TextIOWrapper) -> str:
    """
    Open file as list, throw away first line, join to string and reverse string
    """
    return("".join(f.readlines()[1:])[::-1]) 

def translate(sequence: str) -> str:
    aa_seq = ""
    for i in range(0,len(sequence)-3,3):
        aa_seq += CODES[ sequence[i:i+3] ]

    return aa_seq

def compare_sequences(outfile: TextIOWrapper, seq1: str, seq2: str):
    ms_match = 0
    for index,aas in enumerate(zip(seq1,seq2)):
        if aas[0] != aas[1]:
            ms_match+=1
            outfile.write(f"Mismatch found at index {index}: seq1: {aas[0]}, seq2: {aas[1]} \n")

    if ms_match>0:
        print(f"mismatches (if any) found : {ms_match}")
    

def main():
    with open("mutant_1.txt", "r") as person1, open(
        "mutant_2.txt", "r"
    ) as person2, open("outfile.txt", "w") as outfile:
        
        dna_1 = file_to_reversed_string(person1) 
        dna1_sequence = dna_1.replace("T","U") #transcription
        dna1_sequence = translate(dna1_sequence)

        dna_2 = file_to_reversed_string(person2) 
        dna2_sequence = dna_2.replace("T","U") #transcription
        dna2_sequence = translate(dna2_sequence)

        compare_sequences(outfile,dna_1,dna_2)
    
    return

if __name__ == "__main__":
    main()