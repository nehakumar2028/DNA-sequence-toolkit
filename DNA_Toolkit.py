
#A command-line DNA sequence analysis toolkit: computes GC content, nucleotide
#composition, reverse complement, and detects/translates open reading frames (ORFs).

def read_file(filepath): #Defining a function to read a file.
     try:
         with open(filepath) as f: #Using the try catch block to handle errors and exceptions
             text=f.read()
             return text
     except FileNotFoundError:
         print("File not found")
         return None

def gc_content(filepath):  # To determine the amount of Guanine and Cytosine bases in the FASTA file
    text = read_file(filepath)  # Reuses the read function to read the file.
    if text is None:
        return None

    g_count = 0
    c_count = 0
    total_length = 0

    for line in text.splitlines(): #Since the file has split lines of DNA sequences
        g_count += line.count('G')
        c_count += line.count('C')
        total_length += len(line)

    gc_percent = ((g_count + c_count) / total_length) * 100
    return gc_percent

def nucleotide_composition(filepath):
    text = read_file(filepath)
    if text is None:
        return None

    composition = {     #Dictionary
        'A': 0,
        'T': 0,
        'G': 0,
        'C': 0
    }

    for line in text.splitlines():
        for nucleotide in line:
            if nucleotide in composition: #Checks 'A','T','G','C' nucleotides in each line
                composition[nucleotide] += 1 #Stores the count of nucleotide

    total = sum(composition.values())

    for nucleotide in composition:
        composition[nucleotide] = (composition[nucleotide] / total) * 100 #Gives percentage of nucleotides

    return composition



def reverse_complement(filepath): #Prints the complementary strand of sequences in FASTA file
    text = read_file(filepath) #Reads the file
    if text is None:
        return None

    complement = {   #Dictionary of complement bases
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G'
    }

    some_string = '' #Declaring a variable to store complementary bases

    for nucleotide in text:
        if nucleotide in complement:
            some_string += complement[nucleotide]

    some_string = some_string[::-1] #Reversing the string

    return some_string


def get_codons(sequence,frame):
    """
        Splits a DNA sequence into codons (3-letter chunks),
        starting from a given reading frame (0, 1, or 2).
        """
    codons = [] # empty list to store each 3-letter codon as we find it
    # range(frame, len(sequence) - 2, 3) walks through valid starting
    # positions for a full codon, jumping 3 characters at a time.
    # "- 2" makes sure we stop before running out of characters for
    # a complete 3-letter codon near the end of the sequence.
    for position in range(frame, len(sequence)-2,3):
        codon = sequence[position:position+3] # grab 3 characters starting here
        codons.append(codon) # add this codon to our growing list
    return(codons)

def find_orfs(codons):
    """This function finds ORFs in a list of codons."""
    orfs = []           # Completed ORFs
    current_orf = []    # ORF currently being built

    for codon in codons:

        # We are not currently inside an ORF
        if not current_orf:
            if codon == "ATG":
                current_orf = [codon]

        # We are currently inside an ORF
        else:
            current_orf.append(codon)

            # Check whether we've reached a stop codon
            if codon in ["TAA", "TAG", "TGA"]:
                orfs.append(current_orf)
                current_orf = []

    return orfs

def find_all_orfs(sequence):
    """This The function find_all_orfs(sequence) is used to find all Open Reading Frames (ORFs) in a DNA sequence."""
    all_orfs = [] #An empty list

    for frame in [0, 1, 2]: #Starting from all frames
        codons = get_codons(sequence, frame)
        orfs = find_orfs(codons)
        all_orfs.extend(orfs) #Adding all the ORFs from the previous function into an empty list

    return all_orfs

def translate(codons):
    """
       Translates a list of codons into a string of amino acids
       using the standard codon table. Stops at the first stop codon.
       """
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    protein=""
    for codon in codons:
        protein += codon_table[codon]
    return protein

def main():
    filepath = "data/insulin_gene.fasta"

    # Read and clean the file once
    text = read_file(filepath)
    if text is None:
        return
    text = text.replace('\n', '')

    # Run each analysis, print results
    print("GC Content:", gc_content(filepath))
    print("Nucleotide Composition:", nucleotide_composition(filepath))
    print("Reverse Complement:", reverse_complement(filepath))

    # Find and translate ORFs
    orfs = find_all_orfs(text)
    print(f"\nORFs found: {len(orfs)}")
    for orf in orfs:
        protein = translate(orf)
        print(f"  Codons: {orf}")
        print(f"  Protein: {protein}")


if __name__ == "__main__":
    main()





















