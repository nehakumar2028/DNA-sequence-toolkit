# DNA Sequence Analysis Toolkit

A command-line Python toolkit that computes GC content, nucleotide composition, reverse complements of DNA strands, and detects and translates open reading frames (ORFs) from a FASTA file.

## Why I built this

The motivation behind this project was to learn, hands-on, how programming can be used to analyze biological sequences. This is a self-directed learning project built for academic and educational purposes — every function was written and debugged by me as a way to understand both Python fundamentals and basic bioinformatics concepts.

## What it does

- Reads a DNA sequence from a FASTA file
- Calculates **GC content** (% of the sequence that is Guanine/Cytosine)
- Calculates **nucleotide composition** (% of A, T, G, C in the sequence)
- Computes the **reverse complement** of the DNA strand
- Detects **open reading frames (ORFs)** across all 3 forward reading frames
- **Translates** each ORF into its corresponding protein sequence using the standard codon table

## How to run it

1. Place a FASTA file inside the `data/` folder (a sample file is already included: `data/insulin_gene.fasta`)
2. Update the `filepath` variable in `main()` if using a different file
3. Run the script:

```bash
python DNA_Toolkit.py
```

## Example output

Run on the real human insulin gene (NCBI accession `NM_000207`):

```
GC Content: 55.83
Nucleotide Composition: {'A': 19.74, 'T': 16.52, 'G': 30.26, 'C': 33.48}
Reverse Complement: GCTGGTTCAAGGGCTTTATTCCATCTCTCTCGGTGCAGGAGG... (truncated)

ORFs found: 2
  Protein: MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN*
  Protein: ME*
```

The first ORF correctly corresponds to the real preproinsulin protein sequence — confirming the toolkit correctly identified and translated the actual coding region of the gene, not just a random reading frame.

## What I learned

Building this project taught me how Python functions work in practice — not just syntax, but how to structure a program as a set of small, reusable functions that build on each other. I ran into and debugged several real issues along the way: confusing a function's parameter with the value passed into it, misplacing a `return` statement inside a loop (which caused it to exit early), and understanding exactly how Python's indentation rules define what's "inside" a loop versus "after" it. Working through open reading frame detection specifically taught me about reading frames, codon structure, and how a six-letter genetic alphabet maps to a 20-amino-acid protein alphabet through the codon table.

## Possible extensions

This standalone version is well suited for learning and for analyzing small sequences. Future extensions could include using libraries like **Biopython** to handle full-genome-scale FASTA files, sequence alignment, and querying biological databases directly — moving from a self-contained learning tool toward a more scalable, production-style bioinformatics pipeline.

## Project structure

```
DNA_Toolkit.py       # Main script — all functions + main()
data/
  insulin_gene.fasta  # Sample real gene sequence (human insulin, NM_000207)
```