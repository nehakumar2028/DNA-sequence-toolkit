# DNA Sequence Analysis Toolkit

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/nehakumar2028/DNA-sequence-toolkit/blob/main/v2/DNA_Toolkit_v2_Colab.ipynb)

A Python toolkit for basic DNA sequence analysis: GC content, nucleotide composition, reverse complement, and open reading frame (ORF) detection with protein translation.

It started as a self-directed project to learn Python and core bioinformatics concepts, and was then extended to download sequences from NCBI, validate its input, search both DNA strands, and report ORF coordinates. Version 2 was checked against NCBI's own annotation (see [Validation](#validation)).

## Versions

| | Version 1 | Version 2 |
|---|---|---|
| File | `DNA_Toolkit.py` | `v2/dna_toolkit_v2.py` |
| Input | One local FASTA file | NCBI accession, or a local FASTA file with one or more records |
| ORF search | 3 forward reading frames | Both strands (6 reading frames), with coordinates and a minimum length, ranked longest first |
| Input checks | Basic | Strict (rejects protein, RNA, empty and malformed FASTA) |
| Colab notebook | No | Yes, nothing to upload |

**Version 1** reads a FASTA file and reports GC content, nucleotide composition and the reverse complement, then finds ORFs in the three forward frames and translates them. Edit `filepath` in `main()` to point at a FASTA file (a sample is included: `insulin_gene.fasta`).

## Quick start (Version 2)

### In Google Colab (no installation)

Click the **Open in Colab** badge above, then choose **Runtime > Run all**. Change the accession in the *Run the analysis* cell to analyze a different sequence.

### On your computer

Requires Python 3.8 or newer. There are no third-party dependencies (standard library only).

```bash
git clone https://github.com/nehakumar2028/DNA-sequence-toolkit.git
cd DNA-sequence-toolkit/v2
python dna_toolkit_v2.py
```

This runs the demo on NCBI accession `NM_000618` (human IGF1). The download is saved in a `data/` folder and reused on later runs.

### Use it on your own input

```python
from dna_toolkit_v2 import run_pipeline

run_pipeline("NM_000207")                        # an NCBI accession (downloaded and cached)
run_pipeline("my_sequences.fasta", min_aa=100)   # a local FASTA file, one or more records
run_pipeline("NM_000618", force_refresh=True)    # download again even if cached
```

| Parameter | Default | Meaning |
|---|---|---|
| `filepath_or_accession` | (required) | An NCBI accession or the path of a FASTA file |
| `min_aa` | `50` | Minimum ORF length in amino acids (the stop codon is not counted) |
| `force_refresh` | `False` | Re-download from NCBI even if a cached copy exists |

## Example output

Human IGF1 mRNA (`NM_000618.5`), top five ORFs shown:

```
Downloading NM_000618 from NCBI (Timeout: 15s)...
Successfully saved valid FASTA to data/NM_000618.fasta

======================================================================
Record #1: NM_000618.5 Homo sapiens insulin like growth factor 1 (IGF1), transcri...
======================================================================
Sequence Length: 7277 bp
GC Content (over canonical A/C/G/T): 36.62%
Nucleotide Composition:
  A: 31.83%
  C: 19.53%
  G: 17.09%
  T: 31.55%
Ambiguous Bases Found: None
Reverse Complement (First 50bp): TTAATCAAACATGACTAATTTTAATGTAATTACTAAAGAAAGATATACCA...
ORFs Detected (Minimum 50 aa length): 20
  Rank 1 | Strand [+] | Frame 2 | Coordinates: 183-644 | Size: 153 AA
    Protein Sequence (First 40 AA): MGKISSLPTQLFKCCFCDFLKVKMHTMSSSHLFYLALCLL...
  Rank 2 | Strand [+] | Frame 2 | Coordinates: 252-644 | Size: 130 AA
    Protein Sequence (First 40 AA): MHTMSSSHLFYLALCLLTFTSSATAGPETLCGAELVDALQ...
  Rank 3 | Strand [+] | Frame 2 | Coordinates: 261-644 | Size: 127 AA
    Protein Sequence (First 40 AA): MSSSHLFYLALCLLTFTSSATAGPETLCGAELVDALQFVC...
  Rank 4 | Strand [+] | Frame 1 | Coordinates: 1046-1372 | Size: 108 AA
    Protein Sequence (First 40 AA): MNCFLICTSFYTTRAVCFTVSDNLVSLYPPPPFITFIFAE...
  Rank 5 | Strand [+] | Frame 2 | Coordinates: 1701-1976 | Size: 91 AA
    Protein Sequence (First 40 AA): MPSKRNFRHKSPLMQIGLASPERKLWNGKSRRLGILAVLV...
  ... and 15 more ORFs in rank list.
```

Ranks 2 and 3 are not separate genes. They start at later in-frame ATG codons of the same protein as Rank 1 and share its stop codon (see [How ORFs are defined](#how-orfs-are-defined)).

## Validation

Tested on human IGF1 mRNA (NCBI `NM_000618.5`, 7,277 bp):

- The top-ranked ORF (positions 183-644, 153 aa) matches the coding sequence (CDS) annotated by NCBI for this record exactly.
- BLAST of the translated protein against Swiss-Prot returns IGF1.
- Other ORFs of 50 aa or more fall outside the annotated CDS (for example Ranks 4 and 5, which lie in the long 3' UTR). ORF length alone does not show that a sequence is a gene.

This is a single validated example, not a benchmark.

## What Version 2 does

- **NCBI download:** checks the accession format, saves the result in `data/`, reuses the saved copy, uses a 15 second timeout, and rejects error pages that are not FASTA.
- **Input validation:** gives clear errors for an empty file, a file that does not start with `>`, an empty record, protein-like sequences, RNA (U without T), and invalid characters. Duplicate or empty FASTA headers are kept, not overwritten or dropped.
- **Statistics:** GC content and nucleotide composition are calculated over A, C, G and T only. Ambiguous bases (such as N) are counted and reported separately.
- **Reverse complement:** supports all IUPAC nucleotide codes.
- **ORF detection** on both strands, with translation, 1-based coordinates on the original sequence, and ranking by length.

## How ORFs are defined

- An ORF starts at an `ATG` and ends at the first in-frame stop codon (`TAA`, `TAG` or `TGA`).
- **Every** in-frame `ATG` that reaches a stop codon is reported, so several ORFs in the same frame can share one stop codon (nested ORFs, as in Ranks 1-3 above).
- ORFs that run off the end of the sequence without a stop codon are not reported.
- Frames are numbered 0, 1 and 2 on each strand. Coordinates are 1-based and inclusive, and always refer to the original forward sequence, including for ORFs on the reverse strand.
- Length is in amino acids and excludes the stop codon. Translation uses the standard genetic code.

## Limitations

- ORF scanning is **not gene prediction**. Long ORFs can occur by chance, so results should be checked against annotation or by BLAST.
- It does not handle introns, so it is not suitable for eukaryotic genomic DNA. It works best on mRNA/cDNA and on bacterial sequences.
- Only `ATG` is used as a start codon, and only the standard genetic code is supported. Codons with ambiguous bases are translated as `?`.
- Nested ORFs inflate the ORF count, so the total is not the number of distinct genes.
- Every ORF is kept in memory, so it has not been optimized for whole genomes.
- Results are printed only (the top 5 per record). There is no file export yet.

## Development notes

Version 1 (`DNA_Toolkit.py`) is my own work, written function by function while I was learning Python. Version 2 builds on it and was developed with AI assistance (Claude), which I used to extend the design and review the code. I ran the code, checked the output, and validated it against NCBI annotation as described above.

## Repository structure

```
DNA-sequence-toolkit/
├── README.md
├── DNA_Toolkit.py                  # Version 1
├── insulin_gene.fasta              # sample input for Version 1
└── v2/
    ├── dna_toolkit_v2.py           # Version 2
    └── DNA_Toolkit_v2_Colab.ipynb  # Version 2 as a Colab notebook
```

## Possible next steps

- Compare ORFs with the CDS annotation in a GenBank record automatically
- Add a unit test suite (the logic has been checked against Biopython, but the checks are not in this repository yet)
- Export results to CSV or GFF3
- Support other genetic code tables
- Test on a bacterial genome
