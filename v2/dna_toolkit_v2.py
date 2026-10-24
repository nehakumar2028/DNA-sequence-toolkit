import os
import re
import urllib.request
import urllib.parse
from typing import Dict, List, Optional, Tuple

# --- CONSTANTS & CONFIGURATION ---
NCBI_EMAIL = ""  # Optional: NCBI asks API users for an email. Leave empty in public code.
TOOL_NAME = "ColabSequenceAnalyzer"
TIMEOUT = 15  # seconds

# IUPAC nucleotide complement mapping
IUPAC_COMPLEMENT = {
    'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
    'M': 'K', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D',
    'N': 'N', '-': '-'
}

CODON_TABLE = {
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

# --- NCBI DOWNLOAD & SAFETY ---

def validate_accession(accession: str) -> bool:
    """Validates the format of NCBI accession IDs."""
    return bool(re.match(r'^[A-Z]{1,6}_?[A-Z]{0,6}[0-9]{5,12}(\.[0-9]+)?$', accession.strip().upper()))

def download_ncbi_fasta(accession: str, save_dir: str = "data", force_refresh: bool = False) -> Optional[str]:
    """
    Downloads a FASTA file from the NCBI Nuccore database securely with input safety checking,
    caching, timeouts, error detection, and NCBI metadata compliance.
    """
    accession = accession.strip()
    if not validate_accession(accession):
        print(f"Error: Invalid NCBI Accession format: '{accession}'")
        return None

    os.makedirs(save_dir, exist_ok=True)
    filepath = os.path.join(save_dir, f"{accession}.fasta")

    # Implement caching
    if os.path.exists(filepath) and not force_refresh:
        print(f"Using cached FASTA file for {accession}: {filepath}")
        return filepath

    # Prepare safe URL
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "tool": TOOL_NAME
    }
    if NCBI_EMAIL:
        params["email"] = NCBI_EMAIL
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{urllib.parse.urlencode(params)}"
    print(f"Downloading {accession} from NCBI (Timeout: {TIMEOUT}s)...")

    try:
        req = urllib.request.Request(
            url,
            headers={'User-Agent': f'{TOOL_NAME}/1.0'}
        )
        with urllib.request.urlopen(req, timeout=TIMEOUT) as response:
            content_bytes = response.read()
            content = content_bytes.decode('utf-8')

        # Safety checks on content (Check if we actually got an NCBI error page instead of a FASTA file)
        if "<Error>" in content or "ResourceNotFound" in content or not content.strip().startswith(">"):
            raise ValueError("NCBI returned an error response or invalid structure, not a valid FASTA format.")

        with open(filepath, 'w') as f:
            f.write(content)

        print(f"Successfully saved valid FASTA to {filepath}")
        return filepath
    except Exception as e:
        print(f"Error fetching accession '{accession}': {e}")
        return None

# --- FASTA PARSING & STRICT VALIDATION ---

def read_multi_fasta(filepath: str) -> Dict[str, str]:
    """
    Reads and strictly validates a multi-record FASTA file.
    Supports DNA, raising clean errors for RNA, proteins, empty, or headerless inputs.
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File '{filepath}' does not exist.")

    records = {}
    current_header = None
    current_seq = []

    def add_record(header, seq_parts):
        """Store a record, keeping duplicate or empty headers instead of overwriting/dropping them."""
        name = header or "unnamed"
        unique_name, n = name, 2
        while unique_name in records:
            unique_name = f"{name} (duplicate {n})"
            n += 1
        records[unique_name] = "".join(seq_parts).upper()

    with open(filepath, 'r') as f:
        lines = f.readlines()

    if not lines or all(not line.strip() for line in lines):
        raise ValueError(f"Validation Error: File {filepath} is completely empty.")

    first_line = next(line for line in lines if line.strip())
    if not first_line.strip().startswith('>'):
        raise ValueError("Validation Error: Invalid format. The file does not start with a '>' header line.")

    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if current_header is not None:
                add_record(current_header, current_seq)
            current_header = line[1:].strip()
            current_seq = []
        else:
            current_seq.append(line)
    if current_header is not None:
        add_record(current_header, current_seq)

    # Validate sequence contents for each parsed record
    for header, seq in list(records.items()):
        if not seq:
            raise ValueError(f"Validation Error: Record '{header}' is empty.")

        # Detect Protein sequences by checking for high density of amino acid specific letters
        # E.g., 'E', 'F', 'I', 'L', 'P', 'Q' are common in proteins but highly unlikely to dominate DNA unless ambiguous
        aa_specific = len(re.findall(r'[EFILPQX]', seq))
        if len(seq) > 20 and (aa_specific / len(seq)) > 0.15:
            raise ValueError(f"Validation Error: Record '{header}' appears to be a protein sequence, not DNA.")

        # Detect RNA
        if 'U' in seq and 'T' not in seq:
            raise ValueError(f"Validation Error: Record '{header}' contains Uracil (U) and no Thymine (T); looks like RNA.")

        # Validate that the sequence consists of valid IUPAC characters
        invalid_chars = set(seq) - set(IUPAC_COMPLEMENT.keys())
        if invalid_chars:
            raise ValueError(f"Validation Error: Record '{header}' contains invalid nucleotide characters: {invalid_chars}")

    return records

# --- BIOINFORMATICS PIPELINE MODULES ---

def compute_gc_and_composition(sequence: str) -> Tuple[float, Dict[str, float], Dict[str, int]]:
    """
    Calculates GC content strictly over canonical bases (A, C, G, T) only.
    Identifies and aggregates ambiguous / non-canonical bases.
    """
    canonical_bases = ('A', 'C', 'G', 'T')   # tuple, so the output order is always A, C, G, T
    counts = {b: 0 for b in canonical_bases}
    ambiguous_counts = {}

    for base in sequence:
        if base in canonical_bases:
            counts[base] += 1
        else:
            ambiguous_counts[base] = ambiguous_counts.get(base, 0) + 1

    total_canonical = sum(counts.values())
    gc_pct = ((counts['G'] + counts['C']) / total_canonical * 100) if total_canonical > 0 else 0.0

    composition = {k: (v / total_canonical * 100) if total_canonical > 0 else 0.0 for k, v in counts.items()}
    return gc_pct, composition, ambiguous_counts

def reverse_complement_iupac(sequence: str) -> str:
    """Generates the reverse complement of a sequence accommodating all IUPAC codes."""
    return "".join(IUPAC_COMPLEMENT.get(base, 'N') for base in reversed(sequence))

def get_orfs_in_strand(sequence: str, strand_label: str, min_aa: int) -> List[Dict]:
    """
    Finds all Open Reading Frames (ORFs) on a single strand for all 3 frames.
    Maps correct frame representation (0, 1, 2) and calculates precise genomic coordinates.
    """
    orfs = []
    seq_len = len(sequence)

    for frame in [0, 1, 2]:
        # Loop through codons
        for i in range(frame, seq_len - 2, 3):
            codon = sequence[i:i+3]
            if codon == "ATG":
                # Start codon found, trace downstream in-frame until stop codon
                current_orf_codons = []
                for j in range(i, seq_len - 2, 3):
                    curr_codon = sequence[j:j+3]
                    current_orf_codons.append(curr_codon)
                    if curr_codon in ["TAA", "TAG", "TGA"]:
                        break
                else:
                    # Reached end of sequence without reaching a stop codon
                    continue

                # Translate ORF
                protein = ""
                for c in current_orf_codons:
                    aa = CODON_TABLE.get(c, '?')
                    if aa == '*':
                        break
                    protein += aa

                if len(protein) >= min_aa:
                    # Compute relative sequence coordinates
                    start_pos = i + 1
                    end_pos = i + (len(current_orf_codons) * 3)

                    orfs.append({
                        "strand": strand_label,
                        "frame": frame,
                        "start": start_pos,
                        "end": end_pos,
                        "length_aa": len(protein),
                        "sequence_codons": current_orf_codons,
                        "protein_sequence": protein
                    })
    return orfs

def find_and_rank_all_orfs(sequence: str, min_aa: int = 50) -> List[Dict]:
    """
    Finds all ORFs on both the forward and reverse strands across all 6 frames,
    ranks them longest to shortest, and returns them with assigned rank numbers.
    """
    # Forward strand
    forward_orfs = get_orfs_in_strand(sequence, "+", min_aa)

    # Reverse Strand
    rev_seq = reverse_complement_iupac(sequence)
    reverse_orfs_raw = get_orfs_in_strand(rev_seq, "-", min_aa)

    # Map coordinates back to original forward sequence coordinate frame
    seq_len = len(sequence)
    reverse_orfs = []
    for orf in reverse_orfs_raw:
        # Convert coordinates from the reverse strand perspective back to the forward coordinates
        orig_start = seq_len - orf["end"] + 1
        orig_end = seq_len - orf["start"] + 1

        reverse_orfs.append({
            "strand": "-",
            "frame": orf["frame"],
            "start": orig_start,
            "end": orig_end,
            "length_aa": orf["length_aa"],
            "sequence_codons": orf["sequence_codons"],
            "protein_sequence": orf["protein_sequence"]
        })

    all_orfs = forward_orfs + reverse_orfs
    # Sort by descending amino acid length
    all_orfs.sort(key=lambda x: x["length_aa"], reverse=True)

    # Add a unique rank index to each record
    for rank, orf in enumerate(all_orfs, start=1):
        orf["rank"] = rank

    return all_orfs

# --- COMPREHENSIVE PIPELINE WRAPPER ---

def run_pipeline(filepath_or_accession: str, min_aa: int = 50, force_refresh: bool = False):
    """Main orchestration function that parses and analyzes sequence inputs."""
    # 1. Determine whether input is an NCBI accession or a local file
    if not os.path.exists(filepath_or_accession) and validate_accession(filepath_or_accession):
        filepath = download_ncbi_fasta(filepath_or_accession, force_refresh=force_refresh)
        if not filepath:
            print("Aborted: Could not obtain sequence from NCBI.")
            return
    else:
        filepath = filepath_or_accession

    # 2. Parsing and strict validation
    try:
        records = read_multi_fasta(filepath)
    except Exception as e:
        print(f"Pipeline Interrupted: {e}")
        return

    # 3. Process each record in the multi-FASTA file
    for idx, (header, sequence) in enumerate(records.items(), start=1):
        print(f"\n{'='*70}")
        print(f"Record #{idx}: {header[:70]}...")
        print(f"{'='*70}")

        # Base stats
        seq_len = len(sequence)
        gc_pct, comp, ambiguous = compute_gc_and_composition(sequence)

        print(f"Sequence Length: {seq_len} bp")
        print(f"GC Content (over canonical A/C/G/T): {gc_pct:.2f}%")
        print("Nucleotide Composition:")
        for b, p in comp.items():
            print(f"  {b}: {p:.2f}%")

        if ambiguous:
            print("Ambiguous Bases Found:")
            for b, count in ambiguous.items():
                print(f"  {b}: {count} ({(count/seq_len)*100:.2f}% of sequence)")
        else:
            print("Ambiguous Bases Found: None")

        # IUPAC Reverse Complement
        rev_comp = reverse_complement_iupac(sequence)
        print(f"Reverse Complement (First 50bp): {rev_comp[:50]}...")

        # ORF Search and Rank
        orfs = find_and_rank_all_orfs(sequence, min_aa=min_aa)
        print(f"ORFs Detected (Minimum {min_aa} aa length): {len(orfs)}")

        for orf in orfs[:5]:  # Print the top 5 longest ORFs
            print(f"  Rank {orf['rank']} | Strand [{orf['strand']}] | Frame {orf['frame']} | Coordinates: {orf['start']}-{orf['end']} | Size: {orf['length_aa']} AA")
            print(f"    Protein Sequence (First 40 AA): {orf['protein_sequence'][:40]}...")

        if len(orfs) > 5:
            print(f"  ... and {len(orfs) - 5} more ORFs in rank list.")

# --- PIPELINE DEMONSTRATION ---
if __name__ == "__main__":
    # Demo with Human Insulin Accession NM_000618 (Refresh set to False by default)
    run_pipeline("NM_000618", min_aa=50, force_refresh=False)
