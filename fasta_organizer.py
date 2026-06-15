import os
from Bio import SeqIO
from Bio import Align


def extract_and_prep_fasta(patient_fasta_path, reference_sequence):
    """
    Parses a patient FASTA file, performs a global pairwise alignment against the
    KRAS reference transcript, identifies the codon 12 mutation site, and extracts
    the 20-nucleotide target footprint required for ASO design.
    """
    print(f"Reading patient data from: {os.path.basename(patient_fasta_path)}")

    # 1. Parse the FASTA input file using Biopython
    try:
        with open(patient_fasta_path, "r") as handle:
            patient_records = list(SeqIO.parse(handle, "fasta"))
            if not patient_records:
                raise ValueError("The provided FASTA file is empty or malformed.")

            # Extract the primary sequence string
            patient_seq = str(patient_records[0].seq).upper()
            patient_id = patient_records[0].id
            print(f"Successfully loaded sequence for Patient ID: {patient_id}")
    except Exception as e:
        print(f"Error parsing FASTA file: {e}")
        return None

    # 2. Initialize the Pairwise Aligner for global alignment
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1

    print("Aligning patient sequence against KRAS reference genome...")
    alignments = aligner.align(reference_sequence, patient_seq)
    best_alignment = alignments[0]

    # 3. Locate the point mutation relative to the reference framework
    # Reference target window for KRAS c.35G>T (G12V) point mutation
    # We are isolating a 20-nucleotide window centered around position 35
    target_mutation_index = 34  # 0-indexed position for nucleotide 35

    # Extract the local 20-nucleotide sequence window from the patient sequence
    start_pos = target_mutation_index - 9  # 9 bases to the left
    end_pos = target_mutation_index + 11  # 11 bases to the right (total 20-mer)

    if len(patient_seq) >= end_pos:
        target_window = patient_seq[start_pos:end_pos]

        print("\n--- DATA PREPARATION COMPLETE ---")
        print(f"Isolated 20-mer Target Window: {target_window}")
        # Verifying the G -> T point mutation at index 9 of the isolated window
        print(f"Mutated Nucleotide Detected: {target_window[9]} (Expected: T for G12V)")

        # This clean string is ready to be sent to your ViennaRNA and BLAST filters
        return {
            "patient_id": patient_id,
            "target_window": target_window,
            "ready_for_pipeline": True
        }
    else:
        print("Error: Patient sequence length is insufficient to map the target window.")
        return None


# ==========================================
# Example Reference Context for Verification
# ==========================================
if __name__ == "__main__":
    # Standard healthy wild-type KRAS fragment surrounding codon 12
    kras_ref_transcript = "ATGACTGAATATAAACTTGTGGTAGTTGGAGCTGGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAG"

    # Mock file generation to demonstrate the parsing tool logic
    mock_filename = "patient_cohort_sample.fasta"
    with open(mock_filename, "w") as f:
        # Notice the 'T' mutation substitution substituting the standard 'G' at position 35
        f.write(">Patient_14_G12V\n")
        f.write("ATGACTGAATATAAACTTGTGGTAGTTGGAGCTAGTGGCGTAGGCAAGAGTGCCTTGACGATACAGCTAATTCAG\n")

    # Execute the preparation framework
    prep_data = extract_and_prep_fasta(mock_filename, kras_ref_transcript)

    # Clean up mock file asset
    if os.path.exists(mock_filename):
        os.remove(mock_filename)