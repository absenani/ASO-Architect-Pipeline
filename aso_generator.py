import ssl
import os
import time
import tkinter as tk
from tkinter import filedialog
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML

# --- 1. SSL BYPASS ---
if (not os.environ.get('PYTHONHTTPSVERIFY', '') and getattr(ssl, '_create_unverified_context', None)):
    ssl._create_default_https_context = ssl._create_unverified_context

# --- 2. TRANSCRIPT INDEX ---
TRANSCRIPT_MAP = {
    "NM_001369786.1": "KRAS Transcript Variant 4",
    "NM_001369787.1": "KRAS Transcript Variant 5",
    "XM_054372015.1": "KRAS Transcript Variant X1 (Predicted)",
    "M54968.1": "KRAS Transforming Protein mRNA (Wild-type)"
}


def get_variant_name(accession_id):
    return TRANSCRIPT_MAP.get(accession_id, "KRAS Transcript (Clinical/Predicted)")


# --- 3. ACCESSIBILITY LOGIC ---
def calculate_accessibility(target_seq):
    """
    Predicts if the mRNA target site is 'Open' or 'Structured'.
    High GC in the target means it likely forms a stable hairpin (Lower Accessibility).
    """
    gc_content = (target_seq.count('G') + target_seq.count('C')) / len(target_seq) * 100

    if gc_content < 45:
        return "HIGH (Open Loop)", gc_content
    elif 45 <= gc_content <= 60:
        return "MEDIUM (Partial Fold)", gc_content
    else:
        return "LOW (Deep Stem)", gc_content


# --- 4. PIPELINE CONFIG ---
TARGET_START = 215
TARGET_END = 235
STABILITY_THRESHOLD = -15.0


class ASOArchitect:
    def __init__(self, ref_fasta, patient_fasta):
        self.ref_record = SeqIO.read(ref_fasta, "fasta")
        self.ref_seq = str(self.ref_record.seq)
        self.patient_fasta = patient_fasta
        self.blast_cache = {}

    def run_safety_check(self, aso_seq, patient_id):
        if aso_seq in self.blast_cache:
            return self.blast_cache[aso_seq]
        try:
            print(f"    > Pinging NCBI BLAST API for {patient_id}...")
            result_handle = NCBIWWW.qblast("blastn", "nt", aso_seq, entrez_query="txid9606[ORGN]")
            blast_record = NCBIXML.read(result_handle)
            off, kras = 0, 0
            for alignment in blast_record.alignments:
                is_kras = "KRAS" in alignment.title.upper()
                for hsp in alignment.hsps:
                    if hsp.identities == len(aso_seq):
                        if is_kras:
                            kras += 1
                        else:
                            off += 1
            self.blast_cache[aso_seq] = (off, kras)
            return (off, kras)
        except Exception as e:
            print(f"    ⚠️ API Error: {e}")
            return None

    def run_pipeline(self):
        print("\n--- ASO-Architect: Accessibility-Aware Pipeline ---")
        healthy_target = self.ref_seq[TARGET_START:TARGET_END]

        for patient in SeqIO.parse(self.patient_fasta, "fasta"):
            patient_target = str(patient.seq[TARGET_START:TARGET_END])
            if patient_target != healthy_target:
                variant_name = get_variant_name(patient.id)

                # DESIGN ASO
                aso_seq = str(Seq(patient_target).reverse_complement().transcribe())

                # THERMODYNAMICS
                gc = aso_seq.count('G') + aso_seq.count('C')
                au = aso_seq.count('A') + aso_seq.count('U')
                dg = (gc * -1.2) + (au * -0.7)

                # NEW: ACCESSIBILITY PREDICTION
                access_label, gc_perc = calculate_accessibility(patient_target)

                print(f"\n[!] Mutation Identified: {patient.id} ({variant_name})")
                print(f"    - ASO Design: {aso_seq}")
                print(f"    - Binding ΔG: {dg:.2f} kcal/mol")
                print(f"    - Target Accessibility: {access_label} ({gc_perc:.1f}% GC)")

                if dg <= STABILITY_THRESHOLD:
                    safety = self.run_safety_check(aso_seq, patient.id)
                    if safety:
                        off, variants = safety
                        print(f"    - Specificity: {variants} KRAS variants | {off} non-target hits.")
                        if off == 0:
                            print(f"    ✅ VERIFIED SAFE")
                time.sleep(1)


if __name__ == "__main__":
    root = tk.Tk();
    root.withdraw()
    selected_path = filedialog.askopenfilename(title="Select Patient FASTA")
    root.destroy()

    if selected_path:
        pipeline = ASOArchitect("reference.fasta", selected_path)
        pipeline.run_pipeline()