import glob
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def translate_fasta(input_fasta: str,
                    output_fasta: str,
                    translation_table: int = 2,
                    trim_to_codon: bool = True) -> None:
    """
    Translate DNA sequences in a FASTA file into amino acid sequences.

    :param input_fasta: Path to input DNA FASTA file.
    :param output_fasta: Path to output protein FASTA file.
    :param translation_table: NCBI genetic code table ID.
    :param trim_to_codon: If True, trim to multiple of 3.
    """

    protein_records = []

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = record.seq.upper()

        # Trim sequence to full codons if requested
        if trim_to_codon:
            codon_len = len(seq) - (len(seq) % 3)
            seq = seq[:codon_len]

        protein_seq = seq.translate(table=translation_table, to_stop=False)

        protein_records.append(
            SeqRecord(
                protein_seq,
                id=record.id,
                description=f"translated (table={translation_table})"
            )
        )

    SeqIO.write(protein_records, output_fasta, "fasta")
    print(f"Translated {len(protein_records)} sequences -> {output_fasta}")


if __name__ == "__main__":

    # Directory containing FASTA files
    input_dir = "/Users/lyy/Downloads"

    # Find all .fasta files in the directory
    fasta_files = glob.glob(os.path.join(input_dir, "*.fasta"))

    if not fasta_files:
        print("No FASTA files found in the directory.")
        exit()

    print(f"Found {len(fasta_files)} FASTA files")

    # Use genetic code table 2
    translation_table_id = 2

    for fasta_file in fasta_files:
        base, ext = os.path.splitext(fasta_file)
        output_file = f"{base}_protein.fasta"

        translate_fasta(
            input_fasta=fasta_file,
            output_fasta=output_file,
            translation_table=translation_table_id,
            trim_to_codon=True
        )

