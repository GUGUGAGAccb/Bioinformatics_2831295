# Step 1: Import libraries
import glob
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# Step 2: Define a function that reads a FASTA file and translates the DNA
def translate_fasta_file(input_file, output_file, table_id=2):
    """
    Read a FASTA file, translate each DNA sequence into protein,
    then save the result into a new FASTA file.
    """

    protein_list = []  # this will store all translated records

    # Step 2.1: Go through every record in the FASTA file
    for record in SeqIO.parse(input_file, "fasta"):

        # make the sequence uppercase (A, T, C, G)
        dna_seq = record.seq.upper()

        # trim the sequence so its length is divisible by 3
        if len(dna_seq) % 3 != 0:
            trim_size = len(dna_seq) - (len(dna_seq) % 3)
            dna_seq = dna_seq[:trim_size]

        # Step 2.2: Translate using Biopython
        protein_seq = dna_seq.translate(table=table_id, to_stop=False)

        # Step 2.3: Create a new sequence record for the translated protein
        new_record = SeqRecord(
            protein_seq,
            id=record.id,
            description="translated using table {}".format(table_id)
        )
        protein_list.append(new_record)

    # Step 2.4: Write all translated sequences to a new FASTA file
    SeqIO.write(protein_list, output_file, "fasta")

    print("Translated {} sequences into {}".format(len(protein_list), output_file))


# Step 3: Main program
if __name__ == "__main__":

    # Step 3.1: Tell the script where to find the input FASTA files
    folder = "/Users/lyy/Desktop/Bioinformatics/Ref_Seq"

    # Step 3.2: Look for files that end with .fasta
    fasta_files = glob.glob(os.path.join(folder, "*.fasta"))

    # Step 3.3: Check if we found anything
    if len(fasta_files) == 0:
        print("No FASTA files found.")
    else:
        print("Found {} FASTA files".format(len(fasta_files)))

        # Step 3.4: Translate every file we found
        for f in fasta_files:
            name, ext = os.path.splitext(f)
            out_name = name + "_protein.fasta"

            translate_fasta_file(
                input_file=f,
                output_file=out_name,
                table_id=2  # use genetic code table 2
            )
