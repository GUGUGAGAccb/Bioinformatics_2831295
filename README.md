# README

## Analysis Pipeline

This README describes how to reproduce the sequence processing and translation workflow for samples A–D. All scripts and required reference sequences are included in this package.

---

## 1. Folder Structure

Organise the working directory as follows:

```
project/
│
├── sampleA_part1.FASTQ
├── sampleA_part2.FASTQ
├── sampleA_part3.FASTQ
│
├── sampleB_part1.FASTQ
├── sampleB_part2.FASTQ
├── sampleB_part3.FASTQ
│
├── sampleC_part1.FASTQ
├── sampleC_part2.FASTQ
├── sampleC_part3.FASTQ
│
├── sampleD_part1.FASTQ
├── sampleD_part2.FASTQ
├── sampleD_part3.FASTQ
│
├── sampleA.sh
├── sampleB.sh
├── sampleC.sh
├── sampleD.sh
│
├── Translate.py
│
└── Ref_Seq/       <- extracted from Ref_Seq.zip
```

The `Ref_Seq` directory contains the downloaded reference FASTA files.

---

## 2. Software Requirements

### System

* Unix-like environment (Linux / macOS)

### Tools

* `bash`
* `awk`
* `sed`
* Python 3.8+
* Biopython installed:

```
pip install biopython
```

---

## 3. Generate Sample FASTA Files

Each sample has a dedicated bash script.
Run them one by one from the terminal:

```
bash sampleA.sh
bash sampleB.sh
bash sampleC.sh
bash sampleD.sh
```

Each script:

* extracts sequence lines (`NR % 4 == 2`)
* concatenates three FASTQ fragments
* removes line breaks
* writes output single-line FASTA

Example from `sampleA.sh` :

```
awk 'NR % 4 == 2' "$file" >> "$output"
sequence=$(grep -v '^>' sampleA.fasta | tr -d '\n')
echo ">sampleA" > "$output"
echo "$sequence" >> "$output"
sed 's/ //g' sampleA.fasta
```

After running, you should obtain:

```
sampleA.fasta
sampleB.fasta
sampleC.fasta
sampleD.fasta
```

---

## 4. Combine with Reference Sequences

All reference sequences are located inside `Ref_Seq/` (from `Ref_Seq.zip`).
For each sample, create a working directory and copy:

* 1 sample FASTA
* 4 reference FASTA files

Example:

```
mkdir sampleA_analysis
cp sampleA.fasta sampleA_analysis/
cp Ref_Seq/A_refs/*.fasta sampleA_analysis/
```

Repeat for B, C and D.

---

## 5. Translate DNA to Protein

Run the Python script:

```
python Translate.py
```

The script automatically:

* scans the working directory
* finds all `.fasta` files
* trims sequence length to full codons
* translates using **genetic code table 2**
* writes output as `<name>_protein.fasta`
### AI was using here to fix some bug when running
---

## 6. Multiple Sequence Alignment

For each sample, upload the translated protein FASTA files to **Clustal Omega**:

URL:
[https://www.ebi.ac.uk/Tools/msa/clustalo/](https://www.ebi.ac.uk/Tools/msa/clustalo/)

Settings:

* Input type: protein
* Parameters: default

Download alignment results.

---

## 7. Phylogenetic Tree Reconstruction

Use **Simple Phylogeny** tool:

[https://www.ebi.ac.uk/Tools/phylogeny/simple_phylogeny/](https://www.ebi.ac.uk/Tools/phylogeny/simple_phylogeny/)

Settings:

* Tree format: **Phylip**
* Method: **Neighbour-Joining**
* Exclude gaps: **ON**
* Remaining parameters: default

Export the tree and visualise with **iTOL** if required.

---

## 8. Expected Outputs

For each sample:

* combined FASTA (e.g., `sampleA.fasta`)
* protein FASTA (e.g., `sampleA_protein.fasta`)
* alignment result (Clustal Omega)
* phylogenetic tree (Simple Phylogeny / iTOL)

---
