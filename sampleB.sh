#!/bin/bash

# There are 2 header line in sampleB part1, so delete the first line
sed '1d' sampleB_part1.fastq > sampleB_part1_1.fastq

# Output file name
output="sampleB.fasta"

# Remove old output file if it exists
rm -f "$output"

# Loop through the three FASTQ files
for file in sampleB_part1_1.FASTQ sampleB_part2.FASTQ sampleB_part3.FASTQ
do
    # Extract only read lines:
    # In FASTQ format, the sequence is always on every 2nd line starting from line 2.
    # Use awk to print every 2nd line.
    awk 'NR % 4 == 2' "$file" >> "$output"
done

# Remove newlines inside the merged sequence (make one continuous line)
# Store as a temporary file then overwrite output
sequence=$(grep -v '^>' "$output" | tr -d '\n')
echo ">sampleB" > "$output"
echo "$sequence" >> "$output"

# Remove all spaces in the FASTA
sed 's/ //g' sampleB.fasta

echo "Done. Output written to sampleB.fasta"
