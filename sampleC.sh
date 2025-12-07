#!/bin/bash

# Output file name
output="sampleC.fasta"

# Remove old output file if it exists
rm -f "$output"

# Loop through the three FASTQ files
for file in sampleC_part1.FASTQ sampleC_part2.FASTQ sampleC_part3.FASTQ
do
    # Extract only read lines:
    # In FASTQ format, the sequence is always on every 2nd line starting from line 2.
    # Use awk to print every 2nd line.
    awk 'NR % 4 == 2' "$file" >> "$output"
done

# Remove newlines inside the merged sequence (make one continuous line)
# Store as a temporary file then overwrite output
sequence=$(grep -v '^>' "$output" | tr -d '\n')
echo ">sampleC" > "$output"
echo "$sequence" >> "$output"

# Remove all spaces in the FASTA
sed 's/ //g' sampleA.fasta

echo "Done. Output written to sampleC.fasta"
