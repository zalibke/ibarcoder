from Bio import SeqIO
import sys

input_file = sys.argv[1]
output_file = "../data/temp_nodup.fasta"

# Dictionary to store unique headers
unique_records = {}

# Read the input FASTA file and filter out duplicates based on headers
for record in SeqIO.parse(input_file, "fasta"):
    if record.id not in unique_records:
        unique_records[record.id] = record

# Write the filtered sequences to a new FASTA file
with open(output_file, "w") as output_handle:
    SeqIO.write(unique_records.values(), output_handle, "fasta")

print(f"Filtered out duplicated headers, saved to {output_file}")
