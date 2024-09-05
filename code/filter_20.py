# Simple python script to filter consensus sequences for >20 supporting reads
# IMPORTANT: Currently takes the value from supporting_reads = int(header.split("_")[-2]) (AKA second to last character when split by "_")
# Not the most robust method, but works as long as the header is constant. Test rigorously in the future, maybe come up with more rigorous solution?
# Zane Libke
# Last edit - September 5, 2024

import sys
from Bio import SeqIO
import os

def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <fasta_file>")
        sys.exit(1)  # Exit the script with an error code

    # Get the input FASTA file from command-line arguments
    fasta_file = sys.argv[1]

    # Extract the file prefix (base name without extension) for output file naming
    file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]

    # Initialize lists to store sequences
    filtered_seqs = []
    discard_seqs = []

    try:
        # Parse the FASTA file
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract the header and sequence
            header = record.id
            sequence = str(record.seq)

            # Split the header and extract the number of supporting reads
            # Assuming the format is: >ZL052_Boana_nigra_16-SarF_16S-BrR-O_consensus_cl_id_14_total_supporting_reads_201_segment0
            # The number of reads is in the 2nd to last position after splitting by "_"
            # I'm assuming the segment stuff will be consistent if there are multiple?
            supporting_reads = int(header.split("_")[-2])

            # Check the number of supporting reads
            if supporting_reads < 20:
                discard_seqs.append(record)
            else:
                filtered_seqs.append(record)

    except FileNotFoundError:
        print(f"Error: The file '{fasta_file}' was not found.")
        sys.exit(1)

    except ValueError as e:
        print(f"Error parsing supporting reads: {e}")
        sys.exit(1)
    # Define output file names based on the input file prefix
    filtered_output_file = f"../data/{file_prefix}_filtered.fasta"
    discard_output_file = f"../data/{file_prefix}_discard.fasta"

    # Save the filtered sequences to a new FASTA file
    with open(filtered_output_file, "w") as output_file:
        SeqIO.write(filtered_seqs, output_file, "fasta")

    # Save the discarded sequences to another FASTA file
    with open(discard_output_file, "w") as output_file:
        SeqIO.write(discard_seqs, output_file, "fasta")

    # Output the results
    print(f"Total filtered sequences: {len(filtered_seqs)}")
    print(f"Total discarded sequences: {len(discard_seqs)}")
    print(f"Filtered sequences saved to: {filtered_output_file}")
    print(f"Discarded sequences saved to: {discard_output_file}")

if __name__ == "__main__":
    main()
