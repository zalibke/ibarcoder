import sys
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML

def blast_sequence(sequence):
    """
    Query the NCBI BLAST service for a given DNA sequence.
    """
    try:
        # Perform BLAST search against the nt database
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

        # Parse the BLAST results
        blast_record = NCBIXML.read(result_handle)

        # Check if there are any matches
        if len(blast_record.alignments) > 0:
            # Return the first match (top hit)
            first_alignment = blast_record.alignments[0]
            first_hsp = first_alignment.hsps[0]
            return first_alignment, first_hsp
        else:
            return None, None

    except Exception as e:
        print(f"Error querying BLAST: {e}")
        return None, None

def main(fasta_file):
    # Read sequences from the FASTA file
    try:
        sequences = list(SeqIO.parse(fasta_file, "fasta"))
    except FileNotFoundError:
        print(f"Error: The file '{fasta_file}' was not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)

    # List to store the top hits
    blast_hits = []

    # Iterate over each sequence
    for record in sequences:
        sequence = str(record.seq)
        print(f"Querying BLAST for sequence ID: {record.id}")

        # Query BLAST for the sequence
        top_alignment, top_hsp = blast_sequence(sequence)

        if top_alignment and top_hsp:
            print(f"Top hit for {record.id}:")
            print(f"Title: {top_alignment.title}")
            print(f"Length: {top_alignment.length}")
            print(f"E-value: {top_hsp.expect}")
            print(f"Score: {top_hsp.score}")
            print(f"Identity: {top_hsp.identities}/{top_hsp.align_length}")
            print(f"Query Strand: {top_hsp.strand[0]}, Subject Strand: {top_hsp.strand[1]}")
            print(f"Sequence: {top_hsp.sbjct}")

            # Create a SeqRecord for the top hit
            hit_record = SeqRecord(
                Seq(top_hsp.sbjct),
                id=top_alignment.hit_id,
                description=f"{top_alignment.title} E-value: {top_hsp.expect}"
            )

            # Append the top hit to the list
            blast_hits.append(hit_record)
        else:
            print(f"No hit found for {record.id}")

    # Determine output file name based on input file prefix
    file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    output_file = f"../data/{file_prefix}_blasthits.fasta"

    # Save the top hits to the new FASTA file
    with open(output_file, "w") as f:
        SeqIO.write(blast_hits, f, "fasta")

    print(f"Top hits saved to {output_file}")

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <fasta_file>")
        sys.exit(1)

    # Get the input FASTA file from command-line arguments
    fasta_file = sys.argv[1]

    # Run the BLAST query for each sequence
    main(fasta_file)
