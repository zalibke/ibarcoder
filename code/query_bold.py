import sys
import requests
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os

def query_bold(sequence):
    """
    Query the BOLD Systems API for a given DNA sequence.
    """
    # BOLD API endpoint for sequence matching
    url = "http://v3.boldsystems.org/index.php/API_Public/sequence?"

    # Set up the payload for the POST request
    params = {
        'sequence': sequence
    }

    try:
        # Send a request to BOLD API
        response = requests.post(url, data=params)
        response.raise_for_status()  # Raise an error for bad responses

        # Parse the JSON response
        data = response.json()

        # Check if there are any matches
        if data and 'bold_records' in data and len(data['bold_records']) > 0:
            # Return the first match (top hit)
            return data['bold_records'][0]
        else:
            return None

    except requests.exceptions.RequestException as e:
        print(f"Error querying BOLD: {e}")
        return None

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
    bold_hits = []

    # Iterate over each sequence
    for record in sequences:
        sequence = str(record.seq)
        print(f"Querying BOLD for sequence ID: {record.id}")

        # Query BOLD for the sequence
        top_hit = query_bold(sequence)

        if top_hit:
            print(f"Top hit for {record.id}:")
            print(f"Taxonomy: {top_hit.get('taxonomy')}")
            print(f"Similarity: {top_hit.get('similarity')}")
            print(f"Specimen ID: {top_hit.get('specimen_id')}")
            print(f"Process ID: {top_hit.get('processid')}")
            print(f"GenBank ID: {top_hit.get('genbank_id')}")
            print(f"Sequence: {top_hit.get('nucleotides')}")

            # Create a SeqRecord for the top hit
            hit_record = SeqRecord(
                Seq(top_hit.get('nucleotides')),
                id=top_hit.get('processid', 'unknown_id'),
                description=f"{top_hit.get('taxonomy', 'unknown_taxonomy')} Similarity: {top_hit.get('similarity', 'N/A')}"
            )

            # Append the top hit to the list
            bold_hits.append(hit_record)
        else:
            print(f"No hit found for {record.id}")

    # Determine output file name based on input file prefix
    file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]
    output_file = f"{file_prefix}_boldhits.fasta"

    # Save the top hits to the new FASTA file
    with open(output_file, "w") as f:
        SeqIO.write(bold_hits, f, "fasta")

    print(f"Top hits saved to {output_file}")

if __name__ == "__main__":
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <fasta_file>")
        sys.exit(1)

    # Get the input FASTA file from command-line arguments
    fasta_file = sys.argv[1]

    # Run the BOLD query for each sequence
    main(fasta_file)
