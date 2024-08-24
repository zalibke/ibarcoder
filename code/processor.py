import os
import glob
import subprocess
import pandas as pd
from Bio import SeqIO
import retrieve_inat
import sys

#samples = sys.argv[1]
#primers_db = sys.argv[2]
#inat = FALSE

## functions ##

def list_matching_files(directory, pattern):
    return glob.glob(os.path.join(directory, pattern))
    
# Function to read a FASTA file and return the header and sequence as a tuple
def read_fasta(fasta_file):
    headers = []
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq))
    return headers, sequences

samples = pd.read_csv("../data/samples.csv")
primers_db = pd.read_csv("../data/primers_db.txt")

#retrieve missing inat data from samples with an inat_link field
#if sys.argv[3] = TRUE - make this an optional parameter
#if inat:
#	samples = retrieve_inat.retrieve_inat_all(samples)

#init lists and dicts to store  headers and sequences
all_headers = []
all_sequences = []
consensus_dict = {}


for idx, row in samples.iterrows():
    index_number = row['index_number']
    barcode_folder = "barcode" + str(index_number).zfill(2)
    sample_code = row['sample_code']
    primerF = row['primerF']
    primerR = row['primerR']
    genus_prior = row['genus_prior']
    species_prior = row['species_prior']
    
    # make temp primer file from ../data/primers_db.txt with primerF/primerR combo 
    # code for this here
    # for testing, premade ../data/primers_temp.txt
    input = "../data/" + barcode_folder + "/input"
    if not os.path.exists(input):
        os.makedirs(input)

    cat_files = "cat ../data/" + barcode_folder + "/*.fastq > ../data/" + barcode_folder + "/input/temp.fastq"
    subprocess.run(cat_files, shell=True)
    
    command = "NGSpeciesID --ont --consensus --medaka --primer_file ../data/primers_temp.txt --fastq ../data/" + barcode_folder + "/input/temp.fastq --outfolder ../data/" + barcode_folder + " --m 750 --s 150 --sample_size 1000"
    print(command)
    
    try:
        # Run the full command
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # Print the output and errors
        print(result.stdout)
        print(result.stderr)

        # Remove the temporary file
        remove_temp = "rm ../data/" + barcode_folder + "/input/temp.fastq"
        subprocess.run(remove_temp, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing the command: {e}")
        print(f"Output: {e.output}")
        print(f"Error: {e.stderr}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
   

    ############################################################
    
    #extract consensus fastas
    pattern = 'medaka_cl_id_*/consensus.fasta'
    directory = "../data/" + barcode_folder
    consensus_files = list_matching_files(directory, pattern)
    
    # Read each FASTA file and store the data
    headers_temp = []
    sequences_temp = []
    for consensus_file in consensus_files:
        headers, sequences = read_fasta(consensus_file)
        headers_temp.extend(headers)
        sequences_temp.extend(sequences)
        
    # Store the headers and sequences in the dictionary
    consensus_dict[sample_code] = (headers_temp, sequences_temp)
    all_headers.extend(headers_temp)
    all_sequences.extend(sequences_temp)
    
# Find the maximum number of consensus sequences for any sample
max_sequences = max(len(sequences) for headers, sequences in consensus_dict.values())

# Add columns to the dataframe for each consensus sequence
for i in range(max_sequences):
    samples[f'consensus_header_{i+1}'] = None
    samples[f'consensus_sequence_{i+1}'] = None

# Populate the dataframe with consensus sequences
for idx, row in samples.iterrows():
    sample_code = row['sample_code']
    if sample_code in consensus_dict:
        headers, sequences = consensus_dict[sample_code]
        for i, (header, sequence) in enumerate(zip(headers, sequences)):
            samples.at[idx, f'consensus_header_{i+1}'] = header
            samples.at[idx, f'consensus_sequence_{i+1}'] = sequence
            
samples.to_csv("../results/samples_results.csv")

    
