#!/bin/bash
# Zane Libke
# CREATED : August 17, 2024
# LAST EDIT: September 5, 2024
# script to take sample and blasthit fasta and create a tree from them

#conda install mafft
#conda install seqkit
#conda install iqtree

# Get the input FASTA file and BLAST XML file from command-line arguments
sample_fasta="$1"
#min_src="$2" # allow user to decide min supporting reads to discard consensuses
#retrieve_inat="$3" # retrieve metadata from inat flag

#Names of various necessary files
base_name=$(basename "$sample_fasta" .fasta)
filtered_fasta="../data/${base_name}_filtered.fasta"
hits_fasta="../data/${base_name}_filtered_blasthits.fasta"
samples_aligned="../results/${base_name}_samples_alignment.fasta"
sampleshits_aligned="../results/${base_name}_sampleshits_alignment.fasta"

echo "$filtered_fasta"
echo "now basename next:"
echo "$base_name"
echo "now blasthits fasta:"
echo "$hits_fasta"

#filter sequences for >=20 reads
python3 filter_20.py "$sample_fasta"

# Blast filtered sequences, store hits
if [ -e "$hits_fasta" ] && [ -s "$hits_fasta" ]; then
    echo "Blast hits file '$hits_fasta' exists and is not empty, skipping blast query..."
else
    echo "Blast hits file '$hits_fasta' either does not exist or is empty, querying blast..."
    python3 query_blast.py "$filtered_fasta"
fi

#remove duplicates from blasthits file
python3 remove_duplicates.py "$hits_fasta"

mafft --adjustdirectionaccurately --quiet "$filtered_fasta" > ../data/temp.fasta
mafft --adjustdirectionaccurately --quiet ../data/temp_nodup.fasta > ../data/temp_blast.fasta
cat ../data/temp.fasta ../data/temp_blast.fasta > ../data/temp_2.fasta
mafft --adjustdirectionaccurately --quiet ../data/temp_2.fasta > ../data/temp_3.fasta

#iqtree -s temp_3.fasta -st DNA -m JC69 -nt AUTO -bb 1000

#cleanup 
mv ../data/temp.fasta "$samples_aligned"
rm ../data/temp_2.fasta
mv ../data/temp_3.fasta "$sampleshits_aligned"
rm ../data/temp_blast.fasta
rm ../data/temp_nodup.fasta

