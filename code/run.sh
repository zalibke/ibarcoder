#!/bin/bash

#bash script to run simple ibarcoder pipeline on example dataset
#Zane Libke - 5 September 2024

# Step 1 : input samples.txt and primers_db.txt

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <samples.txt> <primers_db.txt>"
    exit 1
fi

# Step 1: Input samples.txt and primers_db.txt
samples=$1
primers_db=$2

# Step 2: run processor.py to create consensus sequences for all samples
python3 processor.py "$samples" "$primers_db"
echo "consensus sequences produced, saved to ../results/consensus_sequences.fasta"

# Step 3: run align_maketree to discard consensus <20, blast, and create mafft alignment
bash align_maketree.sh ../results/consensus_sequences.fasta
echo "blasthits and alignment produced, paste ../results/*_sampleshits_alignment.fasta into mafft.cbrc.jp to view tree"


