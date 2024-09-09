# R script to get alignment seeds
# Zane Libke - 9 September 2024

library(dplyr)
library(Biostrings)

# Read the CSV file
df <- read.csv("../data/samples.csv")

# Extract unique combinations of "Order", "primerF", and "primerR"
unique_combinations <- df %>% select(Order, primerF, primerR) %>% distinct()

# Define the FASTA file and output file
fasta_file <- "../data/seeds_db.fasta" #seeds with header of format order_primerF_primerR
output_fasta <- "../data/temp_seeds.fasta"

# Read the sequences from the FASTA file
sequences <- readDNAStringSet(fasta_file)

# Open the output file
output_handle <- file(output_fasta, "w")

# Loop over each unique combination
for (i in 1:nrow(unique_combinations)) {
  order <- unique_combinations$Order[i]
  primerF <- unique_combinations$primerF[i]
  primerR <- unique_combinations$primerR[i]
  
  # Find matching sequences based on the combination
  matching_seqs <- sequences[grepl(paste0("^", order, "_", primerF, "_", primerR), names(sequences))] 
   
  # Write matching sequences to the output file
  writeXStringSet(matching_seqs, output_handle, format = "fasta") #CHECK THIS!!
}

# Close the output file
close(output_handle)

cat("Alignment seeds written to ", output_fasta, "\n")
