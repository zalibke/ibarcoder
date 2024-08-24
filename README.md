---

# Seamless Bioinformatics Pipeline for Nanopore Amplicon Sequencing

**Version** Alpha

> **Note:** Do not attempt to use `#3 "processor.py"` yet - still working out some bugs.  
> `#4 align_maketree.sh` should work fine. Usage:
> ```bash
> bash align_maketree.sh samples.fasta
> ```

## Pipeline Roadmap

1. **Basecall (POD5 -> FASTQ):**  
   - Currently not integrated, testing has been carried out on data basecalled with Dorado.

2. **Demultiplex:**  
   - Currently not integrated, testing has been carried out on data demuxed with Dorado.  
   - *Note*: The current format requires demuxed reads to be placed into their own directories with the syntax "barcode01, barcode02, up to barcode96".

3. **Consensus Generation:**  
   - *Script*: `code/processor.py`  
   - processor runs NGspeciesID on all barcodes and extracts output consensus sequences into a convenient dataframe with metadata.  
   - *Input*:
     - **Demultiplexed FASTQ files**: Should be placed in the "data" directory. Each barcode should have its own directory in the format "barcode01, barcode02, up to barcode96".
     - **`samples.txt`**: A database of samples. An ideal file contains at least the following columns:
       - `index (int 1:96)`
       - `sample_code (chr)` - Sample ID from field collection, etc.
       - `primerF (chr)` - Forward primer used to amplify.
       - `primerR (chr)` - Reverse primer used to amplify.
       - `genus_prior (chr)` - Preliminary genus ID.
       - `species_prior (chr)` - Preliminary species ID.
       - `lat (float)` - Latitude of sample collection.
       - `long (float)` - Longitude of sample collection.
       - `acc (float)` - Accuracy of coordinates (in meters).
       - `inat (chr)` - Link to iNaturalist observation. This could also consider an iNat ID.
       - `collection_date (chr)` - Date of collection in a format like `2024:06:28:17:21` or similar.
     - **Retrieve inat? - Optional - Retrieve data from iNaturalist observation!**:  
       - If this flag is TRUE, the program will call the script `retrieve_inat.py`, which takes only an iNat observation ID or link and queries the iNat API to retrieve coordinates, date/time, preliminary genus/species ID, and associated notes.

4. **Filter, Blast, and Make Tree:**  
   - `code/align_maketree.sh`  
   - This script should be run once you have a FASTA of consensus sequences. It filters out any consensus sequences with <20 supporting reads and saves them to a FASTA called `discard.fasta`. It then blasts the passed sequences via the BLAST API and downloads the top hit for each sequence. Once this has been done for every sequence, it removes duplicate BLAST hits and aligns both the samples and BLAST hits with MAFFT, accounting for reverse complements. This alignment can then be uploaded to MAFFT's online GUI to visualize the tree.
   - *Input*:
     - `samples.fasta` - A FASTA file with the samples you want to make a tree of together.

--- 