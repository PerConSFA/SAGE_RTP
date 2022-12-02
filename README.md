# Serine Recombinase-Assisted Genome Engineering - Relative Transcriptional Profiling (SAGE-RTP)
# [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7388406.svg)](https://zenodo.org/badge/latestdoi/7388406)
This code for relative transcriptional profiling (RTP) is largely sourced from [DRAFTS](https://github.com/ssyim/DRAFTS) with modifications adapted to accept nucleotide sequence and Illumina read length differences for the workflows required to process DNA and RNA barcode counts related to **S**erine-recombinase **A**ssisted **G**enome **E**ngineering (SAGE). These scripts are in bash & Python.

Once barcode counts and normalized expression levels were generated, we generated R scripts to create statistical analyses of promoter expression levels and to visualize the data.

### Workflow

* Run Pre-processing.ipynb notebook to create folder structure compatible with DRAFTS scripts
  * Requires 'home_dir' path to Illumina DNA/RNA barcode reads
  * Requires path to 01_DRAFTS_process_raw.sh
  * Required bbmap installation (prerequisite for DRAFTS)
  * Outputs Paired sequence data
* Run run_DRAFTS_DNA_v2.py to count DNA barcodes
* Run run_DRAFTS_RNA_v2.py to count RNA barcodes
* Run Build_tx_read_dataframe.ipynb notebook to calculate promoter-specific relative transcriptional activity (rta)
