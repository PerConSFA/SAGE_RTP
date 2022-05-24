# SAGE_RTP (Relative Transcriptional Profiling)
This code for relative transcriptional profiling (RTP) is largely sourced from [DRAFTS](https://github.com/ssyim/DRAFTS) with modifications adapted to accept nucleotide sequence and Illumina read length differences for the workflows required to process DNA and RNA barcode counts related to **S**erine-recombinase **A**ssisted **G**enome **E**ngineering (SAGE). These scripts are in bash & Python.

Once barcode counts and normalized expression levels were generated, we generated R scripts to create statistical analyses of promoter expression levels and to visualize the data.

### Workflow

* Run Pre-processing.ipynb script
  * Requires path to 01_DRAFTS_process_raw.sh
  * Required bbmap installation (prerequisite for DRAFTS)
  * Outputs Paired sequence data
* Run 
