# **S**erine-recombinase **A**ssisted **G**enome **E**ngineering - Relative Transcriptional Profiling (SAGE-RTP) Toolkit
# [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7388406.svg)](https://zenodo.org/badge/latestdoi/7388406)
The SAGE-RTP toolkit leverages source code from previously reported DNA regulatory element analysis package [DRAFTS](https://doi.org/10.15252/msb.20198875), providing new modification features adapted to accept nucleotide sequences and Illumina read length differences, for data process workflows utilizing DNA and RNA barcode counts for statistical analysis of sequence promoter expression levels with data visualizations.


### Workflow

* Run Pre-processing.ipynb notebook to create folder structure compatible with DRAFTS scripts
  * Requires 'home_dir' path to Illumina DNA/RNA barcode reads
  * Requires path to [01_DRAFTS_process_raw.sh](https://github.com/ssyim/DRAFTS/blob/master/code/01_DRAFTS_process_raw.sh)
  * Required [bbmap](https://sourceforge.net/projects/bbmap/) installation (prerequisite for [DRAFTS](https://github.com/ssyim/DRAFTS))
  * Outputs Paired sequence data
* Run run_DRAFTS_DNA_v2.py to count DNA barcodes
* Run run_DRAFTS_RNA_v2.py to count RNA barcodes
* Run Build_tx_read_dataframe.ipynb notebook to calculate promoter-specific relative transcriptional activity (rta)


### Citing SAGE-RTP
Robert G. Egbert, & Joshua R. Elmore. (2022). Serine Recombinase-Assisted Genome Engineering - Relative Transcriptional Profiling (SAGE-RTP) Toolkit (1.0). Zenodo. https://doi.org/10.5281/zenodo.7388406
