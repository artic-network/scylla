# scylla

This pipeline performs classification and taxonomic read binning for metagenomic reads

For each sample provided, the following analyses are performed:
 - Input reads are QC filtered with fastp by length, quality, ambiguous bases and adaptors and polyX runs are trimmed (default settings)
 - Read length and quality stats are collected and single read files are checked for duplicate read names
 - Reads are classified using Kraken2 and the PlusPF database to get per-read taxon classifications
 - For species level taxa meeting the threshold percentage and counts (higher threshold for bacteria, lower for viruses etc), a file of binned taxon reads is extracted
 - The human read count in classifications is checked to make sure it is less than the threshold (1000). It is not possible to extract human-classified reads in any output.
 - Reads are screened for high consequence infectious diseases (HCID) using both classified counts and read mapping to a reference genome panel (including close confounders) and checking the proportion of the genome covered. 
 - Reads are screened for the specified spike ins and a summary collected.
 - Separate files of the unclassified read fraction, the viral+unclassified read and all non-human reads are extracted.
 - A single sample report is generated

This pipeline can optionally perform a reclassification of the viral+unclassified read fraction with a custom viral database.

There is the option to perform _de novo_ viral assembly and classification using virbot or genomad.

This pipeline can be run on a single fastq/pair of fastq files, or a directory of fastq files (which will be concatenated), or a demultiplexed directory of barcode subdirectories. 

On personal computers we recommend running with the `--local` flag to ensure more reasonable resource requirements. 

### Example command
These examples use the `--local` flag to set resource requirements for a laptop. This will use a smaller kraken database, and therefore when working on a cluster we recommend you remove this flag.

1. Run in mSCAPE ingest mode (assumes input is a single sample, allows single/paired fastq file input).
```
nextflow run main.nf --fastq test/test_data/barcode01/barcode01.fq.gz -profile docker --local
```
or
```
nextflow run main.nf --fastq1 test/test_data/illumina/barcode02_R1.fq.gz --fastq2 test/test_data/illumina/barcode02_R2.fq.gz --paired -profile docker --local```
```

2. Run on a demultiplexed run directory (providing a `run_dir`). This must either contain a subdirectory per barcode, or pairs of files with names in the format `*_R{1,2}*.f*q*`.
```
nextflow run main.nf --run_dir test/test_data [--paired] -profile docker --local
```

3. Run a module independently of the main workflow. Supported modules are [preprocess, qc_checks, centrifuge_classification, kraken_classification, sourmash_classification, check_hcid_status, check_spike_status, extract_all, classify_novel_viruses].
```
nextflow run main.nf --module preprocess --fastq test/test_data/barcode01/barcode01.fq.gz -profile docker --local
```