# scylla

This pipeline performs classification and taxonomic read binning for metagenomic reads

For each sample provided, the following analyses are performed:
 - Input reads are QC filtered with fastp by length, quality, ambiguous bases and adaptors and polyX runs are trimmed (default settings)
 - Reads are classified using Kraken2 and the PlusPF database to get per-read taxon classifications
 - For species level taxa meeting the threshold percentage and counts (higher threshold for bacteria, lower for viruses etc), a file of binned taxon reads is extracted
 - The human read count in classifications is checked to make sure it is less than the threshold (1000). It is not possible to extract human-classified reads in any output.
 - Reads are screened for high consequence infectious diseases (HCID) using both classified counts and read mapping to a reference genome panel (including close confounders) and checking the proportion of the genome covered. 
 - Separate files of the unclassified read fraction, the viral+unclassified read and all non-human reads are extracted.
 - A single sample report is generated

This pipeline can optionally perform _de novo_ viral classification using virbot or genomad.

This pipeline can be run on a single fastq/pair of fastq files, or a directory of fastq files (which will be concatenated), or a demultiplexed directory of barcode subdirectories. 

On personal computers we recommend running with the `--local` flag to ensure more reasonable resource requirements. 

### Example command
```
nextflow run main.nf --run_dir test/test_data -profile docker --local
```
