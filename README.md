# Workflow files for UCSF Research

## bulk-rna-seq.sh and run-bulk-rna-seq.sh
Designed to run together. run-bulk-rna-seq.sh loops through a text file and submits bulk-rna-seq.sh jobs for parallelization. For bulk RNA-seq analysis

## featurecounts.sh
Creates counts matrix from directory of bams. Bams outputs come from bulk-rna-seq.sh

## run-chipseq.sh and chipseq.sh
Designed to run together. run-chipseq.sh loops through a text file and submits chipseq.sh jobs for parallelization. For ChIP-seq analysis
