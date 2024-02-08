# Workflow files for UCSF Research

## bulk-rna-seq.sh and run-bulk-rna-seq.sh
Designed to run together. run-bulk-rna-seq.sh loops through a text file and submits bulk-rna-seq.sh jobs for parallelization. For bulk RNA-seq analysis.

## featurecounts.sh
Creates counts matrix from directory of bams. Bams outputs come from bulk-rna-seq.sh.

## run-chipseq.sh and chipseq.sh
Designed to run together. run-chipseq.sh loops through a text file and submits chipseq.sh jobs for parallelization. For ChIP-seq analysis.

## pseudobulk.sh
Runs pseudobulk analysis on single-cell data using Seurat. Post-pseudobulk file is run through DESeq pipeline for QC purposes (compare with known bulk datasets).

## macs3-pipeline.sh
Runs the macs package on ChIP-seq datasets to get peaks for downstream analysis.

## idr.sh
Runs the IDR package for QC purposes, examining reproducibility between replicates in ChIP-seq datasets. Uses post-macs3 narrowPeak files. 

## homer.sh 
Runs the homer package for motif analysis for ChIP-seq datasets.

## doubletfinder.Rmd
Runs the doubletfinder package for single-cell data. Finds doublets that can be removed as a QC step.

## soupx.Rmd
Runs the soupX package for single-cell data. Detects and removes ambient RNA as a QC step.
