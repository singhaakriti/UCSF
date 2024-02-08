#!/bin/bash -l
#$ -S /bin/bash   
#$ -cwd  
#$ -j y 

# This file creates counts matrix from bam files 
# Example job submit command: 
# qsub -cwd -m be -l h_rt=01:00:00 -v gtf="/path/to/gtf",countsFile="/path/to/counts/file/being/created/name.csv",inputs="/path/to/directory/of/bams/*.bam" featurecounts.sh
# -p = paired
# -a = genome annotation file 

featureCounts -p -a $gtf -o $countsFile $inputs 
