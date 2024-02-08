#!/bin/bash -l
#$ -S /bin/bash   
#$ -cwd  
#$ -j y 

# This pipeline is to run macs3 on treated and input .bam  files from ChIP seq
# In the .txt file, the first path should point to the treated.bam file [R1] and the second path should point to the input / control .bam file [R2]

# -f = file type
# -g = genome size
# -q = q value
# -n = name
# outdir = output directory


macs3 callpeak -t /path/to/treated.bam -c /path/to/control.bam -f BAM -g 2.7e9 -q 0.05 -n name --outdir /wynton/home/singh/aakriti/ChIPseq/KP4_TFE3

