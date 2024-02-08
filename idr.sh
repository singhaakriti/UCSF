#!/bin/bash -l
#$ -S /bin/bash   
#$ -cwd  
#$ -j y 

# This file runs IDR on macs narrowPeak outputs as a QC step

# Sort peaks first
sort -k8,8nr copied.NA_peaks.narrowPeak > sorted.copied.NA_peaks.narrowPeak

# Run IDR
idr --samples /path/to/replicate1 /path/to/replicate2 \
--input-file-type narrowPeak \
--rank p.value \
--output-file KP4_TFE3 \
--plot \
--log-output-file kp4_tfe3.idr.log
