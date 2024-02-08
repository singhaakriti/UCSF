#!/bin/bash -l
#$ -S /bin/bash   
#$ -cwd  
#$ -j y 

# This file runs the homer program to find motifs for ChIP-seq analysis  

cd ~/homer

cat peaks/$R1 | grep -n ^chr # tells you at what point the information homer needs starts 
cat peaks/"$c" | tail -n +30 | cut -f 1-3,10 > peaks/homer_"$c"
findMotifsGenome.pl peaks/homer_"$c" human_genome "$c"_OUT/ -size 200
