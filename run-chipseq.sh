#!/bin/bash -l
#$ -S /bin/bash   
#$ -cwd  
#$ -j y 

# Designed to run with chipseq.sh
# Loops through a text file and submits job commands for chipseq.sh for parallelization
# Example command:
# qsub -cwd -m be -l h_rt=01:00:00 -v list="/path/to/list.txt" run-chipseq.sh

# This is for single end 

for line in $(cat $list); do
	# while IFS='\n' read -ra fields; do
	R1=${line};
	b=${R1#*KP4_TFE3/};
	c=${b%.fq.gz};
	echo "$c";

		qsub -cwd -l h_rt=24:00:00 -v R1="$R1",R2="$R2",c="$c",e="$e",out_dir="/wynton/home/singh/aakriti/ChIPseq/published/OUT",genome_dir="/wynton/home/singh/aakriti/ChIPseq/index/GRCh38",home="/wynton/home/singh/aakriti/ChIPseq/published" chipseq.sh

	# done < $list
	# break
done 
