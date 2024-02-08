#!/bin/bash -l
#$ -S /bin/bash   
#$ -cwd  
#$ -j y 

# This file is designed to run bulk-rna-seq.sh
# Pass variables to this file while running in job scheduler (i.e. through '-v')
# Needs a list.txt file with each line the path to both paired reads with a semi colon between; for example: /path/to/read1;/path/to/read2
# Other variables created (R1/R2/b/c/etc.) are passed to bulk-rna-seq.sh for directory creating, file naming, etc.. purposes
# Code loops through the list.txt file and submits a job for each line, allowing for parallelization

for line in $(cat $list); do
	while IFS=';' read -ra fields; do
		R1=${fields[0]};
		R2=${fields[1]};
		b=${R1#*$name/};
		#echo "$b";
		c=${b%.$end};
		#echo "$c";
		d=${R2#*$name/};
		e=${d%.$end};
		f=${c%_*};
		echo "$f";

    # Submit bulk-rna-seq.sh job
    # Pass necessary variables
    # qsub = submit job
    # cwd = in current working directory
    # -m be = send email when job is begun and end
    # -l h_rt = max time for job to run, job will abort if takes longer than this (this command helps the job scheduler assign the job faster)
		qsub -cwd -m be -l h_rt=24:00:00 -v R1="$R1",R2="$R2",c="$c",e="$e",f="$f",out_dir="/wynton/group/singh/aakriti/Perera_files/inputs/RNA_data/$name/outputs/trim_out",genome_dir="/wynton/home/singh/aakriti/RNAseq/fastq/human/index",gtf="/wynton/home/singh/aakriti/RNAseq/fastq/human/GFF/Homo_sapiens.GRCh38.110.gtf",home="/wynton/group/singh/aakriti/Perera_files/inputs/RNA_data/$name/outputs",fasta="/wynton/home/singh/aakriti/RNAseq/fastq/human/index/Homo_sapiens.GRCh38.dna.toplevel.fa" bulk-rna-seq.sh

	done < $list
	break
done 
