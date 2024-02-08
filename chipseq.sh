#!/bin/bash -l
#$ -S /bin/bash   
#$ -cwd  
#$ -j y 

# This file processes  ChIP seq data, created to run with run_chipseq.sh 
# This is for single-ended data

# Load modules
module load CBI
module load bowtie2/2.5.1
module load samtools
#module load macs3

#TRIM FILES
~/TrimGalore-0.6.10/trim_galore --fastqc \
					--fastqc_args "--outdir $out_dir/fastqc" \
					--gzip \
					--output_dir $out_dir/"$c" \
					$R1

#ALIGN READS TO REFERENCE GENOME: BOWTIE2
bowtie2 -q --no-unal -x $genome_dir \
					-U $out_dir/"$c"/"$c"_trimmed.fq.gz \
					-S $out_dir/"$c"/"$c".sam

# #MAKE BAM FILES, MAKE INDEX
# #Remove mitochondrial dna, or unassigned/random reads 
sed '/chrM/d;/random/d;/chrUn/d' $out_dir/"$c"/"$c".sam > $out_dir/"$c"/"$c".filtered.sam
samtools view -b $out_dir/"$c"/"$c".filtered.sam | samtools sort -o $out_dir/"$c"/"$c".bam
samtools view -c $out_dir/"$c"/"$c".bam
samtools index $out_dir/"$c"/"$c".bam

#bamCoverage
bamCoverage -b $out_dir/"$c"/"$c".bam \
					-o $out_dir/"$c"/"$c".bw

# samtools view -b -f 4 K_InpB_3.bam > unmapped_K_InpB_3.bam

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
