#!/bin/bash -l
#$ -S /bin/bash   
#$ -cwd  
#$ -j y 

# This file processes human RNA seq data, ***created to run with run-bulk-rna-seq.sh***

module load CBI
module load fastqc/0.12.1 
module load star/2.7.10b
module load samtools/1.18


# Trim with FASTP 
fastqc $R1 $R2
~/fastp -i $R1 -o $out_dir/OUT_"$c".fq.gz -I $R2 -O $out_dir/OUT_"$e".fq.gz -j $c.json -h $c.html

# ---- CREATE GENOME ----  
# Only need to run once 
# STAR --runMode genomeGenerate --genomeDir $genome_dir \
#             --genomeFastaFiles $fasta \
#             --sjdbGTFtagExonParentTranscript Parent $gtf \
#             --genomeSAindexNbases 10
# outputs index with a Genome, genomeparameters.txt, etc

# ---- ALIGN ----
STAR --runThreadN 2 --genomeDir $genome_dir \
            --readFilesIn $out_dir/OUT_"$c".fq.gz $out_dir/OUT_"$e".fq.gz \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix $home/"$c"_read_alignment_star/"$c"_

# Convert bam
samtools view -q255 -b $home/"$c"_read_alignment_star/"$c"_Aligned.sortedByCoord.out.bam > $home/"$c"_read_alignment_star/"$c".bam

# ----- COUNTS --------
featureCounts -p -a $gtf -o $home/"$c"_read_alignment_star/"$c".featureCounts $home/"$c"_read_alignment_star/"$c".bam
cp $home/"$c"_read_alignment_star/"$c".bam $home/bams/"$c".bam
# #featurecounts output.featureCounts input.bam
# # -g $7

# Index
samtools index $home/"$c"_read_alignment_star/"$c".bam
# # gives me a .bai file (indexed bam)

# Bigwig
bamCoverage -b $home/"$c"_read_alignment_star/"$c".bam \
            -o $home/"$c"_read_alignment_star/"$c".bw
            
# ---- LIST JOB ID ----
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"
