#!/bin/bash
#PBS -N Alignment
##PBS -j oe
#PBS -m e
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -e ./Alignment_$PBS_JOBID.err           # stderr file
#PBS -o ./Alignment_$PBS_JOBID.out
#PBS -V
#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
index=~/genome
data=~/Fasta_file
out=~/alignment

# bwa index [-a bwtsw|is] index_prefix reference.fasta
##downlod refernce genome from ensemble
#bwa index -p hg38index -a bwtsw hg38.fa

# -p index name (change this to whatever you want)
# -a index algorithm (bwtsw for long genomes and is for short genomes)
####alignment 
bwa mem -m ${index}/hg38index ${data}/${name}_R1.fastq.gz ${data}/${name}_R2.fastq.gz | samtools view -buS - | samtools sort - -o ${out}/${name}_map_pe.sorted.bam
###remove PCR duplicate 


