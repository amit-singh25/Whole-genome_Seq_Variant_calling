#!/bin/bash
#PBS -N Variant_calling
##PBS -j oe
#PBS -m e
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -e ./Variant_calling_$PBS_JOBID.err           # stderr file
#PBS -o ./Variant_calling_$PBS_JOBID.out
#PBS -V
#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
gtf=~/genome
data=~/alignment
out=~/Variant_files


#Indexing genome fasta 
samtools faidx hs37d5.fa
bcftools mpileup -O b -o ${out}/${name}.bwa_raw.bcf -f ${gtf}/hs37d5.fa ${data}/${name}.bwa.bam
bcftools call -m -v -o ${out}/${name}.bwa_raw.vcf -O v ${out}/${name}.bwa_raw.bcf
vcfutils.pl varFilter ${out}/${name}.bwa_raw.vcf >${out}/${name}.bwa_final_variant.vcf
bgzip ${out}/${name}.bwa_final_variant.vcf
tabix -p vcf ${out}/${name}.bwa_final_variant.vcf.gz


