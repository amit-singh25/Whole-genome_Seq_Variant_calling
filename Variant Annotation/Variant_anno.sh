#!/bin/bash
#PBS -N Variant_anno
##PBS -j oe
#PBS -m e
##PBS -l file=100GB
#PBS-l mem=100G
#PBS -l walltime=50:13:59
#PBS -l nodes=1:ppn=8   
#PBS -e ./Variant_anno_$PBS_JOBID.err           # stderr file
#PBS -o ./Variant_anno_$PBS_JOBID.out
#PBS -V
#echo $PBS_JOBNAME
#echo $PBS_JOBID

name=$1
gtf=~/genome
data=~/alignment
out=~/Variant_files

vep --cache -i ${out}/${name}.vcf.gz --output_file ${out}/${name}_variant_ann.txt
