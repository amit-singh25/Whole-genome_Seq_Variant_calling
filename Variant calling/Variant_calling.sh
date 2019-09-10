


#Indexing genome fasta 
samtools faidx hs37d5.fa
bcftools mpileup -O b -o file1.chr19.bwa_raw.bcf -f hs37d5.fa file1.chr19.bwa.bam
bcftools call -m -v -o file1.chr19.bwa_raw.vcf -O v file1.chr19.bwa_raw.bcf
vcfutils.pl varFilter file1.chr19.bwa_raw.vcf >file1.chr19.bwa_final_variant.vcf
bgzip file1.chr19.bwa_final_variant.vcf
tabix -p vcf file1.chr19.bwa_final_variant.vcf.gz
