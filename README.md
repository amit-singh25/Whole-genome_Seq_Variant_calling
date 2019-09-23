# Whole genome Sequencing and Variant calling and Variant Annotation



# Table of content
* [Software Requirement](#QRequired-packages)
* [Quality control](#Quality-control)
* [Mapping and Filtering](#Alignment-filtering )
* [Variant Calling](#Variant-calling)
* [Variant Annotation](#Variant-calling)
* [Variant filtering](#Variant-filtering)
* [Visualization](#Visualization)

## Software Requirement
Required packages for processing the ATAC-seq pipeline.

Following softwear can be install in the cluster either from source code or from [conda](https://conda.io/en/latest/) platform 

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Flexbar](https://github.com/seqan/flexbar)

[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

[samtools](http://samtools.sourceforge.net/)

[PICARD](https://broadinstitute.github.io/picard)

[deepTools](https://deeptools.readthedocs.io/en/develop/)

[Bcftools](https://samtools.github.io/bcftools/bcftools.html)


R environment 
## Quality control

### FastQC
It is generally a good idea to generate some quality metrics for raw sequence data using [FastQC]( https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). 

Quality-based trimming as well as Adapter removal can be done in [Flexbar](https://github.com/seqan/flexbar)

## Mapping reads onto a reference genome 
The next step is to align the reads to a reference genome. There are many programs available to perform the alignment. The high efficiency and accuracy alines without allowing large gaps, such as splice junctions is [BWA](http://bio-bwa.sourceforge.net/bwa.shtml). 

#### Genome indexing
The genome indexes can build using bwa index
bwa index [-a bwtsw|is] index_prefix reference.fasta
bwa index -p hg19bwaidx -a bwtsw hg19.fa
-p index name (change this to whatever you want)
-a index algorithm (bwtsw for long genomes and is for short genomes)

#### Align to Reference Genome
bwa mem -M 
-M to flag shorter split hits as secondary. This is optional for Picard compatibility as MarkDuplicates can directly process BWA's alignment

#### Remove duplicated mapped reads
Sequenceing error occurs sometimes in sample or library preparation. When reads come from the exact same input DNA template and accumulate at the same start position on the reference genome due to the PCR dulication. Any sequencing error could lead to artefacts in the downstream variant calling. It is impossible to distinguish the true DNA materials and the PCR artifacts. To decrease this harmful effect of duplicates prior the Picard tol of [MarkDuplicates](https://broadinstitute.github.io/picard/) is useful. The Picard MarkDuplicates uses the start coordinates and orientations of both reads of a read pair. Based on the identical 5’mapping coordinates it discards all duplicates with the exception of the “best” copy.


#### Validating BAM files
It is best practice to validate the BAM file, to make sure there were not issues or mistakes with previous analyses. This is done with [ValidateSamFile](https://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile).

## Variant calling

Nucleotide difference vs.some reference at a given position in an individual genome or transcriptome usually referred to as a Single Nucleotide Polymorphism (SNP). The variant call is usually accompanied by an estimate of variant frequency and some measure of confidence. There are a number of tools available for variant calling, however, bcftools is one of the widely used software. Following steps were performed to for variant calling.

1. Calculate the read coverage of positions in the genome
2. Detect the single nucleotide polymorphisms (SNPs)
3. Filter and report the SNP variants in variant calling format (VCF)
4. Vcf indexing and merging to one file

## Variant Annotation

Annotate the SNP variants connecting to different functional aspects using variant effect predictor[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html). 
The VCF [file header information](https://software.broadinstitute.org/gatk/documentation/article.php?id=1268)


## Genetic variants filtering

## Visualization

#### Creating browser tracks

Create a bigWig file for visualizing the peak covrage using bamCoverage in deepTools. 
An alternative visualization tool is [Integrative Genomics Viewer](https://software.broadinstitute.org/software/igv/).The Peak files can be loaded directly (File → Load from File). Viewing BAM files with IGV requires sorted (by coordinate) and indexed using SAMtools.
For making plot BAM file can be converted to bed (bam to bed) using [bedtools](https://bedtools.readthedocs.io/en/latest/content/tools/bamtobed.html) and load to IGV.  

# BAM and BED file Manipulation 

 ##### First identify the depth at each locus from a bam file.

``` samtools depth test.bam > test.coverage ```

###### genomeCoverageBed
``` genomeCoverageBed -ibam test.bam -g genome.fasta > coverage.txt```

##### Select the coverage (depth) by locus for each chromosome 

##### To select the coverage for a particular chromosome 
``` awk '$1 == 10 {print $0}' test.coverage > chr10_test.coverage```

##### If the chrosomosome has string then

``` awk '$1 == "chr10" {print $0}' test.coverage > chr10_test.coverage ```

##### samtools flagstat test.bam 

```  samtools stats test.bam |grep ^SN | cut -f 2-  ```

```  samtools view test_sorted.bam | wc -l  ```

##### Convert BAM To BED file 

``` bedtools bamtobed -i ${out}/${name}_sort.bam >${out}/${name}_sort.bed ```

##### Convert BAM To bigWig file 

`` bamCoverage -b ${out}/${name}_sort.bam -o ${out}/${name}_sort.bigWig `` 

##### Convert BAM file to bedgraph

`` bam2bed < ${out}/${name}_sort.bam | cut -f1-3,5 > ${out}/${name}_sort.bg `` 

##### Convert GTF file to bed file

``  gtf2bed < foo.gtf > foo.bed `` 
##### bedtools "intersect"
It compares two or more BED/BAM/VCF/GFF files and identifies all the regions in the gemome where the features in the two files overlap

``  bedtools intersect -a test1.bed -b test2.bed | head -5 `` 

##### Reporting the original feature in each file of intersect

`` bedtools intersect -a test1.bed -b test2.bed  -wa -wb | head -15 ``

##### Reporting How many base pairs of overlap between two bed files

`` bedtools intersect -a test1.bed -b test2.bed  -wo | head -10 `` 

##### Counting the number of overlapping features in two bed files

`` bedtools intersect -a test1.bed -b test2.bed  -c | head -10 ``

##### Find features that DO NOT overlap between two bed files 

`` bedtools intersect -a test1.bed -b test2.bed -v >final.bed ``
##### Intersecting multiple bed files at a time

``  bedtools intersect -a test1.bed -b test2.bed test3.bed test4.bed -sorted >final.bed `` 
`` bedtools intersect -a test1.bed -b test2.bed test3.bed test4.bed -sorted -wa -wb >final.bed ``
`` bedtools intersect -a test1.bed -b test2.bed test3.bed test4.bed -sorted -wa -wb -names test test2 chrom ``

##### Bed file sorted
`` sort -k1,1 -k2,2n test.bed > test.sort.bed ``

##### Get genomecovrage 
``  bedtools genomecov -i test.bed -g genome.txt `` 

##### Merging features that are close to one another in a bed file
d is the distance in the. For example, to merge features that are no more than 100bp apart, one would run:
``  bedtools merge -i test.bed -d 100 -c 1 -o count | head -10 `` 
##### Jaccard statistic to meaure the similarity of two datasets in the bed file
`` bedtools jaccard  -a test1.bed  -b test2.bed ``

###### CollectAlignmentSummaryMetrics
 ``` picard CollectAlignmentSummaryMetrics \ ```
          ``` R=genome.fasta \ ```
          ``` I=input.bam \ ```
          ``` O=output.txt  ```


##### Plot the data in R this coverage file 
 ``` test.chr10 <- read.table("~/data/test.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE) ```
``` library(reshape) ```

``` test.chr10<-rename(test.chr10,c(V1="Chr", V2="locus", V3="depth")) # renames the header ```
``` plot(test.chr10$locus,test.chr10$depth) ```
``` library(lattice, pos=10) xyplot(depth ~ locus, type="p", pch=16, auto.key=list(border=TRUE), par.settings=simpleTheme(pch=16), ```  ```scales=list(x=list(relation='same'), y=list(relation='same')), data=test.chr10, main="depth by locus - Chr10)") ```

http://rstudio-pubs-static.s3.amazonaws.com/334574_1329d2c1f7274328a6309cf61a43feb4.html










