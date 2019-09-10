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

## mapping reads onto a reference genome 
#### Genome indexing

For many model organisms, the genome and pre-built reference indexes are available from [iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). Bowtie2 indexes can be made directly from [FASTA](ftp://ftp.ensembl.org/pub/release-97/fasta/)genome file using bowtie2-buid. 

#### Alignment

The next step is to align the reads to a reference genome. There are many programs available to perform the alignment. Two of the most popular are [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) and [Bowtie2](http://bowtie-bio.sourceforge.net/index.shtml). Here focus more on Bowtie2.

#### Remove duplicated mapped reads

During library preparation procedure some PCR artifacts may arise that might interfere with the biological signal of interest 
Therefore, they should be removed as part of the analysis pipeline before peak calling. 
One commonly used program for removing PCR duplicates is Picard’s [MarkDuplicates](https://broadinstitute.github.io/picard/). Removal of PCR duplicates may not necessary in Chip seq data. To undertand the flag number and [samtool format](https://www.samformat.info/sam-format-flag) look her. 

#### Quality control of mapped reads


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

