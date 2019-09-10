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

