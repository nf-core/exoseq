# NGI-ExoSeq Output

NGI-ExoSeq is the new Exo-Seq Best Practice pipeline used by the [National Genomics Infrastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab](https://www.scilifelab.se/platforms/ngi/) in Stockholm, Sweden.

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [TrimGalore](#trimgalore) - adapter trimming
* [BWA](#bwa) - alignment
* [Picard](#picard) - mark PCR duplicates
* [QualiMap2](#qualimap2) - quality control metrics
* [GATK](#gatk) - Genome Analysis Toolkit
	- [BAM recalibration](#bam-recal) - recalibrate BAM file
	- [BAM realignment](#bam-realign) - realing BAM file around indels
	- [HaplotypeCaller](#haplotypecaller) - call variants in GVCF mode
	- [GenotypeGVCFs](#gvcf) - genotype generate GVCFs
	- [SelectVariants](#selectvariants) - select variants
	- [SNP recalibration](#snp-recal) - recalibrate SNPs using Omni, 1000G and dbSNP databases 
	- [Indel recalibration](#indel-recal) - recalibrate INDELS using the Mills golden dataset
	- [CombineVariants](#combine-var) - combine recalibrated files
	- [VariantAnnotator](#annotate-var) - annotate variants
	- [VariantEval](#eval-var) - evaluate variants
* [SnpEff](#snpeff) - annotate variants
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## TrimGalore
The NGI-ExoSeq BP 2.0 pipeline uses [TrimGalore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for removal of adapter contamination and trimming of low quality regions. TrimGalore uses [Cutadapt](https://github.com/marcelm/cutadapt) for adapter trimming and runs FastQC after it finishes.

MultiQC reports the percentage of bases removed by TrimGalore in the _General Statistics_ table, along with a line plot showing where reads were trimmed.

**Output directory: `results/trim_galore`**

Contains FastQ files with quality and adapter trimmed reads for each sample, along with a log file describing the trimming.

* `sample_val_1.fq.gz`, `sample_val_2.fq.gz`
  * Trimmed FastQ data, reads 1 and 2.
  * NB: Only saved if `--saveTrimmed` has been specified.
* `logs/sample_val_1.fq.gz_trimming_report.txt`
  * Trimming report (describes which parameters that were used)
* `FastQC/sample_val_1_fastqc.zip`
  * FastQC report for trimmed reads

Single-end data will have slightly different file names and only one FastQ file per sample.

## BWA
<!-- to do-->

BWA is a read aligner designed for RNA sequencing.  BWA stands for Spliced Transcripts Alignment to a Reference, it produces results comparable to TopHat (the aligned previously used by NGI for RNA alignments) but is much faster.

The BWA section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as _Uniquely mapped_ and few _Unmapped_ reads.

![BWA](images/BWA_alignment_plot.png)

**Output directory: `results/BWA`**

* `Sample_Aligned.sortedByCoord.out.bam`
  * The aligned BAM file
* `Sample_Log.final.out`
  * The BWA alignment report, contains mapping results summary
* `Sample_Log.out` and `Sample_Log.progress.out`
  * BWA log files, containing a lot of detailed information about the run. Typically only useful for debugging purposes.
* `Sample_SJ.out.tab`
  * Filtered splice junctions detected in the mapping


## Picard
[Picard MarkDuplicates](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR.

The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).

The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, please see the following [blog post](https://software.broadinstitute.org/gatk/blog?id=7019) for additional information.

MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads. If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.

**Output directory: `results/picard`**

* `to_do`

Picard MarkDuplicates documentation: [picard docs](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

## QualiMap2

<!-- to do -->
## GATK

### BAM recalibration

### BAM realignment

### HaplotypeCaller

### GenotypeGVCFs

### SelectVariants

### SNP recalibration

### Indel recalibration

### CombineVariants

### VariantAnnotator

### VariantEval
<!-- to do -->

## SnpEff




## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
