# nfcore/ExoSeq Output

nfcore/ExoSeq is the new Exo-Seq Best Practice pipeline used by the [National Genomics Infrastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab](https://www.scilifelab.se/platforms/ngi/) in Stockholm, Sweden.

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [TrimGalore](#trimgalore) - adapter trimming
* [BWA-MEM](#bwa-mem) - alignment
* [Picard](#picard) - mark PCR duplicates
* [QualiMap2](#qualimap2) - quality control metrics
* [GATK](#gatk) - Genome Analysis Toolkit
	- [BAM recalibration](#bam-recalibration) - recalibrate BAM file
	- [BAM realignment](#bam-realignment) - realing BAM file around indels
	- [HaplotypeCaller](#haplotypecaller) - call variants in GVCF mode
	- [GenotypeGVCFs](#genotypegvcfs) - genotype generate GVCFs
	- [SelectVariants](#selectvariants) - select variants
	- [SNP recalibration](#snp-recalibration) - recalibrate SNPs using Omni, 1000G and dbSNP databases 
	- [Indel recalibration](#indel-recalibration) - recalibrate INDELS using the Mills golden dataset
	- [CombineVariants](#combinevariants) - combine recalibrated files
	- [VariantAnnotator](#variantannotator) - annotate variants
	- [VariantEval](#varianteval) - evaluate variants
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
The nfcore/ExoSeq BP 2.0 pipeline uses [TrimGalore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for removal of adapter contamination and trimming of low quality regions. TrimGalore uses [Cutadapt](https://github.com/marcelm/cutadapt) for adapter trimming and runs FastQC after it finishes.

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

## BWA-MEM
<!-- to do-->

BWA-MEM is an alignment algorithm for aligning sequence reads or long query sequences against a large reference genome, e.g. human. It automatically chooses between local and end-to-end alignments, supports paired-end reads. The algorithm is applicable to a wide range of sequence lengths from 70bp to a few megabases. BWA-MEM is robust to sequencing errors and shows better performance than several state-of-art read aligners to date. 

**Output directory: `results/BWAmem`**

## Picard
[Picard MarkDuplicates](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR.

The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method). For more information visit the [picard docs](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates).

**Output directory: `results/Picard_Markduplicates/metrics`**

## QualiMap2

[Qualimap 2](http://qualimap.bioinfo.cipf.es/) represents a next step in the QC analysis of HTS data. Along with comprehensive single-sample analysis of alignment data, it includes new modes that allow simultaneous processing and comparison of multiple samples.

**Output directory: `results/qualimap`**

## GATK

### BAM recalibration

The GATK BaseRecalibrator detects systematic errors in base quality scores. More information can be found here: [HATK BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)

**Output directory: `results/GATK_Recalibration`**

### BAM realignment

The GATK [RealignerTargetCreator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_indels_RealignerTargetCreator.php) is used to define intervals to target for local realignment and the GATK [IndelRealigner](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php) is used to perform local realignment of reads around indels.

**Output directory: `results/GATK_IndelRealigner`**

### HaplotypeCaller

The GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) calls germline SNPs, insertions and deletions via local re-assembly of haplotypes.

**Output directory: `results/GATK_VariantCalling`**

### GenotypeGVCFs

The GATK [GenotypeGVCFs]() performs joint genotyping on gVCF files produced by the GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php).

**Output directory: `results/GATK_GenotypeGVCFs`**

### SelectVariants

The GATK [SelectVariants](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_SelectVariants.php) selects a subset of variants from a larger callset.

**Output directory: `results/GATK_VariantSelection`**

### SNP recalibration

The GATK [VariantRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php) builds a recalibration model only for SNPs (emitting indels untouched in the output VCF) to score variant quality for filtering purposes.

**Output directory: `results/GATK_RecalibrateSNPs`**

### Indel recalibration

The GATK [VariantRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php) builds a recalibration model only for indels (emitting SNPs untouched in the output VCF) to score variant quality for filtering purposes.

**Output directory: `results/GATK_RecalibrateIndels`**

### CombineVariants

The GATK [CombineVariants](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_CombineVariants.php) combines variant records from different sources.

**Output directory: `results/GATK_CombineVariants`**

### VariantAnnotator

The GATK [VariantAnnotator](https://software.broadinstitute.org/gatk/gatkdocs/3.6-0/org_broadinstitute_gatk_tools_walkers_annotator_VariantAnnotator.php) annotates variant calls with context information.

**Output directory: `results/GATK_AnnotatedVariants`**

### VariantEval

The GATK [VariantEval](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php) general-purpose tool for variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more).

**Output directory: `results/GATK_VariantEvaluate`**

## SnpEff

[SnpEff](http://snpeff.sourceforge.net/) is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes). 

**Output directory: `results/SNPEFF_AnnotatedVariants`**

## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

**Output directory: `results/MultiQC`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info
