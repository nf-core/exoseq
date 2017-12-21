## Processes

The current implementation of the pipeline relies on Singularity containers for each pipeline step to allow easier orchestration. 

Current processes:

* FastQC for initial QC of data
* TrimGalore! for trimming of reads
* BWA mem for read mapping
* Samtools sort for read sorting
* Picard MarkDuplicates for read deduplication
* GATK BaseRecalibrator and PrintReads for read recalibration
* GATK RealignerTargetCreator and IndelRealigner for Indel realignment (probably optional soon)
* QualiMap2 for determining metrics of the datasets investigated
* GATK HaplotypeCaller in GVCF mode to create GVCFS
* GATK GenotypeVCFs to create a multi-sample GVCF file
* GATK VariantSelect to extract both SNPs and Indels in separate files to perform recalibration of found variants
* GATK CombineVariants to merge recalibrated files again
* SNPeff Annotate to annotate the variants accordingly
* GATK Annotator to annotate the variants with the SNPEff Profile created
* GATK VariantEvaluate to determine important metrics for downstream analysis (TITV and others)
* MultiQC to collect results from all steps 