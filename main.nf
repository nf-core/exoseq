#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

========================================================================================
               QBIC - E X O S E Q    B E S T    P R A C T I C E
========================================================================================
 This is based on previous work at NGI, see webpage for details
 #### Homepage / Documentation
 https://github.com/apeltzer/QBIC-ExoSeq
 #### Authors
 Senthilkumar Panneerselvam @senthil10 <senthilkumar.panneerselvam@scilifelab.se>
 Phil Ewels @ewels <phil.ewels@scilifelab.se>
 Alex Peltzer @alex_peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
----------------------------------------------------------------------------------------
Developed based on GATK's best practise, takes set of FASTQ files and performs:
 - alignment (BWA)
 - recalibration (GATK)
 - realignment (GATK)
 - variant calling (GATK)
 - vairant evaluation (SnpEff)
*/

// Package version
version = '1.0'

// Help message
helpMessage = """
===============================================================================
QBIC-ExoSeq : Exome/Targeted sequence capture best practice analysis v${version}
===============================================================================

Usage: nextflow apeltzer/QBIC-ExoSeq --reads 'P1234*_R{1,2}.fastq.gz' --genome GRCh37

This is a typical usage where the required parameters (with no defaults) were
given. The all avialable paramaters are listed below based on category

Required parameters:
--reads                        Absolute path to project directory
--genome                       Name of iGenomes reference

Output:
--outdir                       Path where the results to be saved [Default: './results']

Kit files:
--kit                          Kit used to prep samples [Default: 'agilent_v5']
--bait                         Absolute path to bait file
--target                       Absolute path to target file
--target_bed                   Absolute path to target bed file (snpEff compatible format)

Genome/Variation files:
--dbsnp                        Absolute path to dbsnp file
--thousandg                    Absolute path to 1000G file
--clinvar                      Absolute path to ClinVar file
--exac                         Absolute path to ExAC file
--gnomad                       Absolute path to GnomAD
--gfasta                       Absolute path to genome fasta file
--bwa_index                    Absolute path to bwa genome index

Other options:
--project                      Uppnex project to user for SLURM executor

For more detailed information regarding the parameters and usage refer to package
documentation at https:// github.com/apeltzer/QBIC-ExoSeq
"""

// Variables and defaults
params.help = false
params.reads = false
params.singleEnd = false
params.genome = false
params.run_id = false

//Clipping options
params.notrim = false
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0

// Kit options
params.kit = 'agilent_v5'
params.bait = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].bait ?: false : false
params.target = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].target ?: false : false
params.target_bed = params.kitFiles[ params.kit ] ? params.kitFiles[ params.kit ].target_bed ?: false : false
params.dbsnp = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].dbsnp ?: false : false
params.thousandg = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].thousandg ?: false : false
params.clinvar = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].clinvar ?: false : false
params.exac = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].exac ?: false : false
params.gnomad = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].gnomad ?: false : false
params.gfasta = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].gfasta ?: false : false
params.bwa_index = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].bwa_index ?: false : false


//Output configuration
params.outdir = './results'
params.saveAlignedIntermediates = false

//Configuration parameters
params.clusterOptions = false
params.project = false
params.cpus = 2


// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// Nextflow version check
nf_required_version = '0.25.0'
try {
    if( ! workflow.nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
        }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

// Check blocks for ceratin required parameters, to see they are given and exists
if (!params.reads || !params.genome){
    exit 1, "Parameters '--reads' and '--genome' are required to run the pipeline"
}
if (!params.kitFiles[ params.kit ] && ['bait', 'target'].count{ params[it] } != 2){
    exit 1, "Kit '${params.kit}' is not available in pre-defined config, so " +
            "provide all kit specific files with option '--bait' and '--target'"
}
if (!params.metaFiles[ params.genome ] && ['gfasta', 'bwa_index', 'dbsnp', 'thousandg', 'clinvar', 'exac', 'gnomad'].count{ params[it] } != 7){
    exit 1, "Genome '${params.genome}' is not available in pre-defined config, so you need to provide all genome specific " +
            "files with options '--gfasta', '--bwa_index', '--dbsnp', '--thousandg', '--clinvar', '--exac' and '--gnomad'"
}

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming }

/*
 * STEP 1 - trim with trim galore
 */

if(params.notrim){
    trimmed_reads = read_files_trimming
    trimgalore_results = []
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy'

        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file('*fq.gz') into trimmed_reads
        file '*trimming_report.txt' into trimgalore_results

        script:
        single = reads instanceof Path
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (params.singleEnd) {
            """
            trim_galore --gzip $c_r1 $tpc_r1 $reads
            """
        } else {
            """
            trim_galore --paired --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
            """
        }
    }
}


/*
 * STEP 2 - Map with BWA Mem
 */

process bwamem {
    tag "$name"

    input:
    set val(name), file(reads) from trimmed_reads

    output:
    set val(name), file("${name}_bwa.sam") into raw_aln_sam

    script:
    rg="\'@RG\\tID:${params.run_id}\\tSM:${params.run_id}\\tPL:illumina\'"
    if(params.singleEnd){
        """
         bwa mem \\
        -R $rg \\
        -t ${task.cpus} \\
        -k 2 \\
        $params.bwa_index \\
        $reads \\
            > ${name}_bwa.sam
    """
    } else {
    """
        bwa mem \\
        -R $rg \\
        -t ${task.cpus} \\
        -k 2 \\
        $params.bwa_index \\
        $reads\\
            > ${name}_bwa.sam
    """
}
}

/*
*  STEP 3 - Convert to BAM, sort BAM  
*/ 

process sortSam {
    tag "${name}"
    publishDir "${params.outdir}/BWAmem", mode: 'copy',
        saveAs: {filename -> params.saveAlignedIntermediates ? "aligned_sorted"/$filename : null }
    
    input:
    set val(name), file(raw_sam) from raw_aln_sam
    

    output: 
    set val(name), file("${raw_sam}.sorted.bam") into samples_sorted_bam

    script:
    def avail_mem = task.memory == null ? '' : "-m ${task.memory.toBytes() / task.cpus}"
    """
    samtools sort \\
        $raw_sam \\
        -@ ${task.cpus} $avail_mem \\
        -o ${raw_sam}.sorted.bam
    """
}

/*
*  STEP 4 - Mark PCR duplicates in sorted BAM file
*/ 

process markDuplicates {
    tag "${name}"
    publishDir "${params.outdir}/${name}/metrics", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".dup_metrics") > 0 ? filename : null }

    input:
    set val(name), file(sorted_bam) from samples_sorted_bam

    output:
    set val(name), file("${name}_markdup.bam") into samples_markdup_bam
    file("${name}.dup_metrics") into markdup_results

    script:
    """
        picard MarkDuplicates \\
        INPUT=$sorted_bam \\
        OUTPUT=${name}_markdup.bam \\
        METRICS_FILE=${name}.dup_metrics \\
        TMP_DIR=tmp \\
        VALIDATION_STRINGENCY=SILENT \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=false \\
        MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \\
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \\
        SORTING_COLLECTION_SIZE_RATIO=0.25 \\
        READ_NAME_REGEX=\"[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*\" \\
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \\
        VERBOSITY=INFO \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        CREATE_INDEX=false \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''
    """
}

// Recalibrate the bam file with known variants
process recal_bam_files {
    tag "${name}"

    input:
    set val(name), file(markdup_bam) from samples_markdup_bam

    output:
    set val(name), file("${name}_recal.bam"), file("${name}_recal.bai") into samples_recal_bam

    script:
    """
    gatk -T BaseRecalibrator \\
        -I $markdup_bam \\
        -R $params.gfasta \\
        -o ${name}_table.recal \\
        -cov ReadGroupCovariate \\
        -cov QualityScoreCovariate \\
        -cov CycleCovariate \\
        -cov ContextCovariate \\
        -U \\
        -OQ \\
        --default_platform illumina \\
        --knownSites $params.dbsnp \\
        -l INFO

    gatk -T PrintReads \\
        -BQSR ${name}_table.recal \\
        -I $markdup_bam \\
        -R $params.gfasta \\
        -o ${name}_recal.bam \\
        -baq RECALCULATE \\
        -U \\
        -OQ \\
        -l INFO
    """
}


// Realign the bam files based on known variants
process realign {
    tag "${name}"
    publishDir "${params.outdir}/${name}/alignment", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/realign/, "sorted_dupmarked_recalibrated_realigned")}

    input:
    set val(name), file(recal_bam), file(recal_bam_ind) from samples_recal_bam

    output:
    set val(name), file("${name}_realign.{bam,bai}") into (bam_vcall, bam_metrics)

    script:
    """
    gatk -T RealignerTargetCreator \\
        -I $recal_bam \\
        -R $params.gfasta \\
        -o ${name}_realign.intervals \\
        --known $params.dbsnp \\
        -l INFO

    gatk -T IndelRealigner \\
        -I $recal_bam \\
        -R $params.gfasta \\
        -targetIntervals ${name}_realign.intervals \\
        -o ${name}_realign.bam \\
        -l INFO
    """
}

// Calculate certain metrics
process calculateMetrics {
    tag "${name}"
    publishDir "${params.outdir}/${namesample}/metrics", mode: 'copy'

    input:
    set val(name), file(aligned_bam), file(aligned_bam_ind) from bam_metrics

    output:
    file("*{metrics,pdf}") into metric_files

    script:
    """
    picard CollectAlignmentSummaryMetrics \\
        INPUT=$aligned_bam \\
        OUTPUT=${name}.align_metrics \\
        REFERENCE_SEQUENCE=$params.gfasta \\
        VALIDATION_STRINGENCY=SILENT \\
        MAX_INSERT_SIZE=100000 \\
        ASSUME_SORTED=true \\
        IS_BISULFITE_SEQUENCED=false \\
        METRIC_ACCUMULATION_LEVEL="ALL_READS" \\
        STOP_AFTER=0 \\
        VERBOSITY=INFO \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''

    picard CollectInsertSizeMetrics \\
        HISTOGRAM_FILE=${name}_insert.pdf \\
        INPUT=$aligned_bam \\
        OUTPUT=${name}.insert_metrics \\
        VALIDATION_STRINGENCY=SILENT \\
        DEVIATIONS=10.0 \\
        MINIMUM_PCT=0.05 \\
        STOP_AFTER=0 \\
        METRIC_ACCUMULATION_LEVEL="ALL_READS" \\
        ASSUME_SORTED=true \\
        VERBOSITY=INFO \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''

    picard CalculateHsMetrics \\
        BAIT_INTERVALS=$params.bait \\
        TARGET_INTERVALS=$params.target \\
        INPUT=$aligned_bam \\
        OUTPUT=${name}.hs_metrics \\
        METRIC_ACCUMULATION_LEVEL="ALL_READS" \\
        VERBOSITY=INFO \\
        VALIDATION_STRINGENCY=SILENT \\
        QUIET=false \\
        COMPRESSION_LEVEL=5 \\
        MAX_RECORDS_IN_RAM=500000 \\
        CREATE_INDEX=false \\
        CREATE_MD5_FILE=false \\
        GA4GH_CLIENT_SECRETS=''
    """
}

// Call variants
process variantCall {
    tag "${name}"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/variants/, "raw_variants")}

    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_vcall

    output:
    set val(name), file("${name}_variants.vcf"), file("${name}_variants.vcf.idx") into raw_variants

    script:
    """
    gatk -T HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -o ${name}_variants.vcf \\
        --annotation HaplotypeScore \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --standard_min_confidence_threshold_for_calling 30.0 \\
        --dbsnp $params.dbsnp -l INFO
    """
}

// Select variants
process variantSelect {
    tag "${name}"

    input:
    set val(name), file(raw_vcf), file(raw_vcf_idx) from raw_variants

    output:
    set val(name), file("${name}_snp.vcf"), file("${name}_snp.vcf.idx") into raw_snp
    set val(name), file("${name}_indels.vcf"), file("${name}_indels.vcf.idx") into raw_indels

    script:
    """
    gatk -T SelectVariants \\
        -R $params.gfasta \\
        --variant $raw_vcf \\
        --out ${name}_snp.vcf \\
        --selectTypeToInclude SNP

    gatk -T SelectVariants \\
        -R $params.gfasta \\
        --variant $raw_vcf \\
        --out ${name}_indels.vcf \\
        --selectTypeToInclude INDEL \\
        --selectTypeToInclude MIXED \\
        --selectTypeToInclude MNP \\
        --selectTypeToInclude SYMBOLIC \\
        --selectTypeToInclude NO_VARIATION
    """
}

// Filter SNP
process filterSnp {
    tag "${name}"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set val(name), file(raw_snp), file(raw_snp_idx) from raw_snp

    output:
    set val(name), file("${sample}_filtered_snp.vcf"), file("${sample}_filtered_snp.vcf.idx") into filtered_snp

    script:
    """
    gatk -T VariantRecalibrator \\
        -R $params.gfasta \\
        --input $raw_snp \\
        --maxGaussians 4 \\
        --recal_file ${name}_snp.recal \\
        --tranches_file ${name}_snp.tranches \\
        -resource:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $params.hapmap \\
        -resource:omni,VCF,known=false,training=true,truth=false,prior=12.0 $params.omni \\
        -resource:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $params.dbsnp \\
        --mode SNP \\
        -an QD \\
        -an FS \\
        -an MQ

    gatk -T ApplyRecalibration \\
        -R $params.gfasta \\
        --out ${prefix}_filtered_snp.vcf \\
        --input $raw_snp \\
        --mode SNP \\
        --tranches_file ${prefix}_snp.tranches \\
        --recal_file ${prefix}_snp.recal
    """
}

// Filter indels
process filterIndel {
    tag "${name}"
    publishDir "${params.outdir}/${prefix}/variants", mode: 'copy'

    input:
    set val(name), file(raw_indel), file(raw_indel_idx) from raw_indels

    output:
    set val(name), file("${prefix}_filtered_indels.vcf"), file("${prefix}_filtered_indels.vcf.idx") into filtered_indels

    script:
    """
    gatk -T VariantFiltration \\
        -R $params.gfasta \\
        --variant $raw_indel \\
        --out ${name}_filtered_indels.vcf \\
        --filterName GATKStandardQD \\
        --filterExpression "QD < 2.0" \\
        --filterName GATKStandardReadPosRankSum \\
        --filterExpression "ReadPosRankSum < -20.0" \\
        --filterName GATKStandardFS \\
        --filterExpression "FS > 200.0"
    """
}

// Group filted snp and indels for each sample
filtered_snp
    .cross(filtered_indels)
    .map{ it -> [it[0][0], it[0][1], it[0][2], it[1][1], it[1][2]] }
    .set{ variants_filtered }

// Combine filtered snp and indels for each sample
process combineVariants {
    tag "$sample"

    input:
    set file(fsnp), file(fsnp_idx), file(findel), file(findel_idx) from variants_filtered

    output:
    set file("${sample}_combined_variants.vcf"), file("${sample}_combined_variants.vcf.idx") into (combined_variants_evaluate,combined_variants)

    script:
    """
    java -jar \$GATK_HOME/GenomeAnalysisTK.jar -T CombineVariants \\
        -R $params.gfasta \\
        --out ${sample}_combined_variants.vcf \\
        --genotypemergeoption PRIORITIZE \\
        --variant:${sample}_SNP_filtered $fsnp \\
        --variant:${sample}_indels_filtered $findel \\
        --rod_priority_list ${sample}_SNP_filtered,${sample}_indels_filtered
    """
}

// Evaluate variants
process variantEvaluate {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set file("${sample}_combined_variants.vcf"), file("${sample}_combined_variants.vcf.idx") from combined_variants_evaluate

    output:
    file "${sample}_combined_phased_variants.eval"

    script:
    """
    gatk -T VariantEval \\
        -R $params.gfasta \\
        --eval $phased_vcf \\
        --dbsnp $params.dbsnp \\
        -o ${sample}_combined_phased_variants.eval \\
        -L $params.target \\
        --doNotUseAllStandardModules \\
        --evalModule TiTvVariantEvaluator \\
        --evalModule CountVariants \\
        --evalModule CompOverlap \\
        --evalModule ValidationReport \\
        --stratificationModule Filter \\
        -l INFO
    """
}

// Annotate variants
process variantAnnotate {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/variants", mode: 'copy'

    input:
    set file(phased_vcf), file(phased_vcf_ind) from combined_variants

    output:
    file "*.{vcf,idx,snpeff}"

    script:
    """
    java -Xmx6g -jar \$SNPEFF_HOME/snpEff.jar \\
        -c \$SNPEFF_HOME/snpEff.config \\
        -i vcf \\
        -o gatk \\
        -o vcf \\
        -filterInterval $params.target_bed GRCh37.75 $phased_vcf \\
            > ${sample}_combined_phased_variants.snpeff

    gatk -T VariantAnnotator \\
        -R $params.gfasta \\
        -A SnpEff \\
        --variant $phased_vcf \\
        --snpEffFile ${sample}_combined_phased_variants.snpeff \\
        --out ${sample}_combined_phased_annotated_variants.vcf
    """
} */
