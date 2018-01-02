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
version = '0.8'

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
params.exome = true
params.aligner = 'bwa' //Default, but stay tuned for later ;-) 
params.saveReference = true

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
params.mills = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].mills ?: false : false
params.omni = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].omni ?: false : false
params.gfasta = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].gfasta ?: false : false
params.bwa_index = params.metaFiles[ params.genome ] ? params.metaFiles[ params.genome ].bwa_index ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"


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

/**
Validate Input indices for BWA Mem and GATK
* 
*/ 
if(params.aligner == 'bwa' ){
    bwaId = Channel
        .fromPath("${params.gfasta}.bwt")
        .ifEmpty { exit 1, "BWA index not found: ${params.gfasta}.bwt" }
}

/**
* Build index for BWA if non exists
* 
*/


if(params.aligner == 'bwa' && !params.bwa_index){
    //Create Channels
    fasta_for_bwa_index = Channel
        .fromPath("${params.gfasta}")
    fasta_for_samtools_index = Channel
        .fromPath("${params.gfasta}")
    //Create a BWA index for non-indexed genomes
    process makeBWAIndex {
        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_for_bwa_index

        output:
        file "*.{amb,ann,bwt,pac,sa}" into bwa_index

        script:
        """
        bwa index $fasta
        """
    }
    //Create a FastA index for non-indexed genomes
    process makeFastaIndex {
        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
        
        input:
        file fasta from fasta_for_samtools_index

        output:
        file "*.fai" into samtools_index

        script:
        """
        samtools faidx $fasta
        """
    }
} else {
    bwa_index = Channel.fromPath("${params.bwa_index}")
}


/*
* STEP 0 - FastQC
*
*/

process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results
    file '.command.out' into fastqc_stdout

    script:
    """
    fastqc -q $reads
    """
}


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
        file '*trimming_report.txt' into trimgalore_results, trimgalore_logs

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
    file(bwa_index) from bwa_index

    output:
    set val(name), file("${name}_bwa.sam") into raw_aln_sam
    file '.command.log' into bwa_stdout


    script:
    rg="\'@RG\\tID:${params.run_id}\\tSM:${params.run_id}\\tPL:illumina\'"
    if(params.singleEnd){
    """
         bwa mem \\
        -R $rg \\
        -t ${task.cpus} \\
        -k 2 \\
        $params.gfasta \\
        $reads \\
            > ${name}_bwa.sam

    # Print version number to standard out
    echo "BWA Version:"\$(bwa 2>&1)
    """
    } else {
    """
        bwa mem \\
        -R $rg \\
        -t ${task.cpus} \\
        -k 2 \\
        $params.gfasta \\
        $reads\\
            > ${name}_bwa.sam

    # Print version number to standard out
    echo "BWA Version:"\$(bwa 2>&1)
    """

}
}

/*
*  STEP 3 - Convert to BAM, sort BAM  
*/ 

process sortSam {
    tag "${name}"
    publishDir "${params.outdir}/BWAmem", mode: 'copy',
        saveAs: {filename -> params.saveAlignedIntermediates ? "$filename" : null }
    
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
    file '.command.log' into markDuplicates_stdout

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

    # Print version number to standard out
    echo "File name: $bam_markduplicates Picard version "\$(picard  MarkDuplicates --version 2>&1)
    """
}

// Recalibrate the bam file with known variants
process recal_bam_files {
    tag "${name}"

    input:
    set val(name), file(markdup_bam) from samples_markdup_bam

    output:
    set val(name), file("${name}_recal.bam"), file("${name}_recal.bai") into samples_recal_bam
    file '.command.log' into gatk_stdout

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

    # Print version number to standard out
    echo "GATK version "\$(gatk --version 2>&1)
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
    set val(name), file("${name}_realign.bam"), file("${name}_realign.bai") into (bam_vcall, bam_metrics)

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

//Use QualiMap2, which is much nicer.
process qualiMap {
    tag "${name}"
    publishDir "${params.outdir}/${name}/qualimap", mode: 'copy'
    
    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_metrics

    output:
    file "${name}_qualimap" into qualimap_results
    file '.command.log' into qualimap_stdout

    script:
    gcref = ''
    if(params.genome == 'GRCh37') gcref = '-gd HUMAN'
    if(params.genome == 'GRCm38') gcref = '-gd MOUSE'
    """
    qualimap bamqc $gcref \\
    -bam $realign_bam \\
    -outdir ${name}_qualimap \\
    --skip-duplicated \\
    --collect-overlap-pairs \\
    -nt ${task.cpus} \\
    --java-mem-size=${task.memory.toGiga()}G \\

    # Print version number to standard out
    echo "QualiMap version "\$(qualimap --version 2>&1)
    """
}

// Call variants
process variantCall {
    tag "${name}"
    publishDir "${params.outdir}/${name}/variants", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/variants/, "raw_variants")}

    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_vcall

    output:
    set val(name), file("${name}_variants.vcf"), file("${name}_variants.vcf.idx") into raw_variants

    script:
    if(params.exome){
    """
    gatk -T HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -o ${name}_variants.vcf \\
        -ERC GVCF \\
        -L $params.target \\
        --annotation HaplotypeScore \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --standard_min_confidence_threshold_for_calling 30.0 \\
        --dbsnp $params.dbsnp -l INFO \\
        -variant_index_type LINEAR -variant_index_parameter 128000
    """
    } else { //We have a winner (genome)
    """
    gatk -T HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -o ${name}_variants.vcf \\
        -ERC GVCF \\
        --annotation HaplotypeScore \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --standard_min_confidence_threshold_for_calling 30.0 \\
        --dbsnp $params.dbsnp -l INFO \\
        -variant_index_type LINEAR -variant_index_parameter 128000
    """    
    }
}


//Combine Genotypes following Best Practices guidelines

process genotypegvcfs{
    tag "${name}"

    input:
    set val(name), file(raw_vcf), file(raw_vcf_idx) from raw_variants

    output:
    set val(name), file("${name}_gvcf.vcf"), file("${name}_gvcf.vcf.idx") into raw_gvcfs

    script:
    """
    gatk -T GenotypeGVCFs \\
    -R $params.gfasta \\
    --variant $raw_vcf \\
    -nt $task.cpus \\
    -o ${name}_gvcf.vcf \\
    """
}

// Select variants
process variantSelect {
    tag "${name}"

    input:
    set val(name), file(raw_vcf), file(raw_vcf_idx) from raw_gvcfs

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
process recalSNPs {
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
        -resource:omni,known=false,training=true,truth=true,prior=12.0 $params.omni \\
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 $params.thousandg \\
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $params.dbsnp \\
        --mode SNP \\
        -an QD \\
        -an FS \\
        -an MQ

    gatk -T ApplyRecalibration \\
        -R $params.gfasta \\
        --out ${name}_filtered_snp.vcf \\
        --input $raw_snp \\
        --mode SNP \\
        --tranches_file ${name}_snp.tranches \\
        --recal_file ${name}_snp.recal \\
        --ts_filter_level 99.5 \\
        -mode SNP 
    """
}

// Filter indels
process recalIndels {
    tag "${name}"
    publishDir "${params.outdir}/${name}/variants", mode: 'copy'

    input:
    set val(name), file(raw_indel), file(raw_indel_idx) from raw_indels

    output:
    set val(name), file("${name}_filtered_indels.vcf"), file("${name}_filtered_indels.vcf.idx") into filtered_indels

    script:
    """
    gatk -T VariantRecalibrator \\
        -R $params.gfasta \\
        --input $raw_indel \\
        --maxGaussians 4 \\
        --recal_file ${name}_indel.recal \\
        --tranches_file ${name}_indel.tranches \\
        -resource:mills,known=false,training=true,truth=true,prior=12.0 $params.mills \\
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $params.dbsnp \\
        -an QD -an DP -an FS -an SOR \\
        -mode INDEL 

    gatk -T ApplyRecalibration \\
        -R $params.gfasta \\
        --out ${name}_filtered_indels.vcf \\
        --input $raw_indel \\
        --mode SNP \\
        --tranches_file ${name}_indel.tranches \\
        --recal_file ${name}_indel.recal \\
        --ts_filter_level 99.0 \\
        -mode INDEL
    """
}

// Group SNPs and Indel again
filtered_snp
    .cross(filtered_indels)
    .map{ it -> [it[0][0], it[0][1], it[0][2], it[1][1], it[1][2]] }
    .set{ variants_filtered }

// Combine SNP + Indel Files into a single VCF again
process combineVariants {
    tag "$name"

    input:
    set file(fsnp), file(fsnp_idx), file(findel), file(findel_idx) from variants_filtered

    output:
    set file("${name}_combined_variants.vcf"), file("${name}_combined_variants.vcf.idx") into (combined_variants_evaluate,combined_variants_snpEff, combined_variants_gatk)

    script:
    """
    gatk -T CombineVariants \\
        -R $params.gfasta \\
        --out ${name}_combined_variants.vcf \\
        --genotypemergeoption PRIORITIZE \\
        --variant:${name}_SNP_filtered $fsnp \\
        --variant:${name}_indels_filtered $findel \\
        --rod_priority_list ${name}_SNP_filtered,${name}_indels_filtered
    """
}

// Annotate variants
process variantAnnotatesnpEff {
    tag "$name"
    publishDir "${params.outdir}/${name}/variants", mode: 'copy'

    input:
    set file(phased_vcf), file(phased_vcf_ind) from combined_variants_snpEff

    output:
    file "*.{snpeff}" into combined_variants_gatk_snpeff
    file '.command.log' into snpeff_stdout
    file 'SnpEffStats.csv' into snpeff_results

    script:
    """
        snpEff \\
        -c /usr/local/lib/snpEff/snpEff.config \\
        -i vcf \\
        -csvStats SnpEffStats.csv \\
        -o gatk \\
        -o vcf \\
        -filterInterval $params.target_bed GRCh37.75 $phased_vcf \\
            > ${name}_combined_phased_variants.snpeff 
        
        # Print version number to standard out
        echo "GATK version "\$(snpEff -version 2>&1)
    """
}

process variantAnnotateGATK{     
    tag "$name"
    publishDir "${params.outdir}/${name}/variants", mode: 'copy'

    input:
    set file(phased_vcf), file(phased_vcf_ind) from combined_variants_gatk
    file(phased_vcf_snpeff) from combined_variants_gatk_snpeff

    output:
    file "*.{vcd,idx}"

    script:
    """
    gatk -T VariantAnnotator \\
        -R $params.gfasta \\
        -A SnpEff \\
        --variant $phased_vcf \\
        --snpEffFile ${name}_combined_phased_variants.snpeff \\
        --out ${name}_combined_phased_annotated_variants.vcf
    """
}


// Evaluate variants
process variantEvaluate {
    tag "$name"
    publishDir "${params.outdir}/${name}/variants", mode: 'copy'

    input:
    set file("${name}_combined_variants.vcf"), file("${name}_combined_variants.vcf.idx") from combined_variants_evaluate

    output:
    file "${name}_combined_phased_variants.eval"
    file "${name}_combined_phased_variants.eval" into gatk_variant_eval_results

    script:
    """
    gatk -T VariantEval \\
        -R $params.gfasta \\
        --eval $phased_vcf \\
        --dbsnp $params.dbsnp \\
        -o ${name}_combined_phased_variants.eval \\
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

/*
 * Parse software version numbers
 */
software_versions = [
  'FastQC': null, 'Trim Galore!': null, 'BWA': null, 'Picard MarkDuplicates': null, 'GATK': null,
  'SNPEff': null, 'QualiMap': null, 'Nextflow': "v$workflow.nextflow.version"
]

/*
* Generates a YAML file for software versions in the pipeline
* This is then parsed by MultiQC and the report feature to produce a final report with the software Versions in the pipeline.
*
*/ 

process get_software_versions {
    cache false
    executor 'local'

    input:
    val fastqc from fastqc_stdout.collect()
    val trim_galore from trimgalore_logs.collect()
    val bwa from bwa_stdout.collect()
    val markDuplicates from markDuplicates_stdout.collect()
    val gatk from gatk_stdout.collect()
    val qualimap from qualimap_stdout.collect()
    val snpeff from snpeff_stdout.collect()

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    exec:
    software_versions['FastQC'] = fastqc[0].getText().find(/FastQC v(\S+)/) { match, version -> "v$version"}
    software_versions['Trim Galore!'] = trim_galore[0].getText().find(/Trim Galore version: (\S+)/) {match, version -> "v$version"}
    software_versions['BWA'] = bwa[0].getText().find(/Version: (\S+)/) {match, version -> "v$version"}
    software_versions['Picard MarkDuplicates'] = markDuplicates[0].getText().find(/Picard version ([\d\.]+)/) {match, version -> "v$version"}
    software_versions['GATK'] = gatk[0].getText().find(/GATK version ([\d\.]+)/) {match, version -> "v$version"} 
    software_versions['QualiMap'] = qualimap[0].getText().find(/QualiMap v.(\S+)/) {match, version -> "v$version"}
    software_versions['SNPEff'] = snpeff[0].getText().find(/SnpEff (\S+)/) {match, version -> "v$version" }

    def sw_yaml_file = task.workDir.resolve('software_versions_mqc.yaml')
    sw_yaml_file.text  = """
    id: 'ngi-exoseq'
    section_name: 'NGI-ExoSeq Software Versions'
    section_href: 'https://github.com/SciLifeLab/NGI-ExoSeq'
    plot_type: 'html'
    description: 'are collected at run time from the software output.'
    data: |
        <dl class=\"dl-horizontal\">
${software_versions.collect{ k,v -> "            <dt>$k</dt><dd>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</dd>" }.join("\n")}
        </dl>
    """.stripIndent()
}


/*
* Collects results from individual tools and produces a single output report for all the tools in the pipeline.
* 
*/

process multiqc {
    tag "$name"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('picard/*') from markdup_results.collect()
    file ('snpEff/*') from snpeff_results.collect()
    file ('gatk_variant_eval/*') from gatk_variant_eval_results.collect()
    file ('software_versions/*') from software_versions_yaml.collect()


    output:
    file '*multiqc_report.html' into multiqc_report
    file '*_data' into multiqc_data
    file '.command.err' into multiqc_stderr
    val prefix into multiqc_prefix

    script:
    prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config . 2>&1
    """

}

