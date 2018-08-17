#!/usr/bin/env nextflow
/*

========================================================================================
               nf-core/ E X O S E Q    B E S T    P R A C T I C E
========================================================================================
 #### Homepage / Documentation
 https://github.com/scilifelab/NGI-ExoSeq
 #### Authors
 Senthilkumar Panneerselvam @senthil10 <senthilkumar.panneerselvam@scilifelab.se>
 Phil Ewels @ewels <phil.ewels@scilifelab.se>
 Alex Peltzer @alex_peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
 Marie Gauder <marie.gauder@student.uni-tuebingen.de>

 Some code parts were lent from other NGI-Pipelines (e.g. CAW), specifically the error
 handling, logging messages. Thanks for that @CAW guys.
----------------------------------------------------------------------------------------
Developed based on GATK's best practise, takes set of FASTQ files and performs:
 - alignment (BWA)
 - recalibration (GATK)
 - realignment (GATK)
 - variant calling (GATK)
 - variant evaluation (SnpEff)
*/

// Help message
helpMessage = """
===============================================================================
nf-core/ExoSeq : Exome/Targeted sequence capture best practice analysis v${params.version}
===============================================================================

Usage: nextflow nf-core/ExoSeq --reads '*_R{1,2}.fastq.gz' --genome GRCh37 --kitfiles 'kitpath' --metafiles 'metapath'

This is a typical usage where the required parameters (with no defaults) were
given. The available paramaters are listed below based on category

Required parameters:
    --reads                       Absolute path to project directory
    --genome                      Name of iGenomes reference

Output:
    --outdir                      Path where the results to be saved [Default: './results']
    -w/--work-dir                 The temporary directory where intermediate data will be saved

Kit files:
    --kitfiles                    Path to kitfiles defined in metafiles.config
    --metafiles                   Path to metafiles defined in metafiles.config
    --kit                         Kit used to prep samples [Default: 'agilent_v5']

AWSBatch options:
    --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
    --awsregion                   The AWS Region for your AWS Batch job to run on

For more detailed information regarding the parameters and usage refer to package
documentation at https://github.com/nf-core/ExoSeq""".stripIndent()

// Variables and defaults
params.name = false
params.help = false
params.reads = false
params.singleEnd = false
params.run_id = false
params.aligner = 'bwa' //Default, but stay tuned for later ;-)
params.saveReference = true


// Output configuration
params.outdir = './results'
params.saveAlignedIntermediates = false
params.saveIntermediateVariants = false


// Clipping options
params.notrim = false
params.clip_r1 = 0
params.clip_r2 = 0
params.three_prime_clip_r1 = 0
params.three_prime_clip_r2 = 0



// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Kit & Genome options
params.kit = 'agilent_v5'
params.bait = params.genomes [ params.genome ].kits [params.kit].bait ?: false
params.target = params.genomes [ params.genome ].kits [params.kit].target ?: false
params.target_bed = params.genomes [ params.genome ].kits [params.kit].target_bed ?: false
params.dbsnp = params.genomes [ params.genome ].dbsnp ?: false
params.thousandg = params.genomes [ params.genome ].thousandg ?: false
params.mills = params.genomes [ params.genome ].mills ?: false
params.omni = params.genomes [ params.genome ].omni ?: false
params.gfasta = params.genomes [ params.genome ].gfasta ?: false
params.bwa_index = params.genomes [ params.genome ].bwa_index ?: false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// Check blocks for certain required parameters, to see they are given and exist
if ((!params.reads || !params.genome) && !workflow.profile == 'test'){
    exit 1, "Parameters '--reads' and '--genome' are required to run the pipeline"
}
if (!params.kitfiles){
    exit 1, "No Exome Capturing Kit specified!"
}
if (!params.metafiles){
    exit 1, "No Exome Metafiles specified!"
}
//AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

//Channelize gfasta into multiple channels for later
Channel.fromPath(params.gfasta)
       .into{ch_gfasta_for_bwa_mapping; ch_gfasta_for_recal; ch_gfasta_for_bqsr; ch_gfasta_for_variantcall; ch_gfasta_for_multimetrics; ch_gfasta_for_metrics} 

/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { read_files_fastqc; read_files_trimming }
}


// Validate Input indices for BWA Mem and GATK
if(params.aligner == 'bwa' ){
    bwaId = Channel
        .fromPath("${params.gfasta}.bwt")
        .ifEmpty { exit 1, "BWA index not found: ${params.gfasta}.bwt" }
}

// Create a summary for the logfile
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads ?: params.readPaths
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome Assembly']       = params.genome
summary['FastA']   = params.genomes [params.genome].gfasta ?: false
summary['DBSNP']   = params.genomes [params.genome].dbsnp ?: false
summary['Mills']   = params.genomes [params.genome].mills ?: false
summary['Omni']    = params.genomes [params.genome].omni ?: false
summary['1000G']   = params.genomes [params.genome].thousandg ?: false
summary['Exome Kit'] = params.kit
summary['- Bait'] = params.bait 
summary['- Target'] = params.target
summary['- Target Bed'] = params.target_bed
summary['Trim R1'] = params.clip_r1
summary['Trim R2'] = params.clip_r2
summary["Trim 3' R1"] = params.three_prime_clip_r1
summary["Trim 3' R2"] = params.three_prime_clip_r2
if(params.aligner == 'bwa'){
    summary['Aligner'] = "BWA Mem"
    if(params.bwa_index)          summary['BWA Index']   = params.bwa_index
    else if(params.gfasta)          summary['Fasta Ref']    = params.gfasta
}
summary['Save Intermediate Aligned Files'] = params.saveAlignedIntermediates ? 'Yes' : 'No'
summary['Save Intermediate Variant Files'] = params.saveIntermediateVariants ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="

try {
    if( ! workflow.nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
        }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


Channel.fromPath("${params.dbsnp}")
       .into{ch_dbsnp_for_baserecal; ch_dbsnp_for_multimetrics;ch_dbsnp_for_haplotypecaller; ch_dbsnp_for_vcf_index}

Channel.fromPath("${params.target}")
       .into{ch_kit_target_for_recal; ch_kit_target_for_bqsr; ch_kit_target_for_vcall; ch_kit_target_for_metrics}

Channel.fromPath("${params.bait}")
       .into{ch_kit_targetbait_for_multimetrics; ch_kit_targetbait_for_metrics}


// Build BWA Index if this is required
if(params.aligner == 'bwa' && !params.bwa_index){
    // Create Channels
    fasta_for_bwa_index = Channel
        .fromPath("${params.gfasta}")
    fasta_for_samtools_index = Channel
        .fromPath("${params.gfasta}")
    fasta_for_dict_index = Channel
        .fromPath("${params.gfasta}")
    // Create a BWA index for non-indexed genomes
    process makeBWAIndex {
        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_for_bwa_index

        output:
        file "*.{amb,ann,bwt,pac,sa,fasta}" into bwa_index

        script:
        """
        bwa index $fasta
        """
    }
    // Create a FastA index for non-indexed genomes
    process makeFastaIndex {
        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from fasta_for_samtools_index

        output:
        file "*.fai" into samtools_index, ch_gfasta_index_for_realign, ch_gfasta_index_for_bqsr, ch_gfasta_index_for_vcall, ch_gfasta_index_for_metrics

        script:
        """
        samtools faidx $fasta
        """
    }
    process makeSeqDict {
        tag "$params.gfasta"
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
        
        input:
        file fasta from fasta_for_dict_index

        output:
        file "*.dict" into ch_dict_index_for_recal, ch_dict_index_for_bqsr, ch_dict_index_for_vcall

        script:
        """
        gatk CreateSequenceDictionary --REFERENCE $fasta --OUTPUT "${fasta.baseName}.dict"
        """
    }

    process buildVCFIndex {
        publishDir path: { params.saveReference ? "${params.outdir}/reference_genome" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file(f_reference) from ch_dbsnp_for_vcf_index

        output:
        file("${f_reference}.idx") into ch_dbsnp_vcf_index

        script:
        """
        igvtools index ${f_reference}
        """
}

} else {
    bwa_index = file("${params.bwa_index}")
}

/*
 * 
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
    trimgalore_logs = []
} else {
    process trim_galore {
        tag "$name"
        publishDir "${params.outdir}/trim_galore", mode: 'copy', 
            saveAs: {filename -> 
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
                else params.saveTrimmed ? filename : null
            }
        input:
        set val(name), file(reads) from read_files_trimming

        output:
        set val(name), file(reads) into trimmed_reads
        file '*trimming_report.txt' into trimgalore_results, trimgalore_logs


        script:
        single = reads instanceof Path
        c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
        c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
        tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
        tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
        if (params.singleEnd) {
            """
            trim_galore --gzip $c_r1 $tpc_r1 $reads --fastqc
            """
        } else {
            """
            trim_galore --paired --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads --fastqc
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
    file(gfasta) from ch_gfasta_for_bwa_mapping

    output:
    set val(name), file("${name}_bwa.bam") into samples_sorted_bam
    file '.command.log' into bwa_stdout


    script:
    def avail_mem = task.memory ? "-m ${task.memory.toMega().intdiv(task.cpus)}M" : ''
    rg="\'@RG\\tID:${params.run_id}\\tSM:${params.run_id}\\tPL:illumina\'"

    """
    bwa mem \\
    -R $rg \\
    -t ${task.cpus} \\
    $gfasta \\
    $reads \\
    | samtools sort ${avail_mem} -O bam - > ${name}_bwa.bam
    """
}


/*
*  STEP 4 - Mark PCR duplicates in sorted BAM file
*/

process markDuplicates {
    tag "${name}"
    publishDir "${params.outdir}/Picard_Markduplicates/metrics", mode: 'copy',
        saveAs: { filename -> filename.indexOf(".dup_metrics") > 0 ? filename : null }

    input:
    set val(name), file(sorted_bam) from samples_sorted_bam

    output:
    set val(name), file("${name}_markdup.bam"), file("${name}_markdup.bai") into samples_markdup_bam, samples_for_applyBQSR
    file("${name}.dup_metrics") into markdup_results
    file '.command.log' into markDuplicates_stdout

    script:
    """
        mkdir `pwd`/tmp
        gatk MarkDuplicates \\
        --INPUT $sorted_bam \\
        --OUTPUT ${name}_markdup.bam \\
        --METRICS_FILE ${name}.dup_metrics \\
        --REMOVE_DUPLICATES false \\
        --CREATE_INDEX true \\
        --TMP_DIR=`pwd`/tmp \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}




/*
 * Step 5 - Recalibrate BAM file with known variants and BaseRecalibrator
 *
*/
process recal_bam_files {
    tag "${name}"
    publishDir "${params.outdir}/GATK_Recalibration", mode: 'copy'


    input:
    file gfasta from ch_gfasta_for_recal
    file gfasta_index from ch_gfasta_index_for_realign
    file gfasta_dict from ch_dict_index_for_recal
    file dbsnp from ch_dbsnp_for_baserecal
    file target from ch_kit_target_for_recal
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_markdup_bam

    output:
    set val(name), file("${name}_table.recal") into samples_recal_reports
    file '.command.log' into gatk_stdout
    file '.command.log' into gatk_base_recalibration_results

    script:
    """
    gatk BaseRecalibrator \\
        -R $gfasta \\
        -I $markdup_bam \\
        -O ${name}_table.recal \\
        -L $target \\
        --known-sites $dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}

process applyBQSR {
    tag "${name}"
    publishDir "${params.outdir}/GATK_ApplyBQSR", mode: 'copy'

    input:
    file gfasta from ch_gfasta_for_bqsr
    file gfasta_index from ch_gfasta_index_for_bqsr
    file gfasta_dict from ch_dict_index_for_bqsr
    file target from ch_kit_target_for_bqsr
    set val(name), file("${name}_table.recal") from samples_recal_reports
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_for_applyBQSR

    output:
    set val(name), file("${name}.bam"), file("${name}.bai") into bam_vcall, bam_for_multiple_metrics, bam_for_hs_metrics

    script:
    """
    gatk ApplyBQSR \\
        -R $gfasta \\
        -I $markdup_bam \\
        --bqsr-recal-file ${name}_table.recal \\
        -O ${name}.bam \\
        -L $target \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}

/*
 * Generate relevant statistics
 *
*/
//TODO MAYBE STILL USE QUALIMAP BUT UPDATE IT TO USE BED 3 too? 

process picard_multiple_metrics {
	publishDir "${params.outdir}/Picard/MultipleMetrics", mode: 'copy'
 
	input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_for_multiple_metrics
    file vcfidx from ch_dbsnp_vcf_index
    file gfasta from ch_gfasta_for_multimetrics
    file dbsnp from ch_dbsnp_for_multimetrics
    file bait from ch_kit_targetbait_for_multimetrics

	output:
	file("${prefix}*") into ch_collect_multiple_metrics_output

	script:       

	"""
		gatk CollectMultipleMetrics \
		--PROGRAM MeanQualityByCycle \
		--PROGRAM QualityScoreDistribution \
		--PROGRAM CollectAlignmentSummaryMetrics \
		--PROGRAM CollectInsertSizeMetrics\
       	--PROGRAM CollectSequencingArtifactMetrics \
        --PROGRAM CollectQualityYieldMetrics \
	    --PROGRAM CollectGcBiasMetrics \
		--PROGRAM CollectBaseDistributionByCycle \
		--INPUT $realign_bam \
		--DB_SNP $dbsnp \
        --REFERENCE_SEQUENCE $gfasta \
		--INTERVALS $bait \
		--ASSUME_SORTED true \
		--QUIET true \
		--OUTPUT ${name} \
		--TMP_DIR tmp \
        --java-options -Xmx${task.memory.toGiga()}g
	"""
}	

process picard_hc_metrics {
    publishDir "${params.outdir}/Picard/HcMetrics", mode: 'copy'

    input:
    file gfasta from ch_gfasta_for_metrics
    file gfasta_index from ch_gfasta_index_for_metrics
    file target from ch_kit_target_for_metrics
    file bait from ch_kit_targetbait_for_metrics
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_for_hs_metrics

    output:
    file(outfile) into ch_hybrid_capture_metrics

    script:
    outfile = "${name}" + "_"+ ".hybrid_selection_metrics.txt"

    """
        gatk CollectHsMetrics \
               --INPUT $realign_bam \
               --OUTPUT $outfile \
               --TARGET_INTERVALS $target \
               --BAIT_INTERVALS $bait \
               --REFERENCE_SEQUENCE $gfasta \
               --TMP_DIR tmp \
               --java-options -Xmx${task.memory.toGiga()}g
        """
}

/*
 * Step 8 - Call Variants with HaplotypeCaller in GVCF mode (differentiate between exome and whole genome data here)
 *
*/
process variantCall {
    tag "${name}"
    publishDir "${params.outdir}/GATK_VariantCalling/", mode: 'copy',
        saveAs: {filename -> filename.replaceFirst(/variants/, "raw_variants")}

    input:
    file gfasta from ch_gfasta_for_variantcall
    file gfasta_index from ch_gfasta_index_for_vcall
    file gfasta_dict from ch_dict_index_for_vcall
    file dbsnp from ch_dbsnp_for_haplotypecaller
    file target from ch_kit_target_for_vcall
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_vcall

    output:
    set val(name), file("${name}_variants.vcf"), file("${name}_variants.vcf.idx") into raw_variants

    script:
    """
    gatk HaplotypeCaller \\
        -I $realign_bam \\
        -R $gfasta \\
        -O ${name}_variants.vcf \\
        -ERC GVCF \\
        -L $target \\
        --create-output-variant-index \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp $dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}


/*
* Step 9 - Generate a YAML file for software versions in the pipeline
* This is then parsed by MultiQC and the report feature to produce a final report with the software Versions in the pipeline.
*/

process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo "$params.version" &> v_nfcore_exoseq.txt
    echo "$workflow.nextflow.version" &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    samtools --version &> v_samtools.txt
    bwa &> v_bwa.txt 2>&1 || true
    gatk BaseRecalibrator --version &> v_gatk.txt
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/**
* Step 10 - Generate MultiQC config file
*
*/

process GenerateMultiQCconfig {
  publishDir "${params.outdir}/MultiQC/", mode: 'copy'

  input:

  output:
  file("multiqc_config.yaml") into multiQCconfig

  script:
  """
  touch multiqc_config.yaml
  echo "custom_logo_title: 'Exome Analysis Workflow'" >> multiqc_config.yaml
  echo "extra_fn_clean_exts:" >> multiqc_config.yaml
  echo "- _R1" >> multiqc_config.yaml
  echo "- _R2" >> multiqc_config.yaml
  echo "report_header_info:" >> multiqc_config.yaml
  echo "- Exoseq version: ${params.version}" >> multiqc_config.yaml
  echo "- Command Line: ${workflow.commandLine}" >> multiqc_config.yaml
  echo "- Directory: ${workflow.launchDir}" >> multiqc_config.yaml
  echo "- Genome: "${params.gfasta} >> multiqc_config.yaml
  echo "  dbSNP : ${params.dbsnp}" >> multiqc_config.yaml
  echo "  Omni: ${params.omni}" >> multiqc_config.yaml
  echo "  Mills: ${params.mills}" >> multiqc_config.yaml
  echo "top_modules:" >> multiqc_config.yaml
  echo "- 'fastqc'" >> multiqc_config.yaml
  echo "- 'cutadapt'" >> multiqc_config.yaml
  echo "- 'bwa'" >> multiqc_config.yaml
  echo "- 'samtools'" >> multiqc_config.yaml
  echo "- 'gatk'" >> multiqc_config.yaml
  """
}


/*
* Step 12 - Collect metrics, stats and other resources with MultiQC in a single call
*/

process multiqc {
    tag "$name"
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiQCconfig
    file (fastqc:'fastqc/*') from fastqc_results.toList()
    file ('trimgalore/*') from trimgalore_results.toList()
    file ('gatk_base_recalibration/T*') from gatk_base_recalibration_results.toList()
    file ('gatk_picard_duplicates/*') from markdup_results.toList()
    file ('software_versions/*') from software_versions_yaml.toList()
    file ('*') from ch_collect_multiple_metrics_output.toList()
    file ('*') from ch_hybrid_capture_metrics.toList()


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
    multiqc -f $rtitle $rfilename --config $multiQCconfig .
    """
}


/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def exoMessage() {
  // Display ExoSeq message
  log.info "nf-core/ExoSeq ANALYSIS WORKFLOW ~ ${params.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outdir
  log.info "Genome      : " + params.gfasta
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def versionMessage() {
  // Display version message
  log.info "nf-core/ExoSeq ANALYSIS WORKFLOW"
  log.info "  version   : " + version
  log.info workflow.commitId ? "Git info    : ${workflow.repository} - ${workflow.revision} [${workflow.commitId}]" : "  revision  : " + this.grabRevision()
}


workflow.onComplete {
  // Display complete message
  this.nextflowMessage()
  this.exoMessage()
  this.minimalInformationMessage()
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError {
  // Display error message
  this.nextflowMessage()
  this.exoMessage()
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
