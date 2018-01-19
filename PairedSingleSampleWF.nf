#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-

========================================================================================
               NGI - E X O S E Q    B E S T    P R A C T I C E
========================================================================================
 #### Homepage / Documentation
 https://github.com/scilifelab/NGI-ExoSeq
 #### Authors
 Senthilkumar Panneerselvam @senthil10 <senthilkumar.panneerselvam@scilifelab.se>
 Phil Ewels @ewels <phil.ewels@scilifelab.se>
 Alex Peltzer @alex_peltzer <alexander.peltzer@qbic.uni-tuebingen.de>
 Marie Gauder <marie.gauder@student.uni-tuebingen.de>
----------------------------------------------------------------------------------------
Developed based on GATK's best practise, takes set of FASTQ files and performs:
 - alignment (BWA)
 - recalibration (GATK)
 - realignment (GATK)
 - variant calling (GATK)
 - variant evaluation (SnpEff)
*/

// Package version
version = '0.8.1'

// Help message
helpMessage = """
===============================================================================
NGI-ExoSeq : Exome/Targeted sequence capture best practice analysis v${version}
===============================================================================

Usage: nextflow SciLifeLab/NGI-ExoSeq --reads '*_R{1,2}.fastq.gz' --genome GRCh37

This is a typical usage where the required parameters (with no defaults) were
given. The available paramaters are listed below based on category

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
--mills                        Absolute path to Mills file
--omni                         Absolute path to Omni file
--gfasta                       Absolute path to genome fasta file
--bwa_index                    Absolute path to bwa genome index

Other options:
--exome                        Exome data, if this is not set, run as genome data
--project                      Uppnex project to user for SLURM executor

For more detailed information regarding the parameters and usage refer to package
documentation at https:// github.com/scilifelab/NGI-ExoSeq
"""

// Variables and defaults
params.name = false
params.help = false
params.reads = false
params.singleEnd = false
params.genome = false
params.run_id = false
params.exome = false //default genome, set to true to run restricting to exome positions
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

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

// Check blocks for certain required parameters, to see they are given and exist
if (!params.reads || !params.genome){
    exit 1, "Parameters '--reads' and '--genome' are required to run the pipeline"
}
if (!params.kitFiles[ params.kit ] && ['bait', 'target'].count{ params[it] } != 2){
    exit 1, "Kit '${params.kit}' is not available in pre-defined config, so " +
            "provide all kit specific files with option '--bait' and '--target'"
}
 if (!params.metaFiles[ params.genome ] && ['gfasta', 'bwa_index', 'dbsnp', 'thousandg', 'mills', 'omni'].count{ params[it] } != 6){
     exit 1, "Genome '${params.genome}' is not available in pre-defined config, so you need to provide all genome specific " +
             "files with options '--gfasta', '--bwa_index', '--dbsnp', '--thousandg', '--mills' and '--omni'"
 }

// Create a channel for input files

Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc; read_files_trimming }


// Validate Input indices for BWA Mem and GATK

if(params.aligner == 'bwa' ){
    bwaId = Channel
        .fromPath("${params.gfasta}.bwt")
        .ifEmpty { exit 1, "BWA index not found: ${params.gfasta}.bwt" }
}


// Create a summary for the logfile
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']       = params.genome
summary['WES/WGS']      = params.exome ? 'WES' : 'WGS'
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

// Build BWA Index if this is required

if(params.aligner == 'bwa' && !params.bwa_index){
    // Create Channels
    fasta_for_bwa_index = Channel
        .fromPath("${params.gfasta}")
    fasta_for_samtools_index = Channel
        .fromPath("${params.gfasta}")
    // Create a BWA index for non-indexed genomes
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
    // Create a FastA index for non-indexed genomes
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
    bwa_index = file("${params.bwa_index}")
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
    fastqc --version
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
        publishDir "${params.outdir}/trim_galore", mode: 'copy'

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
    def avail_mem = task.memory ? "-m ${task.memory.toBytes().intdiv(task.cpus)}" : ''
    """
    samtools sort \\
        $raw_sam \\
        -@ ${task.cpus}\\
        $avail_mem \\
        -o ${raw_sam}.sorted.bam
    """
}

/*
*  STEP 4 - Mark PCR duplicates in sorted BAM file
*/ 

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
        gatk MarkDuplicates \\
        --INPUT $sorted_bam \\
        --OUTPUT ${name}_markdup.bam \\
        --METRICS_FILE ${name}.dup_metrics \\
        --REMOVE_DUPLICATES false \\
        --CREATE_INDEX true \\
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
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_markdup_bam

    output:
    set val(name), file("${name}_table.recal") into samples_recal_reports
    file '.command.log' into gatk_stdout
    file '.command.log' into gatk_base_recalibration_results

    script:
    if(params.exome){
    """
    gatk BaseRecalibrator \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        -O ${name}_table.recal \\
        -L $params.target \\
        --known-sites $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    # Print version number to standard out
    echo "GATK version "\$(gatk BaseRecalibrator --version 2>&1)
    """
    } else {
    """
    gatk BaseRecalibrator \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        -O ${name}_table.recal \\
        --known-sites $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    # Print version number to standard out
    echo "GATK version "\$(gatk BaseRecalibrator --version 2>&1)
    """    
    }
}

process applyBQSR {
    tag "${name}"
    publishDir "${params.outdir}/GATK_ApplyBQSR", mode: 'copy'

    input:
    set val(name), file("${name}_table.recal") from samples_recal_reports
    set val(name), file(markdup_bam), file(markdup_bam_ind) from samples_for_applyBQSR

    output:
    set val(name), file("${name}.bam"), file("${name}.bai") into bam_vcall, bam_metrics

    script:
    if(params.exome){
    """
    gatk ApplyBQSR \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        --bqsr-recal-file ${name}_table.recal \\
        -O ${name}.bam \\
        -L $params.target \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    } else {
    """
    gatk ApplyBQSR \\
        -R $params.gfasta \\
        -I $markdup_bam \\
        --bqsr-recal-file ${name}_table.recal \\
        -O ${name}.bam \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """    
    }

}

/*
 * Step 7 - Determine quality metrics of mapped BAM files using QualiMap 2
 * 
*/ 
process qualiMap {
    tag "${name}"
    publishDir "${params.outdir}/Qualimap", mode: 'copy'
    
    input:
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_metrics

    output:
    file "${name}" into qualimap_results
    file '.command.log' into qualimap_stdout

    script:
    gcref = ''
    if(params.genome == 'GRCh37') gcref = '-gd HUMAN'
    if(params.genome == 'GRCm38') gcref = '-gd MOUSE'
    """
    qualimap bamqc $gcref \\
    -bam $realign_bam \\
    -outdir ${name} \\
    --skip-duplicated \\
    --collect-overlap-pairs \\
    -nt ${task.cpus} \\
    --java-mem-size=${task.memory.toGiga()}G \\
    # Print version number to standard out
    echo "QualiMap version "\$(qualimap --version 2>&1)
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
    set val(name), file(realign_bam), file(realign_bam_ind) from bam_vcall

    output:
    set val(name), file("${name}_variants.vcf"), file("${name}_variants.vcf.idx") into raw_variants

    script:
    if(params.exome){
    """
    gatk HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -O ${name}_variants.vcf \\
        -ERC GVCF \\
        -L $params.target \\
        --create-output-variant-index \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
    } else { //We have a winner (genome)
    """
    gatk HaplotypeCaller \\
        -I $realign_bam \\
        -R $params.gfasta \\
        -O ${name}_variants.vcf \\
        -ERC GVCF \\
        --create-output-variant-index \\
        --annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """    
    }
}

/*
 * Step 9 - Generate Software Versions Map
 * 
*/
software_versions = [
  'FastQC': null, 'Trim Galore!': null, 'BWA': null, 'GATK': null,
  'QualiMap': null, 'Nextflow': "v$workflow.nextflow.version"
]

/*
* Step 10 - Generate a YAML file for software versions in the pipeline
* This is then parsed by MultiQC and the report feature to produce a final report with the software Versions in the pipeline.
*/ 

process get_software_versions {
    cache false
    executor 'local'

    input:
    val fastqc from fastqc_stdout.collect()
    val trim_galore from trimgalore_logs.collect()
    val bwa from bwa_stdout.collect()
    val gatk from gatk_stdout.collect()
    val qualimap from qualimap_stdout.collect()

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    exec:
    software_versions['FastQC'] = fastqc[0].getText().find(/FastQC v(\S+)/) { match, version -> "v$version"}
    software_versions['Trim Galore!'] = trim_galore[0].getText().find(/Trim Galore version: (\S+)/) {match, version -> "v$version"}
    software_versions['BWA'] = bwa[0].getText().find(/Version: (\S+)/) {match, version -> "v$version"}
    software_versions['GATK'] = gatk[0].getText().find(/Version:([\d\.]+)/) {match, version -> "v$version"} 
    software_versions['QualiMap'] = qualimap[0].getText().find(/QualiMap v.(\S+)/) {match, version -> "v$version"}

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

/**
* Step 11 - Generate MultiQC config file
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
  echo "- Exoseq version: ${version}" >> multiqc_config.yaml
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
  echo "- 'qualimap'" >> multiqc_config.yaml
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
    file (fastqc:'fastqc/*') from fastqc_results.collect()
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('gatk_base_recalibration/T*') from gatk_base_recalibration_results.collect()
    file ('gatk_picard_duplicates/*') from markdup_results.collect()
    file ('qualimap/*') from qualimap_results.collect()
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
    multiqc -f $rtitle $rfilename --config $multiQCconfig . 
    """
}

