# ![NGI-ExoSeq](https://raw.githubusercontent.com/SciLifeLab/NGI-ExoSeq/master/docs/images/NGI-ExoSeq_logo.png)

### Introduction

NGI-ExoSeq is a bioinformatics package that performs best-practice analysis pipeline for Exome Sequencing data at the [National Genomics Infastructure](https://ngisweden.scilifelab.se/)
at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

The pipeline is built based on [GATK](https://software.broadinstitute.org/gatk/best-practices/) best practices using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. The main steps done by pipeline are the following (more information about the processes can be found [here](docs/processes.md)).

* Alignment
* Marking Duplicates
* Recalibration
* Realignment
* Variant Calling
* Variant Filtration
* Variant Evaluation
* Variant Annotation

### Documentation
The NGI-ExoSeq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Pipeline installation and configuration instructions](docs/installation.md)
2. Pipeline configuration
   * [Local installation](docs/configuration/local.md)
   * [Amazon Web Services](docs/configuration/aws.md)
   * [Swedish UPPMAX clusters](docs/configuration/uppmax.md)
   * [Swedish cs3e Hebbe cluster](docs/configuration/c3se.md)
   * [Tübingen QBiC clusters](docs/configuration/qbic.md)
   * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
   * [Preparing custom exome capture kits](docs/kits.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
The pipeline was initally developed by Senthilkumar Panneerselvam ([@senthil10](https://github.com/senthil10)) with a little help from Phil Ewels ([@ewels](https://github.com/ewels)) and has been adapted by Alex Peltzer ([@apeltzer](https://github.com/apeltzer)) and Marie Gauder ([@mgauder](https://github.com/mgauder)) to our needs at QBiC, Tübingen, Germany.

In addition, we would like to recognise:
* Developers at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden 
* Developers at [National Genomics Infrastructure](https://github.com/orgs/NationalGenomicsInfrastructure/people) for their help, supports and suggestions
* [UPPMAX](http://www.uppmax.uu.se/) team
* [Nextflow](https://www.nextflow.io/docs/latest/index.html#) team
* [BWA](http://bio-bwa.sourceforge.net/) team
* [GATK](https://software.broadinstitute.org/gatk/) team
* [PICARD](http://broadinstitute.github.io/picard/) team
* [SnpEff](http://snpeff.sourceforge.net/) team


### Participating Institutes
NGI-ExoSeq is now used by a number of core sequencing and bioinformatics facilities. Some of these are listed below. If you use this pipeline too, please let us know in an issue and we will add you to the list.

<table>
  <tr>
    <td><img src="https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/NGI_logo.png" width="200"></td>
    <td>National Genomics Infrastructure (NGI), Sweden</td>
    <td>https://ngisweden.scilifelab.se/</td>
  </tr>
  <tr>
    <td><img src="https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/QBiC_logo.png" width="200"></td>
    <td>Quantitative Biology Center (QBiC), Germany</td>
    <td>https://portal.qbic.uni-tuebingen.de/portal/</td>
  </tr>
</table>

---

[![SciLifeLab](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](https://raw.githubusercontent.com/SciLifeLab/NGI-RNAseq/master/docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---