# ![QBiC](https://raw.githubusercontent.com/qbicsoftware/QBIC-ExoSeq/master/docs/images/qbic_logo.png) QBiC-ExoSeq

### Introduction

QBiC-ExoSeq is a bioinformatics package that performs best-practice analysis pipeline for Exome Sequencing data at the Quantitative Biology Center 
([QBiC](http://www.uni-tuebingen.de/en/facilities/zentrale-einrichtungen/quantitative-biology-center-qbic.html)), Tübingen, Germany. It is based upon the pipeline NGI-ExoSeq, which was developed at the 
[National Genomics Infastructure](https://ngisweden.scilifelab.se/) at [SciLifeLab Stockholm](https://www.scilifelab.se/platforms/ngi/), Sweden.

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
The QBiC-ExoSeq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Pipeline installation and configuration instructions](docs/installation.md)
2. Pipeline configuration
   * [Local installation](docs/configuration/local.md)
   * [Amazon Web Services](docs/configuration/aws.md)
   * [Swedish UPPMAX clusters](docs/configuration/uppmax.md)
   * [Swedish cs3e Hebbe cluster](docs/configuration/c3se.md)
   * [Tübingen QBIC clusters](docs/configuration/qbic.md)
   * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
   * [Preparing custom exome capture kits](docs/kits.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
The pipeline was initally developed by Senthilkumar Panneerselvam ([@senthil10](https://github.com/senthil10)) with a little help from Phil Ewels ([@ewels](https://github.com/ewels)) and has been adapted by Alex Peltzer ([@apeltzer](https://github.com/apeltzer) to our needs at QBiC, Tübingen, Germany.

In addition, we would like to recognise:
* Developers at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden 
* Developers at [National Genomics Infrastructure](https://github.com/orgs/NationalGenomicsInfrastructure/people) for their help, supports and suggestions
* [UPPMAX](http://www.uppmax.uu.se/) team
* [Nextflow](https://www.nextflow.io/docs/latest/index.html#) team
* [BWA](http://bio-bwa.sourceforge.net/) team
* [GATK](https://software.broadinstitute.org/gatk/) team
* [PICARD](http://broadinstitute.github.io/picard/) team
* [SnpEff](http://snpeff.sourceforge.net/) team

---

[![QBiC](https://raw.githubusercontent.com/apeltzer/QBIC-ExoSeq/master/docs/images/qbic_logo.png)](http://www.uni-tuebingen.de/en/facilities/zentrale-einrichtungen/quantitative-biology-center-qbic.html)
---
