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

1. [Installation and configuration](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)


### Credits
These scripts were written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/) at [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden.
The pipeline was developed by Senthilkumar Panneerselvam ([@senthil10](https://github.com/senthil10))
with a little help from Phil Ewels ([@ewels](https://github.com/ewels)).

In addition, we would like to recognise:
* Developers at [National Genomics Infrastructure](https://github.com/orgs/NationalGenomicsInfrastructure/people) for their help, supports and suggestions
* [UPPMAX](http://www.uppmax.uu.se/) team
* [Nextflow](https://www.nextflow.io/docs/latest/index.html#) team
* [BWA](http://bio-bwa.sourceforge.net/) team
* [GATK](https://software.broadinstitute.org/gatk/) team
* [PICARD](http://broadinstitute.github.io/picard/) team
* [SnpEff](http://snpeff.sourceforge.net/) team

---

[![SciLifeLab](https://raw.githubusercontent.com/SciLifeLab/NGI-ExoSeq/master/docs/images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](https://raw.githubusercontent.com/SciLifeLab/NGI-ExoSeq/master/docs/images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---