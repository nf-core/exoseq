# ![nf-core/ExoSeq](https://raw.githubusercontent.com/nf-core/Exoseq/master/docs/images/ExoSeq_logo.png)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.6-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/Lobby)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/exoseq.svg)](https://hub.docker.com/r/nfcore/exoseq/)
![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)
## Introduction

**nfcore/ExoSeq** is a bioinformatics analysis pipeline that performs best-practice analysis pipeline for Exome Sequencing data.

The pipeline is built based on [GATK](https://software.broadinstitute.org/gatk/best-practices/) best practices using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. The main steps done by pipeline are the following (more information about the processes can be found [here](docs/processes.md)).

* Alignment
* Marking Duplicates
* Recalibration
* Realignment
* Variant Calling
* Variant Filtration
* Variant Evaluation
* Variant Annotation

## Documentation
The nfcore/ExoSeq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

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

## Credits
The pipeline was initally developed by Senthilkumar Panneerselvam ([@senthil10](https://github.com/senthil10)) with a little help from Phil Ewels ([@ewels](https://github.com/ewels)) at the National Genomics Infrastructure, part of SciLifeLab in Stockholm and has been extended by Alex Peltzer ([@apeltzer](https://github.com/apeltzer)), Marie Gauder ([@mgauder](https://github.com/mgauder)) from QBIC Tuebingen/Germany as well as Marc Hoeppner ([@marchoeppner](https://github.com/marchoeppner)) from IKMB Kiel/Germany.

Many thanks also to others who have helped out along the way too, including [@pditommaso](https://github.com/pditommaso), [@colindaven](https://github.com/colindaven).