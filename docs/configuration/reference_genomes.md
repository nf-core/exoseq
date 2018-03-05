# NGI-ExoSeq: Reference Genomes Configuration

The NGI-ExoSeq pipeline needs a reference genome for alignment and annotation. If not already available, start by downloading the relevant reference, for example from the [GATK Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle).

The minimal requirements are a FASTA file. If a BWA reference is also specified, the pipeline won't have to generate them and will run faster.

Reference genome paths can be specified on the command line each time you run with `--bwa_index` and `--fasta`. Fasta is only required if building a BWA-MEM index.

## Adding paths to a config file
Specifying long paths every time you run the pipeline is a pain. To make this easier, the pipeline comes configured to understand reference genome keywords which correspond to preconfigured paths, meaning that you can just specify `--genome ID` when running the pipeline. This method is used on known systems such as UPPMAX.

Note that this genome key can also be specified in a config file if you always use the same genome.

To use this system, add paths to your config file using the following template:

```groovy
params {
  genomes {
    'YOUR-ID' {
      fasta  = '<PATH TO FASTA FILE>/genome.fa'
      bwa    = '<PATH TO BWA INDEX>/BWAIndex/'
    }
    'OTHER-GENOME' {
      // [..]
    }
  }
  // Optional - default genome. Ignored if --genome 'OTHER-GENOME' specified on command line
  genome = 'YOUR-ID'
}
```

You can add as many genomes as you like as long as they have unique IDs.

## GATK Resource Bundle
We are currently making use of the [GATK Resource Bundle](https://software.broadinstitute.org/gatk/download/bundle) human genome reference build `GRCh37`, but we would also be able to support `GRCh38`, if needed.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---