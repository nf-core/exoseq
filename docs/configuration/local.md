## ngi-ExoSeq: Local configuration

If running the pipeline in a local environment, we highly recommend using Singularity.

# Singularity image

Many HPC environments are not able to run Docker due to security issues. Singularity is a tool designed to run on such HPC systems which is very similar to Docker. 

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity images for you. Instead, you'll have to do this yourself manually first, transfer the image file to a folder `work/singularity` and then start the analysis. Nextflow will then find the pre-cached singularity containers and proceed without downloading them directly. 

You will need all these images locally:

```
singularity pull qbicsoftware/qbic-singularity-gatk
singularity pull qbicsoftware/qbic-singularity-snpeff
singularity pull qbicsoftware/qbic-singularity-qualimap2
singularity pull qbicsoftware/qbic-singularity-picard
singularity pull qbicsoftware/qbic-singularity-samtools
singularity pull qbicsoftware/qbic-singularity-bwa
singularity pull qbicsoftware/qbic-singularity-fastqc
singularity pull qbicsoftware/qbic-singularity-trimgalore
```

If you use one of the typical profiles, you'd probably not need to specify anything - the process configuration via the `base.config` profile will automatically fetch the singularity containers for runnning single steps once. 