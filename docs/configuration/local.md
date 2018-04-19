## nf-core/ExoSeq: Local configuration

If running the pipeline in a local environment, we highly recommend using Singularity.

# Singularity image

Many HPC environments are not able to run Docker due to security issues. Singularity is a tool designed to run on such HPC systems which is very similar to Docker. 

If you intend to run the pipeline offline, nextflow will not be able to automatically download the singularity image for you. Instead, you'll have to do this yourself manually first, transfer the image file to a folder `work/singularity` and then start the analysis. Nextflow will then find the pre-cached singularity container and proceed without downloading them directly. 

You will only need this image locally:

```
singularity pull docker://nfcore/exoseq
```

If you use one of the typical profiles, you'd probably not need to specify anything - the process configuration via the `base.config` profile will automatically fetch the singularity container for runnning single steps once. 