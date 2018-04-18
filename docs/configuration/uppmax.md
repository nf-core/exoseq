# nfcore/ExoSeq: UPPMAX Configuration

The pipeline comes bundled with configurations to use the [Swedish UPPMAX](https://www.uppmax.uu.se/) clusters (tested on `milou`, `rackham`, `bianca` and `irma`). As such, you shouldn't need to add any custom configuration - everything _should_ work out of the box.

To use the pipeline on UPPMAX, you **must** specificy `-profile uppmax` when running the pipeline (as of Nov 2017).

Note that you will need to specify your UPPMAX project ID when running a pipeline. To do this, use the command line flag `--project <project_ID>`. The pipeline will exit with an error message if you try to run it pipeline with the UPPMAX config profile without a project.

**Optional Extra:** To avoid having to specify your project every time you run Nextflow, you can add it to your personal Nextflow config file instead. Add this line to `~/.nextflow/config`:

```groovy
params.project = 'project_ID' // eg. b2017123
```

## Running offline
If you are running the pipeline on Bianca or Irma, you will not have an active internet connection and some automated features will not be able to function. Specifically, you'll need to transfer the pipeline files and the singularity images manually.

First, to generate the singularity images, run the following command. Note that you need singularity installed - this is available on the other UPPMAX clusters (Milou and Rackham):

```bash
singularity pull -name docker://nfcore/exoseq
```

The nfcore/ExoSeq pipeline files can be downloaded from https://github.com/nfcore/ExoSeq

Once transferred, you may create a `work` directory for nextflow and copy the singularity containers in there:

```bash
mkdir -p work/singularity/
mv nfcore-exoseq-*.img work/singularity/
```

```bash
nextflow run /path/to/nfcore/ExoSeq -c /path/to/nfcore/ExoSeq/conf/base.config
```

(Note that you'll need the other common flags such as `--reads` and `--genome` in addition to this).

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
