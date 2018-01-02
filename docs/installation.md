Installation Instructions
=========================
## 1) Install NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v7+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

**You need NextFlow version >= 0.24 to run this pipeline.**

See [nextflow.io](https://www.nextflow.io/) and [NGI-NextflowDocs](https://github.com/SciLifeLab/NGI-NextflowDocs) for further instructions on how to install and configure Nextflow.

## 2) Install the Pipeline
This pipeline itself needs no installation - NextFlow will automatically fetch it from GitHub if `scilifelab/NGI-ExoSeq` is specified as the pipeline name.

### Offline use

If you need to run the pipeline on a system with no internet connection, you will need to download the files yourself from GitHub and run them directly:

```bash
wget https://github.com/scilifelab/NGI-ExoSeq/archive/master.zip
unzip master.zip -d /my-pipelines/
cd /my_data/
nextflow run /my-pipelines/NGI-ExoSeq-master
```

## 3) Bundle installation

GATK v3.X relies on bundle files coming that are recommended. We follow the guidelines here strictly, as incorporating other sources such as ExAC, GnomAD or similar can be used for annotation but not for calibrating variants according to guidelines at the Broad Institute.

Download the bundle data here to your local infrastructure:

https://software.broadinstitute.org/gatk/download/bundle

In case you are using Google Cloud, you may also access the storage bucket directly as listed on the webpage. Similar for Amazon AWS.

---

[![SciLifeLab](images/SciLifeLab_logo.png)](http://www.scilifelab.se/)
[![National Genomics Infrastructure](images/NGI_logo.png)](https://ngisweden.scilifelab.se/)

---
