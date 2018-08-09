#!/bin/bash
echo "Getting an example exome capture kit..."
wget https://github.com/nf-core/test-datasets/raw/exoseq/extra/kits.tar.bz2
tar xvjf kits.tar.bz2
echo "Fetching important database files now ..."
mkdir -p b37 && cd b37
wget https://github.com/nf-core/test-datasets/raw/exoseq/reference/1000G_omni2.5.b37.small.vcf.gz
wget https://github.com/nf-core/test-datasets/raw/exoseq/reference/1000G_phase1.snps.highconf.b37.small.vcf.gz
wget https://github.com/nf-core/test-datasets/raw/exoseq/reference/Mills_and_1000G_gold_standard.indels.b37.small.vcf.gz
wget https://github.com/nf-core/test-datasets/raw/exoseq/reference/dbsnp_138.b37.small.vcf.gz
wget https://github.com/nf-core/test-datasets/raw/exoseq/reference/human_g1k_v37_decoy.small.fasta
gunzip *.gz