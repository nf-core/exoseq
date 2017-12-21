Howto install Exome-Kits
========================

The pipeline follows [best-practices](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS) defined by the Broad Institute of MIT and Harvard. You need to specify the correct exome capture kits that were used to sequence your samples and prepare corresponding `interval_list` and `dict` files for that purpose. 

First Step
----------

Prepare sequence dictionary for your reference genome.
```
picard CreateSequenceDictionary O=GRCh37.dict R=GRCh37.fa 
```

Second Step
-----------

Prepare the `interval_list` files required by GATK for the analysis. 

```
picard BedToIntervalList I=S03723314_Covered.bed O=S03723314_Covered.interval_list SD=GRCh37.dict 
picard BedToIntervalList I=S03723314_Regions.bed O=S03723314_Regions.interval_list SD=GRCh37.dict 
``` 

Then specify these in a config file [as a new kit](https://github.com/apeltzer/QBIC-ExoSeq/blob/master/conf/local.config) and you should be fine. Otherwise, you can also use the parameters to specify the required kit files for your exome capture kit.






