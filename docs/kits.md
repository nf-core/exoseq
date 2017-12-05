Howto install kits
==================

First Step
----------

Prepare sequence dictionary for genome
```
picard CreateSequenceDictionary O=GRCh37.dict R=GRCh37.fa 
```

Second Step
-----------

Prepare kit BED files

```
picard BedToIntervalList I=S03723314_Covered.bed O=S03723314_Covered.interval_list SD=GRCh37.dict 
picard BedToIntervalList I=S03723314_Regions.bed O=S03723314_Regions.interval_list SD=GRCh37.dict 
``` 

Then specify these in your config file as a new kit and you should be fine. 






