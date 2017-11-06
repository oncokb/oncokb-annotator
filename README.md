# oncokb-annotator
Annotates variants in MAF with OncoKB annotation.

Please try 
* `python MafAnnotator.py -h`
* `python FusionAnnotator.py -h`
* `python CnaAnnotator.py -h`
* `python ClinicalDataAnnotator.py -h`
* `python OncoKBPlots.py -h`

Example input files are under [data](data). An example script is here: [example.sh](example.sh)

We recommend processing MAF files by [vcf2maf](https://github.com/mskcc/vcf2maf/) with [MSK override isoforms](https://github.com/mskcc/vcf2maf/blob/master/data/isoform_overrides_at_mskcc) before using the `MafAnnotator` here.
