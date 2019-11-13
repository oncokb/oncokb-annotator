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

## UPDATE: API token required
#### OncoKB API
When you run `MafAnnotator.py`, `FusionAnnotator.py` and `CnaAnnotator.py`, you need a token before accessing the OncoKB data via its web API. Please go to [OncoKB Data Access](https://oncokb.org/dataAccess) to learn how to get a OncoKB API token and specify your token  in the command line after you have that.
```
python ${FILE_NAME.py} -i ${INPUT_FILE} -o ${OUTPUT_FILE} -b ${ONCOKB_API_TOKEN}
``` 

## Columns added in the annotation files
| Column          	| Possible Values                                                                                                                                                            	 	 	| Description                                                                                                                                                                                                                      	|
|-----------------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| mutation_effect 	| Gain-of-function, Likely Gain-of-function, Loss-of-function, Likely Loss-of-function, Switch-of-function, Likely Switch-of-function, Neutral, Likely Neutral, Inconclusive, Unknown 	| The biological effect of a mutation/alteration on the protein function that gives rise to changes in the biological properties of cells expressing the mutant/altered protein compared to cells expressing the wildtype protein. 	|
| oncogenic       	| Oncogenic, Likely Oncogenic, Likely Neutral, Inconclusive Unknown                                                                                                             	 	| In OncoKB, “oncogenic” is defined as “referring to the ability to induce or cause cancer” as described in the second edition of The Biology of Cancer by Robert Weinberg (2014).                                                 	|
| Level_*         	| LEVEL_1, LEVEL_2A, LEVEL_2B, LEVEL_3A, LEVEL_3B, LEVEL_4, LEVEL_R1, LEVEL_R2                                                                                                      	| The treatments available for a mutation/alteration giving a tumor type.                                                                                                                                                          	|
| Highest_level   	| LEVEL_1, LEVEL_2A, LEVEL_2B, LEVEL_3A, LEVEL_3B, LEVEL_4, LEVEL_R1, LEVEL_R2                                                                                                      	| The highest level across all available treatments giving a tumor type.                                                                                                                                                           	|
| citations       	| PMID, Abstract, Website Link                                                                                                                                                 	 	 	| All citations related to a mutation/alteration                                                                                                                                                                                   	|

## FAQs
- **How to get the subversion of the ensembl transcript ID?**  
  vcf2maf is designed to work with all Ensembl releases and reference genomes, that is the reason the subversion is not included. 

  Within MSK, we are using GRCh37 which corresponding to the Ensembl release 75. You can get the full list of versioned transcript IDs here http://grch37.ensembl.org/biomart/martview/40d921c6ab6956144cf6fb2e9a8ca093?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id_version&FILTERS=&VISIBLEPANEL=attributepanel.  
