
## UPDATE: We recently updated our Level of Evidence, please see [Levels of Evidence section](#levels-of-evidence) for more information
## UPDATE: API token required, please see [OncoKB API section](#oncokb-api) for more information

# oncokb-annotator <a href="https://ascopubs.org/doi/full/10.1200/PO.17.00011"><img src="https://img.shields.io/badge/DOI-10.1200%2FPO.17.00011-1c75cd" /></a>  

Annotates variants in MAF with OncoKB annotation.

Please try 
* `python MafAnnotator.py -h`
* `python FusionAnnotator.py -h`
* `python CnaAnnotator.py -h`
* `python ClinicalDataAnnotator.py -h`
* `python OncoKBPlots.py -h`

Example input files are under [data](data). An example script is here: [example.sh](example.sh)

We recommend processing MAF files by [vcf2maf](https://github.com/mskcc/vcf2maf/) with [MSK override isoforms](https://github.com/mskcc/vcf2maf/blob/master/data/isoform_overrides_at_mskcc) before using the `MafAnnotator` here.

## Levels of Evidence
We recently updated our [Levels of Evidence](news...). 
The core changes are Level 2A has been changed to Level 2 and Level 2B has been merged to Level 3B.
If you want to annotate your files with previous Levels of Evidence, please checkout to branch levels-v1.

## OncoKB API
When you run `MafAnnotator.py`, `FusionAnnotator.py` and `CnaAnnotator.py`, you need a token before accessing the OncoKB data via its web API. Please visit [OncoKB Data Access Page](https://www.oncokb.org/dataAccess) for more information about how to register an account and get an OncoKB API token.  
With the token listed under [OncoKB Account Settings Page](https://www.oncokb.org/account/settings), you could use it in the following format.
```
python ${FILE_NAME.py} -i ${INPUT_FILE} -o ${OUTPUT_FILE} -b ${ONCOKB_API_TOKEN}
``` 

## Columns added in the annotation files
| Column          	| Possible Values                                                                                                                                                            	 	 	| Description                                                                                                                                                                                                                      	|
|-----------------	|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	|
| mutation_effect 	| Gain-of-function, Likely Gain-of-function, Loss-of-function, Likely Loss-of-function, Switch-of-function, Likely Switch-of-function, Neutral, Likely Neutral, Inconclusive, Unknown 	| The biological effect of a mutation/alteration on the protein function that gives rise to changes in the biological properties of cells expressing the mutant/altered protein compared to cells expressing the wildtype protein. 	|
| oncogenic       	| Oncogenic, Likely Oncogenic, Likely Neutral, Inconclusive Unknown                                                                                                             	 	| In OncoKB, “oncogenic” is defined as “referring to the ability to induce or cause cancer” as described in the second edition of The Biology of Cancer by Robert Weinberg (2014).                                                 	|
| Level_*         	| LEVEL_1, LEVEL_2, LEVEL_3A, LEVEL_3B, LEVEL_4, LEVEL_R1, LEVEL_R2                                                                                                      	| The treatments available for a mutation/alteration giving a tumor type.                                                                                                                                                          	|
| Highest_level   	| LEVEL_1, LEVEL_2, LEVEL_3A, LEVEL_3B, LEVEL_4, LEVEL_R1, LEVEL_R2                                                                                                      	| The highest level across all available treatments giving a tumor type.                                                                                                                                                           	|
| citations       	| PMID, Abstract, Website Link                                                                                                                                                 	 	 	| All citations related to a mutation/alteration                                                                                                                                                                                   	|

## FAQs
- **How to get the subversion of the ensembl transcript ID?**  
  vcf2maf is designed to work with all Ensembl releases and reference genomes, that is the reason the subversion is not included. 

  Within MSK, we are using GRCh37 which corresponding to the Ensembl release 75. You can get the full list of versioned transcript IDs here http://grch37.ensembl.org/biomart/martview/40d921c6ab6956144cf6fb2e9a8ca093?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id_version&FILTERS=&VISIBLEPANEL=attributepanel.  
