
## UPDATE: API token required, please see [OncoKB™ API section](#oncokb-api) for more information
## UPDATE: v3.3 introduces additional annotation columns, please check your pipeline before using

# oncokb-annotator <a href="https://ascopubs.org/doi/full/10.1200/PO.17.00011"><img src="https://img.shields.io/badge/DOI-10.1200%2FPO.17.00011-1c75cd" /></a>

## Status

[![Run all python tests](https://github.com/oncokb/oncokb-annotator/workflows/Run%20all%20python%20tests/badge.svg)](https://github.com/oncokb/oncokb-annotator/actions?query=workflow%3A%22Run+all+python+tests%22) [![Compare Annotation](https://github.com/oncokb/oncokb-annotator/workflows/Compare%20Annotation/badge.svg)](https://github.com/oncokb/oncokb-annotator/actions?query=workflow%3A%22Compare+Annotation%22)

## Install dependencies
For python 3
```
pip install -r requirements/common.txt -r requirements/pip3.txt
```

For python 2.7
```
pip install -r requirements/common.txt -r requirements/pip2.7.txt
```


## Usage
Example input files are under [data](data). An example script is here: [example.sh](example.sh)

### MAF
Annotates variants in MAF(https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) with OncoKB™ annotation. Supports both python2 and python3.  
Get more details on the command line using `python MafAnnotator.py -h`.  

We recommend processing VCF files by [vcf2maf](https://github.com/mskcc/vcf2maf/) with [OncoKB™ isoforms](https://www.oncokb.org/api/v1/utils/allCuratedGenes) before using the `MafAnnotator` here.

#### Atypical Alteration
You can still use MAF format to annotate atypical alterations, such as MSI-H, TMB-H, EGFR vIII. Please see more examples [HERE](data/example_atypical_alterations.txt).  

### Copy Number Alteration
#### Use GISTIC data format
We use GISTIC 2.0 format by default. For more information, please see https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#discrete-copy-number-data, please see examples [HERE](data/example_cna.txt).   
Columns `Locus ID` and `Cytoband` are not required.
#### Individual CNA
You can also list copy number alteration individually by specifying `-f individual`, please see examples [HERE](data/example_individual_cna.txt).

Get more details on the command line using `python CnaAnnotator.py -h`.  

### Fusion
OncoKB™ offers to annotate functional fusions.
The fusion format for intragenic deletion is `GENE-intragenic` or `GENE-GENE`.
For other fusions, please use `GENEA-GENEB` or `GENEA-GENEB Fusion`.  

Get more details on the command line using `python FusionAnnotator.py -h`.  

### Structural Variant
OncoKB™ offers to annotate structural variant.
The types supported are DELETION, TRANSLOCATION, DUPLICATION, INSERTION, INVERSION, FUSION, UNKNOWN.
All other types will be converted to UNKNOWN.

All structural variants with two different gene partners, they will be considered as functional fusions.

Get more details on the command line using `python StructuralVariantAnnotator.py -h`.

### Clinical Data (Combine MAF+CNA+Fusion)
You can combine all annotation on sample/patient level using the clinical data annotator.  

Get more details on the command line using `python ClinicalDataAnnotator.py -h`.  

### Annotate with HGVSp_Short, HGVSp, HGVSg or Genomic Change
OncoKB™ MafAnnotator supports annotating the alteration with HGVSp, HGVSp_Short, HGVSg or Genomic Change format. Please specify the query type with -q parameter.
The acceptable values are HGVSp_Short, HGVSp, HGVSg and Genomic_Change(case-insensitive). Please see data/example.sh for examples.  
If you do not specify query type, the MafAnnotator will try to figure out the query type based on the headers.  

#### For HGVSp_Short
The annotator takes alteration from the column HGVSp_Short or Alteration  

#### For HGVSp
The annotator takes alteration from the column HGVSp or Alteration  

#### For HGVSg
The annotator takes alteration from the column HGVSg or Alteration  

#### For Genomic_Change
The annotator takes genomic change from columns Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1(Optional) and Tumor_Seq_Allele2.  
Typically Tumor_Seq_Allele1 is the reference allele, Tumor_Seq_Allele2 is the variant allele. This is why Tumor_Seq_Allele1 is optional.  
The annotator uses both if the value is different from Reference_Allele. Tumor_Seq_Allele2 has higher priority than Tumor_Seq_Allele1.  

Annotation with Genomic_Change is relatively slow. We need to annotate the variant first with GenomeNexus(https://www.genomenexus.org/) then get annotation one by one. There is a plan to improve this method. If you are annotating a lot of data, please prioritize using other query type if applicable. 


### Annotate with different reference genomes (GRCh37, GRCh38)
OncoKB™ MafAnnotator supports annotating the alteration with reference genome GRCh37 and GRCh38.  

The annotator will get the reference genome from MAF file column NCBI_Build or Reference_Genome.  
If there is no reference genome specified in the file, we will use the default reference genome through -r parameter.  

You can specify the default reference genome using -r parameter (This is only applicable to MafAnnotator.py).  
The acceptable values are GRCh37, GRCh38 (case in-sensitive).  

If both values are not specified, the annotator will use OncoKB™ default reference genome which is GRCh37.


## Levels of Evidence
Introducing [Simplified OncoKB™ Levels of Evidence](https://www.oncokb.org/levels):
- New Level 2, defined as “Standard care biomarker recommended by the NCCN or other expert panels predictive of response to an FDA-approved drug in this indication” (formerly Level 2A).
- Unified Level 3B, defined as “Standard care or investigational biomarker predictive of response to an FDA-approved or investigational drug in another indication” (combination of previous Levels 2B and 3B).

We have implemented these changes for 2 reasons:
- To be consistent with the [Joint Consensus Recommendation by AMP, ASCO and CAP](https://www.sciencedirect.com/science/article/pii/S1525157816302239?via%3Dihub) and the [ESMO Scale for Clinical Actionability of molecular Targets (ESCAT)](https://academic.oup.com/annonc/article/29/9/1895/5076792?searchresult=1)
- To reflect the clinical data that demonstrates patients with investigational predictive biomarkers for a specific tumor type based on compelling clinical evidence (currently Level 3A) are more likely to experience clinical benefit compared to patients with predictive biomarkers that are considered standard care in a different tumor type (previously Level 2B, now combined into Level 3B).


## OncoKB™ API
When you run `MafAnnotator.py`, `FusionAnnotator.py` and `CnaAnnotator.py`, you need a token before accessing the OncoKB™ data via its web API. Please visit [OncoKB™ Data Access Page](https://www.oncokb.org/dataAccess) for more information about how to register an account and get an OncoKB™ API token.  
With the token listed under [OncoKB™ Account Settings Page](https://www.oncokb.org/account/settings), you could use it in the following format.
```
python ${FILE_NAME.py} -i ${INPUT_FILE} -o ${OUTPUT_FILE} -b ${ONCOKB_API_TOKEN}
``` 


## Columns added in the annotation files using MafAnnotator/CnaAnnotator/StructuralVariantAnnotator/FusionAnnotator
| Column                    | Possible Values                                                                                                                                                                     | Description                                                                                                                                                                                                                      |
|---------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| ANNOTATED                 | True, False                                                                                                                                                                         | Whether the variant is annotated by OncoKB successfully.                                                                                                                                                                         |
| GENE_IN_ONCOKB            | True, False                                                                                                                                                                         | Whether the gene has been curated by the OncoKB Team.                                                                                                                                                                            |
| VARIANT_IN_ONCOKB         | True, False                                                                                                                                                                         | Whether the variant has been curated by the OncoKB Team. Note: when a variant does not exist, it may still have annotations.                                                                                                     |
| MUTATION_EFFECT           | Gain-of-function, Likely Gain-of-function, Loss-of-function, Likely Loss-of-function, Switch-of-function, Likely Switch-of-function, Neutral, Likely Neutral, Inconclusive, Unknown | The biological effect of a mutation/alteration on the protein function that gives rise to changes in the biological properties of cells expressing the mutant/altered protein compared to cells expressing the wildtype protein. |
| MUTATION_EFFECT_CITATIONS | PMID, Abstract, Website link                                                                                                                                                        | All citations related to the biological effect.                                                                                                                                                                                  |
| ONCOGENIC                 | Oncogenic, Likely Oncogenic, Likely Neutral, Inconclusive, Unknown, Resistance                                                                                                      | In OncoKB™, “oncogenic” is defined as “referring to the ability to induce or cause cancer” as described in the second edition of The Biology of Cancer by Robert Weinberg (2014).                                                 |
| LEVEL_*                   | Therapeutic implications                                                                                                                                                            | The leveled therapeutic implications.                                                                                                                                                                                            |
| HIGHEST_LEVEL             | LEVEL_1, LEVEL_2, LEVEL_3A, LEVEL_3B, LEVEL_4, LEVEL_R1, LEVEL_R2                                                                                                                   | The highest level of evidence for therapeutic implications. Order: LEVEL_R1 > LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4 > LEVEL_R2                                                                                       |
| HIGHEST_SENSITIVE_LEVEL   | LEVEL_1, LEVEL_2, LEVEL_3A, LEVEL_3B, LEVEL_4                                                                                                                                       | The highest sensitive level of evidence for therapeutic implications. Order: LEVEL_1 > LEVEL_2 > LEVEL_3A > LEVEL_3B > LEVEL_4                                                                                                   |
| HIGHEST_RESISTANCE_LEVEL  | LEVEL_R1, LEVEL_R2                                                                                                                                                                  | The highest resistance level of evidence for therapeutic implications. Order: LEVEL_R1 > LEVEL_R2                                                                                                                                |
| TX_CITATIONS              | PMID, Abstract, Website link                                                                                                                                                        | All citations related to therapeutic implications.                                                                                                                                                                               |
| LEVEL_Dx*                 | Tumor type the level of evidence is assigned to                                                                                                                                     | The leveled diagnostic implications.                                                                                                                                                                                             |
| HIGHEST_DX_LEVEL          | LEVEL_Dx1, LEVEL_Dx2, LEVEL_Dx3                                                                                                                                                     | The highest level of evidence for diagnostic implications.                                                                                                                                                                       |
| DX_CITATIONS              | PMID, Abstract, Website link                                                                                                                                                        | All citations related to diagnostic implications.                                                                                                                                                                                |
| LEVEL_Px*                 | Tumor type the level of evidence is assigned to                                                                                                                                     | The leveled prognostic implications.                                                                                                                                                                                             |
| HIGHEST_PX_LEVEL          | LEVEL_Px1, LEVEL_Px2, LEVEL_Px3                                                                                                                                                     | The highest level of evidence for prognostic implications.                                                                                                                                                                       |
| PX_CITATIONS              | PMID, Abstract, Website link                                                                                                                                                        | All citations related to prognostic implications.                                                                                                                                                                                |

## Columns added in the files using ClinicalDataAnnotator
Please see description above for columns LEVEL_&ast;, HIGHEST_LEVEL, HIGHEST_SENSITIVE_LEVEL, HIGHEST_RESISTANCE_LEVEL, LEVEL_Dx*, HIGHEST_DX_LEVEL, LEVEL_Px*, HIGHEST_PX_LEVEL.   
Beside these columsn, the following columns will also be added.

| Column                                              | Description                                                                 |
|-----------------------------------------------------|-----------------------------------------------------------------------------|
| ONCOGENIC_MUTATIONS                                 | The list of mutations that are Oncogenic or Likely Oncogenic.               |
| #ONCOGENIC_MUTATIONS                                | Number of oncogenic mutations.                                              |
| RESISTANCE_MUTATIONS                                | The list of resistance mutations.                                           |
| #RESISTANCE_MUTATIONS                               | Number of resistance mutations.                                             |
| #MUTATIONS_WITH_SENSITIVE_THERAPEUTIC_IMPLICATIONS  | Number of mutations in the sample with sensitive therapeutic implications.  |
| #MUTATIONS_WITH_RESISTANCE_THERAPEUTIC_IMPLICATIONS | Number of mutations in the sample with resistance therapeutic implications. |
| #MUTATIONS_WITH_DIAGNOSTIC_IMPLICATIONS             | Number of mutations in the sample with diagnostic implications.             |
| #MUTATIONS_WITH_PROGNOSTIC_IMPLICATIONS             | Number of mutations in the sample with prognostic implications.             |
| #MUTATIONS                                          | Number of mutations in the sample.                                          |
## Questions?
The best way is to email contact@oncokb.org so all our team members can help.
