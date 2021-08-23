#!/usr/bin/env bash
IMAF="data/example_maf.txt"
OMAF="data/example_maf.oncokb.txt"

IMAF38="data/example_maf_grch38.txt"
OMAF38="data/example_maf_grch38.oncokb.txt"

OMAFHGVSPSHORT="data/example_maf_hgsp_short.oncokb.txt"
OMAFHGVSP="data/example_maf_hgsp.oncokb.txt"
OMAFHGVSG="data/example_maf_hgsg.oncokb.txt"
OMAFGC="data/example_maf_genomic_change.oncokb.txt"

IATYPICALALT="data/example_atypical_alterations.txt"
OATYPICALALT="data/example_atypical_alterations.oncokb.txt"

IF="data/example_fusions.txt"
OF="data/example_fusions.oncokb.txt"

ICNA="data/example_cna.txt"
OCNA="data/example_cna.oncokb.txt"

IC="data/example_clinical.txt"
OC="data/example_clinical.oncokb.txt"

OCPDF="data/example_clinical.oncokb.pdf"

TOKEN="" #OncoKB API Token
README="data/example_README.txt"

python MafAnnotator.py -i "$IMAF" -o "$OMAF" -c "$IC" -b "$TOKEN"
python MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSPSHORT" -c "$IC" -b "$TOKEN" -q hgvsp_short
python MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSP" -c "$IC" -b "$TOKEN" -q hgvsp
python MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSG" -c "$IC" -b "$TOKEN" -q hgvsg
python MafAnnotator.py -i "$IMAF" -o "$OMAFGC" -c "$IC" -b "$TOKEN" -q genomic_change

python MafAnnotator.py -i "$IMAF38" -o "$OMAF38" -c "$IC" -b "$TOKEN"

python MafAnnotator.py -i "$IATYPICALALT" -o "$OATYPICALALT" -c "$IC" -b "$TOKEN"

python FusionAnnotator.py -i "$IF" -o "$OF" -c "$IC" -b "$TOKEN"
python CnaAnnotator.py -i "$ICNA" -o "$OCNA" -c "$IC" -b "$TOKEN"
python ClinicalDataAnnotator.py -i "$IC" -o "$OC" -a "$OMAF,$OATYPICALALT,$OCNA,$OF"
python OncoKBPlots.py -i "$OC" -o "$OCPDF" -c ONCOTREE_CODE #-n 10
python GenerateReadMe.py -o "$README"
