#!/usr/bin/env bash
IMAF="data/example_maf.txt"
OMAF="data/example_maf.oncokb.txt"
IF="data/example_fusions.txt"
OF="data/example_fusions.oncokb.txt"
ICNA="data/example_cna.txt"
OCNA="data/example_cna.oncokb.txt"
IC="data/example_clinical.txt"
OC="data/example_clinical.oncokb.txt"
OCPDF="data/example_clinical.oncokb.pdf"
TOKEN="" #OncoKB API Token
README="data/example_README.txt"
python MafAnnotator.py -i $IMAF -o $OMAF -c $IC -b $TOKEN
python FusionAnnotator.py -i $IF -o $OF -c $IC -b $TOKEN
python CnaAnnotator.py -i $ICNA -o $OCNA -c $IC -b $TOKEN
python ClinicalDataAnnotator.py -i $IC -o $OC -a $OMAF,$OCNA,$OF
python OncoKBPlots.py -i $OC -o $OCPDF -c ONCOTREE_CODE #-n 10
python GenerateReadMe.py -o $README
