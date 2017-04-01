#!/usr/bin/env bash
IMAF="data/example_maf.txt"
OMAF="data/example_maf.oncokb.txt"
ICNA="data/example_cna.txt"
OCNA="data/example_cna.oncokb.txt"
IC="data/example_clinical.txt"
OC="data/example_clinical.oncokb.txt"
python MafAnnotator.py -i $IMAF -o $OMAF -c $IC
python CnaAnnotator.py -i $ICNA -o $OCNA -c $IC
python ClinicalDataAnnotator.py -i $IC -o $OC -a $OMAF,$OCNA

