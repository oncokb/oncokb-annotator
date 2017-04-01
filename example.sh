#!/usr/bin/env bash
python MafAnnotator.py -i data/example_maf.txt -o data/example_maf.oncokb.txt -c data/example_clinical.txt
python CnaAnnotator.py -i data/example_cna.txt -o data/example_cna.oncokb.txt -c data/example_clinical.txt
python ClinicalDataAnnotator.py -i data/example_clinical.txt -o data/example_clinical.oncokb.txt -a data/example_maf.oncokb.txt#data/example_cna.oncokb.txt

