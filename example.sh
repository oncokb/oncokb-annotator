#!/usr/bin/env bash
IMAF="data/example_maf.txt"
OMAF="data/example_maf.oncokb.txt"

IMAF38="data/example_maf_grch38.txt"
OMAF38="data/example_maf_grch38.oncokb.txt"

OMAFHGVSPSHORT="data/example_maf_hgvsp_short.oncokb.txt"
OMAFHGVSP="data/example_maf_hgvsp.oncokb.txt"
OMAFHGVSG="data/example_maf_hgvsg.oncokb.txt"
OMAFGC="data/example_maf_genomic_change.oncokb.txt"

IATYPICALALT="data/example_atypical_alterations.txt"
OATYPICALALT="data/example_atypical_alterations.oncokb.txt"

IF="data/example_fusions.txt"
OF="data/example_fusions.oncokb.txt"

ISV="data/example_sv.txt"
OSV="data/example_sv.oncokb.txt"

ICNA="data/example_cna.txt"
OCNA="data/example_cna.oncokb.txt"

IICNA="data/example_individual_cna.txt"
OICNA="data/example_individual_cna.oncokb.txt"

IC="data/example_clinical.txt"
OC="data/example_clinical.oncokb.txt"

TOKEN="d5001618-260c-4ae6-9f46-880fec61c4ff" #OncoKB API Token
README="data/example_README.txt"

if command -v python3 &>/dev/null; then
    PYTHON=python3
elif command -v python &>/dev/null; then
    PYTHON=python
else
    echo "Error: Python not found" >&2
    exit 1
fi

$PYTHON MafAnnotator.py -i "$IMAF" -o "$OMAF" -c "$IC" -b "$TOKEN"
$PYTHON MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSPSHORT" -c "$IC" -b "$TOKEN" -q hgvsp_short
$PYTHON MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSP" -c "$IC" -b "$TOKEN" -q hgvsp
$PYTHON MafAnnotator.py -i "$IMAF" -o "$OMAFHGVSG" -c "$IC" -b "$TOKEN" -q hgvsg
$PYTHON MafAnnotator.py -i "$IMAF" -o "$OMAFGC" -c "$IC" -b "$TOKEN" -q genomic_change

$PYTHON MafAnnotator.py -i "$IMAF38" -o "$OMAF38" -c "$IC" -b "$TOKEN"

$PYTHON MafAnnotator.py -i "$IATYPICALALT" -o "$OATYPICALALT" -c "$IC" -b "$TOKEN"

$PYTHON FusionAnnotator.py -i "$IF" -o "$OF" -c "$IC" -b "$TOKEN"
$PYTHON StructuralVariantAnnotator.py -i "$ISV" -o "$OSV" -c "$IC" -b "$TOKEN"
$PYTHON CnaAnnotator.py -i "$ICNA" -o "$OCNA" -c "$IC" -b "$TOKEN"
$PYTHON CnaAnnotator.py -i "$IICNA" -o "$OICNA" -c "$IC" -b "$TOKEN" -f "individual" -z
$PYTHON ClinicalDataAnnotator.py -i "$IC" -o "$OC" -a "$OMAF,$OATYPICALALT,$OCNA,$OF,$OSV"

$PYTHON GenerateReadMe.py -o "$README"
