#!/usr/bin/python

import sys
import argparse
from AnnotatorCore import *


def main(argv):
    if argv.help:
        print 'MafAnnotator.py -i <input MAF file> -o <output MAF file> [-p previous results] [-c <input clinical file>] [-s sample list filter] [-t <default tumor type>] [-u oncokb-base-url] [-b oncokb_api_bear_token] [-a]'
        print '  Essential MAF columns (case insensitive):'
        print '    HUGO_SYMBOL: Hugo gene symbol'
        print '    VARIANT_CLASSIFICATION: Translational effect of variant allele'
        print '    TUMOR_SAMPLE_BARCODE: sample ID'
        print '    AMINO_ACID_CHANGE: amino acid change'
        print '    PROTEIN_START: protein start'
        print '    PROTEIN_END: protein end'
        print '    PROTEIN_POSITION: can be used instead of PROTEIN_START and PROTEIN_END (in the output of vcf2map)'
        print '  Essential clinical columns:'
        print '    SAMPLE_ID: sample ID'
        print '    ONCOTREE_CODE: tumor type code from oncotree (oncotree.mskcc.org)'
        print '  Cancer type will be assigned based on the following priority:'
        print '     1) ONCOTREE_CODE in clinical data file'
        print '     2) ONCOTREE_CODE exist in MAF'
        print '     3) default tumor type (-t)'
        print '  Default OncoKB base url is http://oncokb.org'
        print '  use -a to annotate mutational hotspots'
        sys.exit()
    if argv.input_file == '' or argv.output_file == '' or argv.oncokb_api_bearer_token == '':
        print 'for help: python MafAnnotator.py -h'
        sys.exit(2)

    if argv.sample_ids_filter:
        setsampleidsfileterfile(argv.sample_ids_filter)
    if argv.cancer_hotspots_base_url:
        setcancerhotspotsbaseurl(argv.cancer_hotspots_base_url)
    if argv.oncokb_api_url:
        setoncokbbaseurl(argv.oncokb_api_url)
    if argv.oncokb_api_bearer_token:
        setoncokbapitoken(argv.oncokb_api_bearer_token)

    cancertypemap = {}
    if argv.input_clinical_file:
        readCancerTypes(argv.input_clinical_file, cancertypemap)

    print 'annotating %s ...' % argv.input_file
    processalterationevents(argv.input_file, argv.output_file, argv.previous_result_file, argv.default_cancer_type, cancertypemap, False, argv.annotate_hotspots)

    print 'done!'


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-h', dest='help', action="store_true", default=False)
    parser.add_argument('-i', dest='input_file', default='', type=str)
    parser.add_argument('-o', dest='output_file', default='', type=str)
    parser.add_argument('-p', dest='previous_result_file', default='', type=str)
    parser.add_argument('-c', dest='input_clinical_file', default='', type=str)
    parser.add_argument('-s', dest='sample_ids_filter', default='', type=str)
    parser.add_argument('-t', dest='default_cancer_type', default='cancer', type=str)
    parser.add_argument('-u', dest='oncokb_api_url', default='', type=str)
    parser.add_argument('-a', dest='annotate_hotspots', action="store_true", default=False)
    parser.add_argument('-v', dest='cancer_hotspots_base_url', default='', type=str)
    parser.add_argument('-b', dest='oncokb_api_bearer_token', default='', type=str)
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
