#!/usr/bin/python

import argparse
from AnnotatorCore import *
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('FusionAnnotator')

def main(argv):
    if argv.help:
        log.info('\n'
        'FusionAnnotator.py -i <input Fusion file> -o <output Fusion file> [-p previous results] [-c <input clinical file>] [-s sample list filter] [-t <default tumor type>] [-u <oncokb api url>] [-b <oncokb api bear token>] [-r <structural variant name format, default: [A-Za-z\d]+-[A-Za-z\d]+>]\n'
        '  Essential Fusion columns (case insensitive):\n'
        '    HUGO_SYMBOL: Hugo gene symbol\n'
        '    VARIANT_CLASSIFICATION: Translational effect of variant allele\n'
        '    TUMOR_SAMPLE_BARCODE: sample ID\n'
        '    FUSION: amino acid change, e.g. "TMPRSS2-ERG"\n'
        '  Essential clinical columns:\n'
        '    SAMPLE_ID: sample ID\n'
        '    ONCOTREE_CODE: tumor type code from oncotree (oncotree.mskcc.org)\n'
        '  Cancer type will be assigned based on the following priority:\n'
        '     1) ONCOTREE_CODE in clinical data file\n'
        '     2) ONCOTREE_CODE exist in Fusion\n'
        '     3) default tumor type (-t)\n'
        '  Default OncoKB base url is https://www.oncokb.org')
        sys.exit()
    if argv.input_file == '' or argv.output_file == '' or argv.oncokb_api_bearer_token == '':
        log.info('for help: python FusionAnnotator.py -h')
        sys.exit(2)
    if argv.sample_ids_filter:
        setsampleidsfileterfile(argv.sample_ids_filter)
    if argv.cancer_hotspots_base_url:
        setcancerhotspotsbaseurl(argv.cancer_hotspots_base_url)
    if argv.oncokb_api_url:
        setoncokbbaseurl(argv.oncokb_api_url)
    setoncokbapitoken(argv.oncokb_api_bearer_token)

    cancertypemap = {}
    if argv.input_clinical_file:
        readCancerTypes(argv.input_clinical_file, cancertypemap)

    log.info('annotating %s ...' % argv.input_file)
    processsv(argv.input_file, argv.output_file, argv.previous_result_file, argv.default_cancer_type,
              cancertypemap, argv.structural_variant_name_format)

    log.info('done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    # ArgumentParser doesn't accept "store_true" and "type=" at the same time.
    parser.add_argument('-h', dest='help', action="store_true", default=False)
    parser.add_argument('-i', dest='input_file', default='', type=str)
    parser.add_argument('-o', dest='output_file', default='', type=str)
    parser.add_argument('-p', dest='previous_result_file', default='', type=str)
    parser.add_argument('-c', dest='input_clinical_file', default='', type=str)
    parser.add_argument('-s', dest='sample_ids_filter', default=None, type=str)
    parser.add_argument('-t', dest='default_cancer_type', default='', type=str)
    parser.add_argument('-u', dest='oncokb_api_url', default='', type=str)
    parser.add_argument('-v', dest='cancer_hotspots_base_url', default='', type=str)
    parser.add_argument('-b', dest='oncokb_api_bearer_token', default='', type=str)
    parser.add_argument('-r', dest='structural_variant_name_format', default=None, type=str)
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
