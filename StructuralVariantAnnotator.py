#!/usr/bin/env python

import sys
import argparse
import logging

from AnnotatorCore import setsampleidsfileterfile
from AnnotatorCore import setcancerhotspotsbaseurl
from AnnotatorCore import setoncokbbaseurl
from AnnotatorCore import setoncokbapitoken
from AnnotatorCore import readCancerTypes
from AnnotatorCore import validate_oncokb_token
from AnnotatorCore import process_sv

logging.basicConfig(level=logging.INFO)
log = logging.getLogger('StructuralVariantAnnotator')


def main(argv):
    if argv.help:
        log.info(
            '\n'
            'StructuralVariantAnnotator.py -i <input structural variant file> -o <output structural variant file> '
            '[-p previous results] [-c <input clinical file>] [-s sample list filter] [-t <default tumor type>] '
            '[-u <oncokb api url>] [-b <oncokb api bear token>] [-d include descriptions]\n'
            '  Essential structural variant columns (case insensitive):\n'
            '    GENEA: Hugo gene symbol for gene A\n'
            '    GENEB: Hugo gene symbol for gene B\n'
            '    SV_TYPE: Structural variant type. Available values: DELETION, TRANSLOCATION, DUPLICATION, INSERTION, INVERSION, FUSION, UNKNOWN. Other type will be converted to UNKNOWN\n'
            '    TUMOR_SAMPLE_BARCODE: sample ID\n'
            '  Essential clinical columns:\n'
            '    SAMPLE_ID: sample ID\n'
            '    ONCOTREE_CODE: tumor type code from oncotree (oncotree.mskcc.org)\n'
            '  Cancer type will be assigned based on the following priority:\n'
            '     1) ONCOTREE_CODE in clinical data file\n'
            '     2) ONCOTREE_CODE exist in structural variant\n'
            '     3) default tumor type (-t)\n'
            '  Default OncoKB base url is https://www.oncokb.org'
        )
        sys.exit()
    if argv.input_file == '' or argv.output_file == '' or argv.oncokb_api_bearer_token == '':
        required_params = []
        if argv.input_file == '':
            required_params.append('-i')
        if argv.output_file == '':
            required_params.append('-o')
        if argv.oncokb_api_bearer_token == '':
            required_params.append('-b')

        log.error('The parameter(s) ' + ', '.join(required_params) + ' can not be empty')
        log.info('for help: python StructuralVariantAnnotator.py -h')
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

    validate_oncokb_token()

    log.info('annotating %s ...' % argv.input_file)
    process_sv(argv.input_file, argv.output_file, argv.previous_result_file, argv.default_cancer_type, cancertypemap,
               argv.include_descriptions)

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
    parser.add_argument('-d', dest='include_descriptions', action="store_true", default=False)
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
