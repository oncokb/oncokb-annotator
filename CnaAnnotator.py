#!/usr/bin/python

import sys
import argparse
import logging

from AnnotatorCore import setsampleidsfileterfile
from AnnotatorCore import setoncokbbaseurl
from AnnotatorCore import setoncokbapitoken
from AnnotatorCore import readCancerTypes
from AnnotatorCore import validate_oncokb_token
from AnnotatorCore import process_cna_data
from AnnotatorCore import CNA_FILE_FORMAT_GISTIC

logging.basicConfig(level=logging.INFO)
log = logging.getLogger('CnaAnnotator')


def main(argv):
    if argv.help:
        log.info(
            '\n'
            'CnaAnnotator.py -i <input CNA file> -o <output CNA file> [-p previous results] [-c <input clinical file>] '
            '[-s sample list filter] [-t <default tumor type>] [-u oncokb-base-url] [-b oncokb_api_bear_token] '
            '[-z annotate_gain_loss] [-f CNA file formt, gistic or individual] [-d include descriptions]\n'
            '  Input CNA file uses GISTIC output by default (https://docs.cbioportal.org/5.1-data-loading/data-loading/file-formats#data-file-1). You can also list copy number alteration individually by specifying -f=individual\n'
            '  Essential clinical columns:\n'
            '    SAMPLE_ID: sample ID\n'
            '  Cancer type will be assigned based on the following priority:\n'
            '     1) ONCOTREE_CODE in clinical data file\n'
            '     2) ONCOTREE_CODE exist in MAF\n'
            '     3) default tumor type (-t)\n'
            '  We do not annotate Gain and Loss by default, add -z to include the analysis. See https://github.com/oncokb/oncokb-annotator/issues/51 for more information.\n'
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
        log.info('for help: python CnaAnnotator.py -h')
        sys.exit(2)
    if argv.sample_ids_filter:
        setsampleidsfileterfile(argv.sample_ids_filter)
    if argv.oncokb_api_url:
        setoncokbbaseurl(argv.oncokb_api_url)
    setoncokbapitoken(argv.oncokb_api_bearer_token)

    cancertypemap = {}
    if argv.input_clinical_file:
        readCancerTypes(argv.input_clinical_file, cancertypemap)

    validate_oncokb_token()

    log.info('annotating %s ...' % argv.input_file)
    process_cna_data(argv.input_file, argv.output_file, argv.previous_result_file, argv.default_cancer_type, cancertypemap, argv.include_descriptions, argv.annotate_gain_loss, argv.cna_file_format.lower())

    log.info('done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-h', dest='help', action="store_true", default=False)
    parser.add_argument('-i', dest='input_file', default='', type=str)
    parser.add_argument('-o', dest='output_file', default='', type=str)
    parser.add_argument('-p', dest='previous_result_file', default='', type=str)
    parser.add_argument('-c', dest='input_clinical_file', default='', type=str)
    parser.add_argument('-s', dest='sample_ids_filter', default='', type=str)
    parser.add_argument('-t', dest='default_cancer_type', default='', type=str)
    parser.add_argument('-u', dest='oncokb_api_url', default='', type=str)
    parser.add_argument('-b', dest='oncokb_api_bearer_token', default='', type=str)
    parser.add_argument('-z', dest='annotate_gain_loss', action="store_true", default=False)
    parser.add_argument('-f', dest='cna_file_format', default=CNA_FILE_FORMAT_GISTIC)
    parser.add_argument('-d', dest='include_descriptions', action="store_true", default=False)
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
