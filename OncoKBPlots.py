#!/usr/bin/python

import argparse
from AnnotatorCore import *
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('OncoKBPlots')

def main(argv):
    params = {
        "catogerycolumn": argv.catogery_column,  # -c
        "thresholdcat": argv.threshold_cat,  # -n
    }
    if argv.help:
        log.info('\n'
        'OncoKBPlots.py -i <annotated clinical file> -o <output PDF file> [-c <categorization column, '
        'e.g. CANCER_TYPE>] [-s sample list filter] [-n threshold of # samples in a category] [-l comma separated levels to include]\n'
        '  Essential clinical columns:\n'
        '    SAMPLE_ID: sample ID\n'
        '    HIGHEST_LEVEL: Highest OncoKB levels\n'
        '  Supported levels (-l): \n'
        '    LEVEL_1,LEVEL_2,LEVEL_3A,LEVEL_3B,LEVEL_4,ONCOGENIC,VUS')
        sys.exit()
    if argv.input_file == '' or argv.output_file == '':
        log.info('for help: python OncoKBPlots.py -h')
        sys.exit(2)
    if argv.sample_ids_filter:
        setsampleidsfileterfile(argv.sample_ids_filter)
    if argv.levels:
        params["levels"] = re.split(',', argv.levels)

    log.info('annotating %s ...' % argv.input_file)
    plotclinicalactionability(argv.input_file, argv.output_file, params)

    log.info('done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-h', dest='help', action="store_true", default=False)
    parser.add_argument('-i', dest='input_file', default='', type=str)
    parser.add_argument('-o', dest='output_file', default='', type=str)
    parser.add_argument('-c', dest='catogery_column', default='CANCER_TYPE', type=str)
    parser.add_argument('-s', dest='sample_ids_filter', default='', type=str)
    parser.add_argument('-n', dest='threshold_cat', default=0, type=int)
    parser.add_argument('-l', dest='levels', default='', type=str)
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
