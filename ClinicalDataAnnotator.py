#!/usr/bin/python

import sys
import argparse
from AnnotatorCore import *


def main(argv):
    if argv.help:
        print 'ClinicalDataAnnotator.py -i <input clinical file> -o <output clinical file> -a <annotated alteration files, separate by ,> [-s sample list filter]'
        print '  Essential clinical columns:'
        print '    SAMPLE_ID: sample ID'
        sys.exit()
    if argv.sample_ids_filter:
        setsampleidsfileterfile(argv.sample_ids_filter)

    annotated_alteration_files = re.split(',|, ', argv.annotated_alteration_files)
    if argv.input_file == '' or argv.output_file == '' or len(annotated_alteration_files) == 0:
        print 'for help: python ClinicalDataAnnotator.py -h'
        sys.exit(2)

    print 'annotating %s ...' % argv.input_file
    processclinicaldata(annotated_alteration_files, argv.input_file, argv.output_file)

    print 'done!'


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-h', dest='help', action="store_true", default=False)
    parser.add_argument('-i', dest='input_file', default='', type=str)
    parser.add_argument('-o', dest='output_file', default='', type=str)
    parser.add_argument('-s', dest='sample_ids_filter', default='', type=str)
    parser.add_argument('-a', dest='annotated_alteration_files', default='', type=str)
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
