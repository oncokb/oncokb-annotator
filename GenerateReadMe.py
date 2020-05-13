#!/usr/bin/python

import argparse
from AnnotatorCore import *
import logging
logging.basicConfig(level=logging.INFO)
log = logging.getLogger('GenerateReadMe')


def main(argv):
    if argv.help:
        log.info('\nGenerateReadMe.py -o <output README file> [-u oncokb-base-url]\n'
                 '  Default OncoKB base url is https://www.oncokb.org')
        sys.exit()
    if argv.output_file == '':
        log.info('for help: python GenerateReadMe.py -h')
        sys.exit(2)
    if argv.oncokb_api_url:
        setoncokbbaseurl(argv.oncokb_api_url)

    generateReadme(argv.output_file)
    log.info('done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    # ArgumentParser doesn't accept "store_true" and "type=" at the same time.
    parser.add_argument('-h', dest='help', action="store_true", default=False)
    parser.add_argument('-o', dest='output_file', default='', type=str)
    parser.add_argument('-u', dest='oncokb_api_url', default='', type=str)
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
