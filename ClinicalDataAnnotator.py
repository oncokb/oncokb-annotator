#!/usr/bin/python

import sys
import getopt
from AnnotatorCore import *

def main(argv):

    baseurl = 'http://oncokb.org'
    inputclinicalfile = ''
    outputclinicalfile = ''
    annotatedalterationfiles = []

    try:
        opts, args = getopt.getopt(argv, "hi:o:a:")
    except getopt.GetoptError:
        print 'for help: python ClinicalDataAnnotator.py -h'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'ClinicalDataAnnotator.py -i <input clinical file> -o <output clinical file> -a <annotated alteration files, separate by ,>'
            print '  Essential clinical columns:'
            print '    SAMPLE_ID: sample ID'
            sys.exit()
        elif opt in ("-i"):
            inputclinicalfile = arg
        elif opt in ("-o"):
            outputclinicalfile = arg
        elif opt in ("-a"):
            annotatedalterationfiles = arg.split(',')

    if inputclinicalfile == '' or outputclinicalfile=='' or len(annotatedalterationfiles)==0:
        print 'for help: python ClinicalDataAnnotator.py -h'
        sys.exit(2)

    print 'annotating '+inputclinicalfile+"..."

    processclinicaldata(annotatedalterationfiles, inputclinicalfile, outputclinicalfile)

    print 'done!'

if __name__ == "__main__":
    # argv = [
    #     '-i', 'data/example_clinical.txt',
    #     '-o', 'data/example_clinical.oncokb.txt',
    #     '-a', 'data/example_maf.oncokb.txt,data/example_cna.oncokb.txt'
    # ]
    # main(argv)

    # print sys.argv[1:]
    main(sys.argv[1:])