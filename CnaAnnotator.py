#!/usr/bin/python

import sys
import getopt
from AnnotatorCore import *

def main(argv):

    baseurl = 'http://oncokb.org'
    inputcnafile = ''
    inputclinicalfile = ''
    outputcnafile = ''
    previousresultfile = ''
    defaultcancertype = 'cancer'

    try:
        opts, args = getopt.getopt(argv, "hi:o:p:c:s:d:t:u:")
    except getopt.GetoptError:
        print 'for help: python CnaAnnotator.py -h'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'CnaAnnotator.py -i <input cNA file> -o <output MAF file> [-p previous results] [-c <input clinical file>] [-s sample list filter] [-t <default tumor type>] [-u base-url]'
            print '  Input CNA file should follow the GISTIC output (https://cbioportal.readthedocs.io/en/latest/File-Formats.html#discrete-copy-number-data)'
            print '  Essential clinical columns:'
            print '    SAMPLE_ID: sample ID'
            print '  Cancer type will be assigned based on the following priority:'
            print '     1) ONCOTREE_CODE in clinical data file'
            print '     2) ONCOTREE_CODE exist in MAF'
            print '     3) default tumor type (-t)'
            print '  Default OncoKB base url is http://oncokb.org'
            sys.exit()
        elif opt in ("-i"):
            inputcnafile = arg
        elif opt in ("-o"):
            outputcnafile = arg
        elif opt in ("-p"):
            previousresultfile = arg
        elif opt in ("-c"):
            inputclinicalfile = arg
        elif opt in ("-s"):
            setsampleidsfileterfile(arg)
        elif opt in ("-t"):
            defaultcancertype = arg
        elif opt in ("-u"):
            setoncokbbaseurl(arg)

    if inputcnafile == '' or outputcnafile=='':
        print 'for help: python MafAnnotator.py -h'
        sys.exit(2)

    cancertypemap = {}
    if inputclinicalfile != '':
        readCancerTypes(inputclinicalfile, cancertypemap)

    print 'annotating '+inputcnafile+"..."

    processcnagisticdata(inputcnafile, outputcnafile, previousresultfile, defaultcancertype, cancertypemap, False)

    print 'done!'

if __name__ == "__main__":
    # argv = [
    #     '-i', 'data/example_cna.txt',
    #     '-o', 'data/example_cna.oncokb.txt',
    #     '-c', 'data/example_clinical.txt',
    # ]
    # main(argv)

    # print sys.argv[1:]
    main(sys.argv[1:])
