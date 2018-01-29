#!/usr/bin/python

import sys
import getopt
from AnnotatorCore import *

def main(argv):

    inputfusionfile = ''
    inputclinicalfile = ''
    outputfusionfile = ''
    previousresultfile = ''
    defaultcancertype = 'cancer'

    try:
        opts, args = getopt.getopt(argv, "hi:o:p:c:s:t:u:")
    except getopt.GetoptError:
        print 'for help: python FusionAnnotator.py -h'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'FusionAnnotator.py -i <input Fusion file> -o <output Fusion file> [-p previous results] [-c <input clinical file>] [-s sample list filter] [-t <default tumor type>] [-u oncokb-base-url]'
            print '  Essential Fusion columns (case insensitive):'
            print '    HUGO_SYMBOL: Hugo gene symbol'
            print '    VARIANT_CLASSIFICATION: Translational effect of variant allele'
            print '    TUMOR_SAMPLE_BARCODE: sample ID'
            print '    FUSION: amino acid change, e.g. "TMPRSS2-ERG fusion"'
            print '  Essential clinical columns:'
            print '    SAMPLE_ID: sample ID'
            print '    ONCOTREE_CODE: tumor type code from oncotree (oncotree.mskcc.org)'
            print '  Cancer type will be assigned based on the following priority:'
            print '     1) ONCOTREE_CODE in clinical data file'
            print '     2) ONCOTREE_CODE exist in Fusion'
            print '     3) default tumor type (-t)'
            print '  Default OncoKB base url is http://oncokb.org'
            sys.exit()
        elif opt in ("-i"):
            inputfusionfile = arg
        elif opt in ("-o"):
            outputfusionfile = arg
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
        elif opt in ("-v"):
            setcancerhotspotsbaseurl(arg)

    if inputfusionfile == '' or outputfusionfile=='':
        print 'for help: python FusionAnnotator.py -h'
        sys.exit(2)

    cancertypemap = {}
    if inputclinicalfile != '':
        readCancerTypes(inputclinicalfile, cancertypemap)

    print 'annotating '+inputfusionfile+"..."

    processsv(inputfusionfile, outputfusionfile, previousresultfile, defaultcancertype, cancertypemap, False)

    print 'done!'

if __name__ == "__main__":
    # argv = [
    #     '-i', 'data/example_fusions.txt',
    #     '-o', 'data/example_fusions.oncokb.txt',
    #     '-c', 'data/example_clinical.txt',
    # ]
    # main(argv)

    # print sys.argv[1:]
    main(sys.argv[1:])