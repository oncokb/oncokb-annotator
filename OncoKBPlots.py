#!/usr/bin/python

import sys
import getopt
from AnnotatorCore import *

def main(argv):

    annotatedclinicalfile = ''
    outputpdffile = ''
    params = {
        "catogerycolumn": "CANCER_TYPE",
        "thresholdcat": 0
    }

    try:
        opts, args = getopt.getopt(argv, "hi:o:c:s:n:")
    except getopt.GetoptError:
        print 'for help: python OncoKBPlots.py -h'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'OncoKBPlots.py -i <annotated clinical file> -o <output PDF file> [-c <categorization column, e.g. CANCER_TYPE>] [-s sample list filter] [-n threshold of # samples in a category]'
            print '  Essential clinical columns:'
            print '    SAMPLE_ID: sample ID'
            print '    HIGHEST_LEVEL: Highest OncoKB levels'
            sys.exit()
        elif opt in ("-i"):
            annotatedclinicalfile = arg
        elif opt in ("-o"):
            outputpdffile = arg
        elif opt in ("-c"):
            params["catogerycolumn"] = arg
        elif opt in ("-s"):
            setsampleidsfileterfile(arg)
        elif opt in ("-n"):
            params["thresholdcat"] = int(arg)

    if annotatedclinicalfile == '' or outputpdffile== '':
        print 'for help: python OncoKBPlots.py -h'
        sys.exit(2)

    print 'annotating ' + annotatedclinicalfile + "..."

    plotclinicalactionability(annotatedclinicalfile, outputpdffile, params)

    print 'done!'

if __name__ == "__main__":
    # argv = [
    #     '-i', '/Users/jgao/projects/oncokb-annotator/process/mskimpact/data_clinical_2017-11-01.oncokb.txt',#'data/example_clinical.oncokb.txt',
    #     '-o', '/Users/jgao/projects/oncokb-annotator/process/mskimpact/data_clinical_2017-11-01.oncokb.pdf',#'data/example_clinical.oncokb.pdf',
    #     '-c', 'CANCER_TYPE',
    #     '-n', '100'
    # ]
    # main(argv)

    # print sys.argv[1:]
    main(sys.argv[1:])