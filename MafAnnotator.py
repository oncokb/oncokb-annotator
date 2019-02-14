#!/usr/bin/python

import sys
import getopt
from AnnotatorCore import *

def main(argv):

    inputmaffile = ''
    inputclinicalfile = ''
    outputmaffile = ''
    previousresultfile = ''
    defaultcancertype = 'cancer'
    annotatehotspots = False

    try:
        opts, args = getopt.getopt(argv, "hi:o:p:c:s:t:u:v:a")
    except getopt.GetoptError:
        print 'for help: python MafAnnotator.py -h'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'MafAnnotator.py -i <input MAF file> -o <output MAF file> [-p previous results] [-c <input clinical file>] [-s sample list filter] [-t <default tumor type>] [-u oncokb-base-url] [-a]'
            print '  Essential MAF columns (case insensitive):'
            print '    HUGO_SYMBOL: Hugo gene symbol'
            print '    VARIANT_CLASSIFICATION: Translational effect of variant allele'
            print '    TUMOR_SAMPLE_BARCODE: sample ID'
            print '    AMINO_ACID_CHANGE: amino acid change'
            print '    PROTEIN_START: protein start'
            print '    PROTEIN_END: protein end'
            print '    PROTEIN_POSITION: can be used instead of PROTEIN_START and PROTEIN_END (in the output of vcf2map)'
            print '  Essential clinical columns:'
            print '    SAMPLE_ID: sample ID'
            print '    ONCOTREE_CODE: tumor type code from oncotree (oncotree.mskcc.org)'
            print '  Cancer type will be assigned based on the following priority:'
            print '     1) ONCOTREE_CODE in clinical data file'
            print '     2) ONCOTREE_CODE exist in MAF'
            print '     3) default tumor type (-t)'
            print '  Default OncoKB base url is http://oncokb.org'
            print '  use -a to annotate mutational hotspots'
            sys.exit()
        elif opt in ("-i"):
            inputmaffile = arg
        elif opt in ("-o"):
            outputmaffile = arg
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
        elif opt in ("-a"):
            annotatehotspots = True
        elif opt in ("-v"):
            setcancerhotspotsbaseurl(arg)

    if inputmaffile == '' or outputmaffile=='':
        print 'for help: python MafAnnotator.py -h'
        sys.exit(2)

    cancertypemap = {}
    if inputclinicalfile != '':
        readCancerTypes(inputclinicalfile, cancertypemap)

    print 'annotating '+inputmaffile+"..."

    processalterationevents(inputmaffile, outputmaffile, previousresultfile, defaultcancertype, cancertypemap, False, annotatehotspots)

    print 'done!'

if __name__ == "__main__":
    # argv = [
    #     '-i', 'data/example_maf.txt',
    #     '-o', 'data/example_maf.oncokb.txt',
    #     '-c', 'data/example_clinical.txt',
    #     '-a'
    # ]
    # main(argv)

    # print sys.argv[1:]
    main(sys.argv[1:])