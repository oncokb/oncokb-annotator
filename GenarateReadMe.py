#!/usr/bin/python

import sys
import getopt
from AnnotatorCore import *

def main(argv):

    outputfile = ''

    try:
        opts, args = getopt.getopt(argv, "ho:u:")
    except getopt.GetoptError:
        print 'for help: python GenerateReadMe.py -h'
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print 'GenerateReadMe.py -o <output README file> [-u oncokb-base-url]'
            print '  Default OncoKB base url is http://oncokb.org'
            sys.exit()
        elif opt in ("-o"):
            outputfile = arg
        elif opt in ("-u"):
            setoncokbbaseurl(arg)

    if outputfile=='':
        print 'for help: python GenerateReadMe.py -h'
        sys.exit(2)

    generateReadme(outputfile)

    print 'done!'

if __name__ == "__main__":
    # argv = [
    #     '-o', 'data/README.txt'
    # ]
    # main(argv)

    # print sys.argv[1:]
    main(sys.argv[1:])