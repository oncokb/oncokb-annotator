#!/usr/bin/python

import sys
import re
import argparse
import logging
import os
import csv
import matplotlib.pyplot as plt

from AnnotatorCore import setsampleidsfileterfile
from AnnotatorCore import readheaders
from AnnotatorCore import geIndexOfHeader
from AnnotatorCore import sampleidsfilter
from AnnotatorCore import levels
from AnnotatorCore import dxLevels
from AnnotatorCore import pxLevels
from AnnotatorCore import SAMPLE_HEADERS

logging.basicConfig(level=logging.INFO)
log = logging.getLogger('OncoKBPlots')


def plotclinicalactionability(ax, annotatedclinicalfile, outfile, parameters):
    if os.path.isfile(outfile):
        os.remove(outfile)

    extlevels = levels + ["ONCOGENIC", "VUS"]
    if "levels" in parameters:
        extlevels = parameters["levels"]

    with open(annotatedclinicalfile, 'rU') as clinfile:
        reader = csv.reader(clinfile, delimiter='\t')
        headers = readheaders(reader)
        isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
        ilevel = headers['HIGHEST_LEVEL']
        ioncogenic = headers['ONCOGENIC_MUTATIONS']
        icat = headers[parameters["catogerycolumn"].upper()]  # e.g. "CANCER_TYPE"

        catsamplecount = {}
        catactionablesamplecount = {}
        oncogenicsamplecount = {}
        levelcatsamplecount = {}

        for row in reader:
            sample = row[isample]
            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            cat = row[icat]
            if cat not in catsamplecount:
                catsamplecount[cat] = 0
            catsamplecount[cat] += 1

            if cat not in catactionablesamplecount:
                catactionablesamplecount[cat] = 0
                oncogenicsamplecount[cat] = 0

            level = row[ilevel]
            oncogenic = row[ioncogenic]

            exlevel = level

            if level in extlevels:
                catactionablesamplecount[cat] += 1
                oncogenicsamplecount[cat] += 1
            elif len(oncogenic.strip()) > 0:
                oncogenicsamplecount[cat] += 1
                exlevel = "ONCOGENIC"
            else:
                exlevel = "VUS"

            if exlevel not in levelcatsamplecount:
                levelcatsamplecount[exlevel] = {}
            if cat not in levelcatsamplecount[exlevel]:
                levelcatsamplecount[exlevel][cat] = 0
            levelcatsamplecount[exlevel][cat] += 1

    # plot
    catarray = []  # cancer types
    catactionabilityarray = []  # actionabiligy percentages per cancer type
    catoncogenicarray = []  # actionabiligy percentages per cancer type
    for cat in catsamplecount:
        if catsamplecount[cat] >= parameters["thresholdcat"]:
            catarray.append(cat)
            catactionabilityarray.append(catactionablesamplecount[cat] * 100.0 / catsamplecount[cat])
            catoncogenicarray.append(oncogenicsamplecount[cat] * 100.0 / catsamplecount[cat])

    ncat = len(catarray)
    order = reversed(sorted(range(ncat), key=lambda x: (catactionabilityarray[x], catoncogenicarray[x])))
    drawplot(ax, 'OncoKB Actionability', extlevels, levelcatsamplecount, catarray, catsamplecount, order,
             parameters["thresholdcat"])


def plotimplications(ax, header, title, levels, annotatedclinicalfile, outfile, parameters):
    if os.path.isfile(outfile):
        os.remove(outfile)

    extlevels = levels
    if "levels" in parameters:
        extlevels = parameters["levels"]

    with open(annotatedclinicalfile, 'rU') as clinfile:
        reader = csv.reader(clinfile, delimiter='\t')
        headers = readheaders(reader)
        isample = headers['SAMPLE_ID']
        ilevel = headers[header]
        icat = headers[parameters["catogerycolumn"].upper()]

        catsamplecount = {}
        catactionablesamplecount = {}
        levelcatsamplecount = {}

        for row in reader:
            sample = row[isample]
            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            cat = row[icat]
            if cat not in catsamplecount:
                catsamplecount[cat] = 0
            catsamplecount[cat] += 1

            if cat not in catactionablesamplecount:
                catactionablesamplecount[cat] = 0

            level = row[ilevel]

            exlevel = level

            if level in extlevels:
                catactionablesamplecount[cat] += 1
            else:
                exlevel = "Other"

            if exlevel not in levelcatsamplecount:
                levelcatsamplecount[exlevel] = {}
            if cat not in levelcatsamplecount[exlevel]:
                levelcatsamplecount[exlevel][cat] = 0
            levelcatsamplecount[exlevel][cat] += 1

    # plot
    catarray = []  # cancer types
    catactionabilityarray = []  # actionabiligy percentages per cancer type
    for cat in catsamplecount:
        if catsamplecount[cat] >= parameters["thresholdcat"]:
            catarray.append(cat)
            catactionabilityarray.append(catactionablesamplecount[cat] * 100.0 / catsamplecount[cat])

    ncat = len(catarray)
    order = reversed(sorted(range(ncat), key=lambda x: (catactionabilityarray[x])))
    drawplot(ax, title, extlevels, levelcatsamplecount, catarray, catsamplecount, order, parameters["thresholdcat"])


def drawplot(ax, title, extlevels, levelcatsamplecount, catarray, catsamplecount, order, thresholdcat):
    # level colors
    levelcolors = {
        'LEVEL_1': '#33A02C',
        'LEVEL_2': '#1F78B4',
        'LEVEL_3A': '#984EA3',
        'LEVEL_3B': '#BE98CE',
        'LEVEL_4': '#a8a8a8',
        'LEVEL_R1': '#EE3424',
        'LEVEL_R2': '#F79A92',

        'LEVEL_Dx1': '#33A02C',
        'LEVEL_Dx2': '#1F78B4',
        'LEVEL_Dx3': '#984EA3',

        'LEVEL_Px1': '#33A02C',
        'LEVEL_Px2': '#1F78B4',
        'LEVEL_Px3': '#984EA3',

        'ONCOGENIC': '#ffdab9',
        'VUS': '#d1d1d1',
        'Other': 'grey'
    }

    # level legend
    levellegend = {
        'LEVEL_1': 'Level 1',
        'LEVEL_2': 'Level 2',
        'LEVEL_3A': 'Level 3A',
        'LEVEL_3B': 'Level 3B',
        'LEVEL_4': 'Level 4',
        'LEVEL_R1': 'Level R1',
        'LEVEL_R2': 'Level R2',

        'LEVEL_Dx1': 'Level Dx1',
        'LEVEL_Dx2': 'Level Dx2',
        'LEVEL_Dx3': 'Level Dx3',

        'LEVEL_Px1': 'Level Px1',
        'LEVEL_Px2': 'Level Px2',
        'LEVEL_Px3': 'Level Px3',

        'ONCOGENIC': 'Oncogenic, no level',
        'VUS': 'VUS',
        'Other': 'Other'
    }

    ncat = len(catarray)
    if ncat > 0:
        catarray = [catarray[i] for i in order]

        ind = range(ncat)

        legends = []
        plts = []
        accumlevelcancerperc = [0] * ncat
        for level in extlevels:
            if level not in levelcatsamplecount:
                continue

            levelcancerperc = [0] * ncat
            for k in ind:
                cat = catarray[k]
                if catsamplecount[cat] < thresholdcat:
                    continue
                if cat in levelcatsamplecount[level]:
                    levelcancerperc[k] = levelcatsamplecount[level][cat] * 100.0 / catsamplecount[cat]

            width = 0.75
            plts = [ax.bar(ind, levelcancerperc, width, color=levelcolors[level], bottom=accumlevelcancerperc)] + plts
            legends = [levellegend[level]] + legends
            accumlevelcancerperc = list(map(sum, zip(accumlevelcancerperc, levelcancerperc)))

        ax = plt.gca()
        ax.set_axisbelow(True)
        ax.set_aspect(0.1)

        ax.tick_params(axis='y', which='major', labelsize=6)
        ax.set_ylabel('% of samples', fontsize=6)
        ax.set_title(title, fontsize=8)
        ax.set_xticks([i + 0.5 for i in ind])
        ax.set_xticklabels(catarray, rotation=60, ha="right", fontsize=4)
        # plt.yticks(np.arange(0, 81, 10))
        ax.legend(plts, legends, fontsize=6, bbox_to_anchor=(1.01, 1), loc="upper left")


def main(argv):
    params = {
        "catogerycolumn": argv.catogery_column,  # -c
        "thresholdcat": argv.threshold_cat,  # -n
    }
    if argv.help:
        log.info(
            '\n'
            'OncoKBPlots.py -i <annotated clinical file> -o <output PDF file> [-c <categorization column, '
            'e.g. CANCER_TYPE>] [-s sample list filter] [-n threshold of # samples in a category] [-l comma separated levels to include]\n'
            '  Essential clinical columns:\n'
            '    SAMPLE_ID: sample ID\n'
            '    HIGHEST_LEVEL: Highest OncoKB levels\n'
            '  Supported levels (-l): \n'
            '    LEVEL_1,LEVEL_2,LEVEL_3A,LEVEL_3B,LEVEL_4,ONCOGENIC,VUS'
        )
        sys.exit()
    if argv.input_file == '' or argv.output_file == '':
        required_params = []
        if argv.input_file == '':
            required_params.append('-i')
        if argv.output_file == '':
            required_params.append('-o')

        log.error('The parameter(s) ' + ', '.join(required_params) + ' can not be empty')
        log.info('for help: python OncoKBPlots.py -h')
        sys.exit(2)
    if argv.sample_ids_filter:
        setsampleidsfileterfile(argv.sample_ids_filter)
    if argv.levels:
        params["levels"] = re.split(',', argv.levels)

    log.info('annotating %s ...' % argv.input_file)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    plotclinicalactionability(ax1, argv.input_file, argv.output_file, params)

    # ax.yaxis.grid(linestyle="dotted", color="lightgray") # horizontal lines
    # plt.margins(0.01)

    plotclinicalactionability(ax1, args.input_file, args.output_file, params)
    plotimplications(ax2, 'HIGHEST_DX_LEVEL', 'OncoKB Diagnostic Implications', dxLevels, args.input_file,
                     argv.output_file, params)
    plotimplications(ax3, 'HIGHEST_PX_LEVEL', 'OncoKB Prognostic Implications', pxLevels, args.input_file,
                     argv.output_file, params)

    plt.subplots_adjust(left=0.2, bottom=0.3)
    plt.gcf().text(0.90, 0.1, "Generated by OncoKB\n[Chakravarty et al., JCO PO 2017]", fontsize=6,
                   horizontalalignment='right', verticalalignment='bottom')
    fig.tight_layout()
    fig.savefig(argv.output_file, bbox_inches='tight')

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
