#!/usr/bin/python

import sys
import getopt
import csv
import json
import urllib
import os.path
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

csv.field_size_limit(sys.maxsize) # for reading large files

oncokbbaseurl = "http://oncokb.org"
def setoncokbbaseurl(u):
    global oncokbbaseurl
    oncokbbaseurl = u

cancerhotspotsbaseurl = "http://www.cancerhotspots.org"
def setcancerhotspotsbaseurl(u):
    global cancerhotspotsbaseurl
    cancerhotspotsbaseurl = u

_3dhotspotsbaseurl = "http://www.3dhotspots.org"
def set3dhotspotsbaseurl(u):
    global _3dhotspotsbaseurl
    _3dhotspotsbaseurl = u

sampleidsfilter = None
def setsampleidsfileterfile(f):
    global sampleidsfilter
    content = [line.rstrip() for line in open(f)]
    sampleidsfilter = set(content)
    print len(sampleidsfilter)

levels = [
    'LEVEL_1',
    'LEVEL_2A',
    'LEVEL_2B',
    'LEVEL_3A',
    'LEVEL_3B',
    'LEVEL_4',
    'LEVEL_R1'
]

mutationtypeconsequencemap = {
    '3\'Flank': ['any'],
    '5\'Flank ': ['any'],
    'Targeted_Region': ['inframe_deletion', 'inframe_insertion'],
    'COMPLEX_INDEL': ['inframe_deletion', 'inframe_insertion'],
    'ESSENTIAL_SPLICE_SITE': ['feature_truncation'],
    'Exon skipping': ['inframe_deletion'],
    'Frameshift deletion': ['frameshift_variant'],
    'Frameshift insertion': ['frameshift_variant'],
    'FRAMESHIFT_CODING': ['frameshift_variant'],
    'Frame_Shift_Del': ['frameshift_variant'],
    'Frame_Shift_Ins': ['frameshift_variant'],
    'Fusion': ['fusion'],
    'Indel': ['frameshift_variant', 'inframe_deletion', 'inframe_insertion'],
    'In_Frame_Del': ['inframe_deletion'],
    'In_Frame_Ins': ['inframe_insertion'],
    'Missense': ['missense_variant'],
    'Missense_Mutation': ['missense_variant'],
    'Nonsense_Mutation': ['stop_gained'],
    'Nonstop_Mutation': ['stop_lost'],
    'Splice_Site': ['splice_region_variant'],
    'Splice_Site_Del': ['splice_region_variant'],
    'Splice_Site_SNP': ['splice_region_variant'],
    'splicing': ['splice_region_variant'],
    'Translation_Start_Site': ['start_lost'],
    'vIII deletion': ['any']
}


def getsampleid(rawsampleid):
    if rawsampleid.startswith("TCGA"):
        return rawsampleid[:15]
    return rawsampleid


def getcuratedgenes(genelistfile):
    genelist = set()
    with open(genelistfile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            genelist.add(row[1])
    return genelist




def gethotspots(url, type):
    hotspotsjson = json.load(urllib.urlopen(url))
    hotspots = {}
    for hs in hotspotsjson:
        gene = hs['hugoSymbol']
        start = hs['aminoAcidPosition']['start']
        end = hs['aminoAcidPosition']['end']
        if type is None or hs['type'] == type:
            if gene not in hotspots:
                hotspots[gene] = set()
            for i in range(start, end + 1):
                hotspots[gene].add(i)
    return hotspots

missensesinglehotspots = None
indelsinglehotspots = None
_3dhotspots = None
curatedgenes = getcuratedgenes(os.path.join(os.path.dirname(__file__), 'data/curated_genes.txt'))

def inithotspots():
    global missensesinglehotspots
    global indelsinglehotspots
    global _3dhotspots
    global curatedgenes
    missensesinglehotspots = gethotspots(cancerhotspotsbaseurl+"/api/hotspots/single", "single residue")
    indelsinglehotspots = gethotspots(cancerhotspotsbaseurl+"/api/hotspots/single", "in-frame indel")
    _3dhotspots = gethotspots(_3dhotspotsbaseurl+"/api/hotspots/3d", None)
    curatedgenes |= set(missensesinglehotspots.keys())
    curatedgenes |= set(indelsinglehotspots.keys())
    curatedgenes |= set(_3dhotspots.keys())


def processalterationevents(eventfile, outfile, previousoutfile, defaultCancerType, cancerTypeMap, retainonlycuratedgenes):
    inithotspots()
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+', 1000)
    with open(eventfile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)
        outf.write(headers['^-$'])

        outf.write("\tis-a-hotspot")
        outf.write("\tis-a-3d-hotspot")

        # outf.write("\tmutation_effect")
        outf.write("\toncogenic")

        for l in levels:
            outf.write('\t' + l)

        outf.write("\tHighest_level")


        outf.write("\n")

        ihugo = geIndexOfHeader(headers, ['HUGO_SYMBOL', 'HUGO_GENE_SYMBOL'])
        iconsequence = geIndexOfHeader(headers, ['VARIANT_CLASSIFICATION', 'MUTATION_TYPE'])
        ihgvs = geIndexOfHeader(headers, ['ALTERATION', 'HGVSP_SHORT', 'AMINO_ACID_CHANGE', 'FUSION'])
        isample = geIndexOfHeader(headers, ['SAMPLE_ID', 'TUMOR_SAMPLE_BARCODE'])
        istart = geIndexOfHeader(headers, ['PROTEIN_START'])
        iend = geIndexOfHeader(headers, ['PROTEIN_END'])
        iproteinpos = geIndexOfHeader(headers, ['PROTEIN_POSITION'])
        icancertype = geIndexOfHeader(headers, ['ONCOTREE_CODE', 'CANCER_TYPE'])

        posp = re.compile('[0-9]+')

        i = 0
        for row in reader:
            i = i + 1
            if i % 100 == 0:
                print i

            sample = getsampleid(row[isample])
            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            hugo = row[ihugo]
            if retainonlycuratedgenes and hugo not in curatedgenes:
                continue

            consequence = None
            if iconsequence >= 0 and row[iconsequence] != 'NULL':
                consequence = row[iconsequence]
            if consequence in mutationtypeconsequencemap:
                consequence = '%2B'.join(mutationtypeconsequencemap[consequence])

            hgvs = row[ihgvs]
            if hgvs.startswith('p.'):
                hgvs = hgvs[2:]
            if hugo=='TERT' and (row[iconsequence]=='5\'Flank' or row[iconsequence]=='5\'UTR'):
                hgvs = "Promoter Mutations"

            cancertype = defaultCancerType
            if icancertype >= 0:
                cancertype = row[icancertype]
            if sample in cancerTypeMap:
                cancertype = cancerTypeMap[sample]
            if cancertype == "":
                print "Cancer type for all samples must be defined\n"
                print "line "
                print i
                print ": "
                print row
                # continue

            start = None
            if istart >= 0 and row[istart] != 'NULL':
                start = row[istart]

            end = None
            if iend >= 0 and row[iend] != 'NULL':
                end = row[iend]

            if start is None and iproteinpos >= 0 and row[iproteinpos] != "" and row[iproteinpos] != "." and row[iproteinpos] != "-":
                poss = row[iproteinpos].split('/')[0].split('-')
                try:
                    if len(poss) > 0:
                        start = int(poss[0])
                    if len(poss) == 2:
                        end = int(poss[1])
                except ValueError:
                    print "position wrong at line" + str(i) + ": " + row[iproteinpos]

            if start is None and consequence == "missense_variant":
                m = posp.search(hgvs)
                if m:
                    start = m.group()

            if start is not None and end is None:
                end = start

            hotspot = pullsinglehotspots(hugo, hgvs, None, consequence, start, end, cancertype)
            row.append(hotspot)

            _3dhotspot = pull3dhotspots(hugo, hgvs, None, consequence, start, end, cancertype)
            row.append(_3dhotspot)

            oncokblevels = pulloncokb(hugo, hgvs, None, consequence, start, end, cancertype)
            row.append(oncokblevels)

            outf.write('\t'.join(row) + "\n")

    outf.close()


def processsv(svdata, outfile, previousoutfile, defaultCancerType, cancerTypeMap, retainonlycuratedgenes):
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+')
    with open(svdata, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)
        outf.write(headers['^-$'])
        # outf.write("\tmutation_effect")
        outf.write("\toncogenic")
        for l in levels:
            outf.write('\t' + l)
        outf.write("\tHighest_level\n")

        unknownvaraint = pulloncokb("NOT-A-GENE", "Deletion", None, None, None, None, "")

        igene1 = geIndexOfHeader(headers, ['GENE1'])
        igene2 = geIndexOfHeader(headers, ['GENE2'])
        ifusion = geIndexOfHeader(headers, ['FUSION'])
        isample = geIndexOfHeader(headers, ['SAMPLE_ID', 'TUMOR_SAMPLE_BARCODE'])
        icancertype = geIndexOfHeader(headers, ['ONCOTREE_CODE', 'CANCER_TYPE'])

        i = 0
        for row in reader:
            i = i + 1
            if i % 100 == 0:
                print i

            sample = getsampleid(row[isample])

            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            gene1 = None
            gene2 = None
            if igene1>=0:
                gene1 = row[igene1]
            if igene2>=0:
                gene2 = row[igene2]
            if igene1<0 and igene2<0 and ifusion>=0:
                fusion = row[ifusion]
                if fusion.endswith(" fusion"):
                    fusion = fusion[0:len(fusion)-len(" fusion")]
                parts = fusion.split("-")
                gene1 = parts[0]
                gene2 = gene1
                if len(parts)>1 and parts[0]!="intragenic":
                    gene2 = parts[1]

            if retainonlycuratedgenes and gene1 not in curatedgenes and gene2 not in curatedgenes:
                continue

            cancertype = defaultCancerType
            if icancertype >= 0:
                cancertype = row[icancertype]
            if sample in cancerTypeMap:
                cancertype = cancerTypeMap[sample]
            if cancertype == "":
                print "Cancer type for all samples must be defined\n"
                print "line "
                print i
                print ": "
                print row
                # continue

            oncokblevels = None
            if gene1 == gene2:
                oncokblevels = pulloncokb(gene1, "Deletion", None, None, None, None, cancertype)
            else:
                oncokblevels = pulloncokb(gene1+'-'+gene2, '', 'structural_variant', 'fusion', None, None, cancertype)

                # oncokblevels = pulloncokb(gene1, gene1 + '-' + gene2 + " fusion", None, None, None, None, cancertype)
                # if oncokblevels == unknownvaraint:
                #     oncokblevels = pulloncokb(gene1, gene2 + '-' + gene1 + " fusion", None, None, None, None,
                #                               cancertype)
                # if oncokblevels == unknownvaraint:
                #     oncokblevels = pulloncokb(gene2, gene1 + '-' + gene2 + " fusion", None, None, None, None,
                #                               cancertype)
                # if oncokblevels == unknownvaraint:
                #     oncokblevels = pulloncokb(gene2, gene2 + '-' + gene1 + " fusion", None, None, None, None,
                #                               cancertype)
            row.append(oncokblevels)
            outf.write('\t'.join(row) + "\n")

    outf.close()


cnaEventMap = {
    "-2": 'Deletion',
    "-1.5": 'Deletion',
    "-1": 'Loss',
    "0": 'Diploid',
    "1": 'Gain',
    "2": 'Amplification'
}


def processcnagisticdata(cnafile, outfile, previousoutfile, defaultCancerType, cancerTypeMap, retainonlycuratedgenes):
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+', 1000)
    with open(cnafile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')
        headers = readheaders(reader)
        startofsamples = getfirstcolumnofsampleingisticdata(headers['^-$'].split('\t'))
        rawsamples = headers['^-$'].split('\t')[startofsamples:]
        samples = []
        for rs in rawsamples:
            samples.append(getsampleid(rs))

        if defaultCancerType == '' and not set(cancerTypeMap.keys()).issuperset(set(samples)):
            print "Cancer type for all samples must be defined\n"
            print "samples with cancer type:\n"
            print cancerTypeMap.keys()
            print "\nsamples in cna file:\n"
            print samples
            # quit()

        outf.write('SAMPLE_ID\tCANCER_TYPE\tHUGO_SYMBOL\tALTERATION')
        # outf.write("\tmutation_effect")
        outf.write("\toncogenic")
        for l in levels:
            outf.write('\t' + l)
        outf.write("\tHighest_level\n")

        i = 0
        for row in reader:
            i = i + 1
            if i % 100 == 0:
                print i

            hugo = row[0]
            if retainonlycuratedgenes and hugo not in curatedgenes:
                continue

            for rawsample in rawsamples:
                rawsample = rawsample.upper()
                if rawsample in headers:
                    cna = row[headers[rawsample]]
                    if cna in cnaEventMap:
                        alteration = cnaEventMap[cna]
                        if alteration == "Amplification" or alteration == "Deletion":
                            cancertype = defaultCancerType
                            sample = getsampleid(rawsample)

                            if sampleidsfilter and sample not in sampleidsfilter:
                                continue

                            if sample in cancerTypeMap:
                                cancertype = cancerTypeMap[sample]
                            oncokblevels = pulloncokb(hugo, alteration, None, None, None, None, cancertype)
                            outf.write(sample + "\t")
                            outf.write(cancertype + "\t")
                            outf.write(hugo + "\t")
                            outf.write(alteration + "\t")
                            outf.write(oncokblevels)
                            outf.write('\n')
    outf.close()


def getfirstcolumnofsampleingisticdata(headers):
    header0 = headers[0].lower()
    if header0 != "hugo_symbol" and header0 != "gene symbol":
        print "Gistic data should start with Hugo_Symbol"
        quit()

    header1 = headers[1].lower()
    if header1 != "entrez_gene_id" and header1 != "locus id":
        return 1

    header2 = headers[2].lower()
    if header2 != "cytoband":
        return 2

    return 3


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def processclinicaldata(annotatedmutfiles, clinicalfile, outfile):
    samplelevels = {}
    sampleleveltreatments = {}
    sampledrivers = {}
    sampleactionablecount = {}
    for annotatedmutfile in annotatedmutfiles:
        # print annotatedmutfile
        with open(annotatedmutfile, 'rU') as mutfile:
            reader = csv.reader(mutfile, delimiter='\t')
            headers = readheaders(reader)

            igene1 = geIndexOfHeader(headers, ['HUGO_SYMBOL', 'GENE1', 'HUGO_GENE_SYMBOL'])  # fusion
            igene2 = geIndexOfHeader(headers, ['HUGO_SYMBOL', 'GENE2', 'HUGO_GENE_SYMBOL'])  # fusion

            ihugo = geIndexOfHeader(headers, ['HUGO_SYMBOL', 'HUGO_GENE_SYMBOL'])
            iconsequence = geIndexOfHeader(headers, ['VARIANT_CLASSIFICATION', 'MUTATION_TYPE'])
            ihgvs = geIndexOfHeader(headers, ['ALTERATION', 'HGVSP_SHORT', 'AMINO_ACID_CHANGE', 'AA_CHANGE', 'FUSION'])
            isample = geIndexOfHeader(headers, ['SAMPLE_ID', 'TUMOR_SAMPLE_BARCODE'])
            istart = geIndexOfHeader(headers, ['PROTEIN_START'])
            iend = geIndexOfHeader(headers, ['PROTEIN_END'])
            icancertype = geIndexOfHeader(headers, ['ONCOTREE_CODE', 'CANCER_TYPE'])
            # imutationeffect = headers['MUTATION_EFFECT']
            ioncogenic = headers['ONCOGENIC']

            isfusion = igene1 != -1 & igene2 != -1
            ismutorcna = ihugo != -1 & ihgvs != -1

            if not isfusion and not ismutorcna:
                print "missing proper header"
                exit()

            for row in reader:
                sample = getsampleid(row[isample])

                oncogenic = ""
                if ioncogenic < len(row):
                    oncogenic = row[ioncogenic].lower()
                if sample not in samplelevels:
                    samplelevels[sample] = {}
                    sampleleveltreatments[sample] = {}
                    sampledrivers[sample] = []
                    sampleactionablecount[sample] = 0

                hugo = row[ihugo]
                alteration = row[ihgvs]
                gene1 = row[igene1]
                gene2 = row[igene2]

                variant = "NA"
                if ismutorcna:
                    variant = hugo + " " + alteration
                elif isfusion:
                    if gene1 == gene2:
                        variant = gene1 + " intragenic deletion"
                    else:
                        variant = gene1 + "-" + gene2 + " fusion"

                if oncogenic == "oncogenic" or oncogenic == "likely oncogenic" or oncogenic == "predicted oncogenic":
                    sampledrivers[sample].append(variant)

                for l in levels:
                    il = headers[l]
                    if il < len(row) and row[il] != '':
                        if l not in samplelevels[sample]:
                            samplelevels[sample][l] = []
                            sampleleveltreatments[sample][l] = []
                        samplelevels[sample][l].append(row[il] + "(" + variant + ")")
                        sampleleveltreatments[sample][l].extend(row[il].split(","))

                        if not l.startswith('LEVEL_R'):
                            sampleactionablecount[sample] += 1

    outf = open(outfile, 'w+')

    # export to anntoated file
    with open(clinicalfile, 'rU') as clinfile:
        reader = csv.reader(clinfile, delimiter='\t')
        headers = readheaders(reader)
        outf.write(headers['^-$'])
        for l in levels:
            outf.write('\t' + l)
        outf.write('\tHIGHEST_LEVEL\toncogenic_mutations\t#actionable_mutations\t#oncogenic_mutations\n')
        isample = headers['SAMPLE_ID']

        for row in reader:
            sample = row[isample]

            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            # print row
            outf.write('\t'.join(row))

            for l in levels:
                outf.write('\t')
                if sample in samplelevels and l in samplelevels[sample]:
                    outf.write(";".join(samplelevels[sample][l]))

            highestlevel = ''
            if sample in sampleleveltreatments:
                highestlevel = gethighestsensitivitylevel(sampleleveltreatments[sample])
            # if highestlevel == '':
            #     if sample in sampledrivers and len(sampledrivers[sample])>0:
            #         highestlevel = 'Oncogenic, no level'
            #     else:
            #         highestlevel = "VUS"
            outf.write('\t' + highestlevel)


            actionablecount = 0
            if sample in sampleactionablecount:
                actionablecount = sampleactionablecount[sample]

            drivercount = 0
            drivermutations = ""
            if sample in sampledrivers:
                drivercount = len(sampledrivers[sample])
                drivermutations = ";".join(sampledrivers[sample])

            outf.write('\t' + drivermutations)
            outf.write('\t' + str(actionablecount))
            outf.write('\t' + str(drivercount))

            outf.write('\n')

    outf.close()

def plotclinicalactionability(annotatedclinicalfile, outfile, parameters):
    if os.path.isfile(outfile):
        os.remove(outfile)

    extlevels = levels + ["ONCOGENIC", "VUS"]
    if "levels" in parameters:
        extlevels = parameters["levels"]

    with open(annotatedclinicalfile, 'rU') as clinfile:
        reader = csv.reader(clinfile, delimiter='\t')
        headers = readheaders(reader)
        isample = headers['SAMPLE_ID']
        ilevel = headers['HIGHEST_LEVEL']
        ioncogenic = headers['ONCOGENIC_MUTATIONS']
        icat = headers[parameters["catogerycolumn"].upper()] #e.g. "CANCER_TYPE"

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

    # level colors
    levelcolors = {
        'LEVEL_1': '#33A02C',
        'LEVEL_2A': '#1F78B4',
        'LEVEL_2B': '#80B1D3',
        'LEVEL_3A': '#984EA3',
        'LEVEL_3B': '#BE98CE',
        'LEVEL_4': '#a8a8a8',
        'LEVEL_R1': '#EE3424',
        'LEVEL_R2': '#F79A92',
        'LEVEL_R3': '#FCD6D3',
        'ONCOGENIC': '#ffdab9',
        'VUS': '#d1d1d1',
        'Other': 'grey'
    }

    # level legend
    levellegend = {
        'LEVEL_1': 'Level 1',
        'LEVEL_2A': 'Level 2A',
        'LEVEL_2B': 'Level 2B',
        'LEVEL_3A': 'Level 3A',
        'LEVEL_3B': 'Level 3B',
        'LEVEL_4': 'Level 4',
        'LEVEL_R1': 'Level R1',
        'LEVEL_R2': 'Level R2',
        'LEVEL_R3': 'Level R3',
        'ONCOGENIC': 'Oncogenic, no level',
        'VUS': 'VUS',
        'Other': 'Other'
    }

    # plot
    catarray = [] # cancer types
    catactionabilityarray = [] # actionabiligy percentages per cancer type
    catoncogenicarray = [] # actionabiligy percentages per cancer type
    for cat in catsamplecount:
        if catsamplecount[cat] >= parameters["thresholdcat"]:
            catarray.append(cat)
            catactionabilityarray.append(catactionablesamplecount[cat] * 100.0 / catsamplecount[cat])
            catoncogenicarray.append(oncogenicsamplecount[cat] * 100.0 / catsamplecount[cat])

    ncat = len(catarray)
    if ncat > 0:
        # sort categories (cancer type) based on actionability and then oncogenic frequency
        order = reversed(sorted(range(ncat),key=lambda x:(catactionabilityarray[x],catoncogenicarray[x])))
        catarray = [catarray[i] for i in order]

        ind = range(ncat)

        f = plt.figure()

        legends = []
        plts = []
        accumlevelcancerperc = [0] * ncat
        for level in extlevels:
            if level not in levelcatsamplecount:
                continue

            levelcancerperc = [0] * ncat
            for k in ind:
                cat = catarray[k]
                if catsamplecount[cat] < parameters["thresholdcat"]:
                    continue
                if cat in levelcatsamplecount[level]:
                    levelcancerperc[k] = levelcatsamplecount[level][cat] * 100.0 / catsamplecount[cat]

            width = 0.75
            plts = [plt.bar(ind, levelcancerperc, width, color=levelcolors[level], bottom=accumlevelcancerperc)] + plts
            legends = [levellegend[level]] + legends
            accumlevelcancerperc = map(sum, zip(accumlevelcancerperc,levelcancerperc))

        ax = plt.axes()
        ax.set_axisbelow(True)
        ax.set_aspect(0.1)
        # ax.yaxis.grid(linestyle="dotted", color="lightgray") # horizontal lines
        plt.margins(0.01)
        plt.tick_params(axis='y', which='major', labelsize=6)
        plt.ylabel('% of samples')
        plt.title('OncoKB Actionability')
        plt.xticks([i+0.5 for i in ind], catarray, rotation=60, ha="right", fontsize=6)
        plt.subplots_adjust(left=0.2, bottom=0.3)
        # plt.yticks(np.arange(0, 81, 10))
        plt.legend(plts, legends, fontsize=6, bbox_to_anchor=(1.01, 1), loc="upper left")
        plt.gcf().text(0.90, 0.1, "Generated by OncoKB\n[Chakravarty et al., JCO PO 2017]", fontsize=6, horizontalalignment='right', verticalalignment='bottom')

        # plt.show()
        f.savefig(outfile, bbox_inches='tight')

def processmutationdata(mutfile, outfile, clinicaldata):
    outf = open(outfile, 'w+')
    with open(mutfile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')
        headers = readheaders(reader)

        ihugo = headers['Hugo_Symbol']
        iconsequence = headers['Consequence']
        ihgvs = headers['HGVSp_Short']
        isample = headers['Tumor_Sample_Barcode']
        istart = headers['protein_start']
        iend = headers['protein_end']

        i = 0
        for row in reader:
            if i % 100 == 0:
                print i
            i = i + 1

            sample = row[isample]
            hugo = row[ihugo]
            consequence = row[iconsequence]
            if consequence == 'NULL':
                consequence = None
            if consequence in mutationtypeconsequencemap:
                consequence = '%2B'.join(mutationtypeconsequencemap[consequence])

            hgvs = row[ihgvs]
            if hgvs.startswith('p.'):
                hgvs = hgvs[2:]
            cancertype = None

            start = row[istart]
            if start == 'NULL':
                start = None

            end = row[iend]
            if end == 'NULL':
                end = None

            if sample in clinicaldata:
                cancertype = clinicaldata[sample]
            oncokbevidences = pulloncokb(hugo, hgvs, None, consequence, start, end, cancertype)
            annotatedrow = [hugo, consequence, start, end, hgvs, sample, cancertype, oncokbevidences]
            outf.write('\t'.join(annotatedrow) + "\n")

    outf.close()


oncokbcache = {}


def cacheannotated(annotatedfile, defaultCancerType, cancerTypeMap):
    with open(annotatedfile, 'rU') as infile:
        try:
            reader = csv.reader(infile, delimiter='\t')
            headers = readheaders(reader)

            ihugo = geIndexOfHeader(headers, ['HUGO_SYMBOL', 'HUGO_GENE_SYMBOL'])
            iconsequence = geIndexOfHeader(headers, ['VARIANT_CLASSIFICATION', 'MUTATION_TYPE'])
            ihgvs = geIndexOfHeader(headers, ['ALTERATION', 'HGVSP_SHORT', 'AMINO_ACID_CHANGE', 'AA_CHANGE', 'FUSION'])
            isample = geIndexOfHeader(headers, ['SAMPLE_ID', 'TUMOR_SAMPLE_BARCODE'])
            istart = geIndexOfHeader(headers, ['PROTEIN_START'])
            iend = geIndexOfHeader(headers, ['PROTEIN_END'])
            icancertype = geIndexOfHeader(headers, ['ONCOTREE_CODE', 'CANCER_TYPE'])
            # imutationeffect = headers['MUTATION_EFFECT']
            ioncogenic = headers['ONCOGENIC']

            for row in reader:
                try:
                    hugo = row[ihugo]
                    if hugo not in curatedgenes:
                        continue

                    hgvs = row[ihgvs]
                    if hgvs.startswith('p.'):
                        hgvs = hgvs[2:]

                    sample = row[isample]
                    cancertype = defaultCancerType
                    if icancertype >= 0:
                        cancertype = row[icancertype]
                    if sample in cancerTypeMap:
                        cancertype = cancerTypeMap[sample]
                    key = '-'.join([hugo, hgvs, cancertype])
                    #            oncokb = row[ioncokb]

                    oncokbcache[key] = {}
                    # oncokbcache[key]['mutation_effect'] = row[imutationeffect]
                    oncokbcache[key]['oncogenic'] = row[ioncogenic]
                    for l in levels:
                        il = headers[l]
                        if il < len(row):
                            oncokbcache[key][l] = row[il].split(',')
                        else:
                            oncokbcache[key][l] = []
                except Exception:
                    pass
        except Exception:
            pass

def geIndexOfHeader(headers, keywords):
    for k in keywords:
        if k in headers:
            return headers[k]
    return -1


def pullsinglehotspots(hugo, proteinchange, alterationtype, consequence, start, end, cancertype):
    try:
        if hugo in missensesinglehotspots and consequence == "missense_variant":
            for i in range(int(start), int(end) + 1):
                if i in missensesinglehotspots[hugo]:
                    return "Y"
        if hugo in indelsinglehotspots and (consequence == "inframe_insertion" or consequence == "inframe_deletion"):
            for i in range(int(start), int(end) + 1):
                if i in indelsinglehotspots[hugo]:
                    return "Y"
    except TypeError:
        print hugo + ":" + str(start) + "-" + str(end)
    return ""


def pull3dhotspots(hugo, proteinchange, alterationtype, consequence, start, end, cancertype):
    try:
        if hugo in _3dhotspots and consequence == "missense_variant":
            for i in range(int(start), int(end) + 1):
                if i in _3dhotspots[hugo]:
                    return "Y"
    except TypeError:
        print hugo + ":" + str(start) + "-" + str(end)
    return ""


def pulloncokb(hugo, proteinchange, alterationtype, consequence, start, end, cancertype):
    if hugo not in curatedgenes and alterationtype and alterationtype.lower() != 'fusion' and alterationtype.lower() != 'structural_variant':
        return ""

    key = '-'.join([hugo, proteinchange, cancertype])
    if key not in oncokbcache:
        # url = 'http://dashi-dev.cbio.mskcc.org:8080/oncokb/api/indicator.json?source=cbioportal'
        # url = 'http://localhost:8080/oncokb/api/indicator.json?source=cbioportal'
        # url = 'http://dashi.cbio.mskcc.org:38080/internal/legacy-api/indicator.json?source=cbioportal'
        url = oncokbbaseurl+'/legacy-api/indicator.json?source=cbioportal'
        url += '&hugoSymbol=' + hugo
        url += '&alteration=' + proteinchange
        url += '&tumorType=' + cancertype
        if consequence:
            url += '&consequence=' + consequence
        if start and start != '\\N' and start != 'NULL':
            url += '&proteinStart=' + str(start)
        if end and end != '\\N' and end != 'NULL':
            url += '&proteinEnd=' + str(end)
        if alterationtype and alterationtype != 'NULL' and alterationtype != '':
            url += '&alterationType=' + alterationtype

        oncokbdata = {}
        for l in levels:
            oncokbdata[l] = []

        oncokbdata['oncogenic'] = "Unknown"
        # oncokbdata['mutation_effect'] = "Unknown"

        try:
            evidences = json.load(urllib.urlopen(url))
            # if not evidences['geneExist'] or (not evidences['variantExist'] and not evidences['alleleExist']):
            #     return ''

            # # mutation effect
            # for e in evidences[0]:
            #     evidenceType = e['evidenceType']
            #     if evidenceType == "MUTATION_EFFECT" and e['knownEffect'] is not None:
            #         oncokbdata['mutation_effect'] = e['knownEffect']
            #         break

            # oncogenic
            oncokbdata['oncogenic'] = evidences['oncogenic']

            # get treatment
            for treatment in evidences['treatments']:
                level = treatment['level']

                if level not in levels:
                    oncokbdata[level].append('')
                else:
                    drugs = treatment['drugs']
                    if len(drugs) == 0:
                        oncokbdata[level].append('[NOT SPECIFIED]')
                    else:
                        drugnames = []
                        for drug in drugs:
                            drugnames.append(drug['drugName'])
                        oncokbdata[level].append('+'.join(drugnames))
        except:
            print url
            # sys.exit()

        oncokbcache[key] = oncokbdata

    oncokbdata = oncokbcache[key]
    ret = []
    # ret.append(oncokbdata['mutation_effect'])
    ret.append(oncokbdata['oncogenic'])
    for l in levels:
        ret.append(','.join(oncokbdata[l]))
    ret.append(gethighestsensitivitylevel(oncokbdata))

    ret = "\t".join(ret)
    return ret

def gethighestsensitivitylevel(oncokbdata):
    r1 = set()
    if "LEVEL_R1" in oncokbdata:
        r1 = set(oncokbdata["LEVEL_R1"])
    for l in levels:
        if l.startswith("LEVEL_R") or l not in oncokbdata or oncokbdata[l]=='':
            continue
        if not r1.issuperset(set(oncokbdata[l])):
            return l
    return ""


def gettreatments(evidence):
    treatments = []
    for t in evidence['treatments']:
        drugs = []
        for d in t['drugs']:
            drugs.append(d['drugName'])
        treatments.append('+'.join(drugs))
    return treatments


def readCancerTypes(clinicalFile, data):
    with open(clinicalFile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')
        headers = readheaders(reader)

        iSample = geIndexOfHeader(headers, ['SAMPLE_ID'])
        iCancerType = geIndexOfHeader(headers, ['ONCOTREE_CODE', 'CANCER_TYPE'])

        for row in reader:
            data[row[iSample]] = row[iCancerType]

    return data


def readheaders(reader):
    headers = {}
    for row in reader:
        if not row[0].startswith("#"):
            headers["^-$"] = '\t'.join(row)  # the whole line
            i = 0
            for h in row:
                headers[h.upper()] = i
                i = i + 1
            break
    return headers
