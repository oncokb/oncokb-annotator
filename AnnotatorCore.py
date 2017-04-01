#!/usr/bin/python

import sys
import getopt
import csv
import json
import urllib
import os.path
import re

csv.field_size_limit(sys.maxsize) # for reading large files

baseurl = "http://oncokb.org"

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
    with open(genelistfile, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            genelist.add(row[1])
    return genelist


curatedgenes = getcuratedgenes('data/curated_genes.txt')


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


missensesinglehotspots = gethotspots("http://cancerhotspots.mskcc.org/internal/api/hotspots/single", "single residue")
indelsinglehotspots = gethotspots("http://cancerhotspots.mskcc.org/internal/api/hotspots/single", "in-frame indel")
_3dhotspots = gethotspots("http://3dhotspots.org/3d/api/hotspots/3d", None)
curatedgenes |= set(missensesinglehotspots.keys())
curatedgenes |= set(indelsinglehotspots.keys())
curatedgenes |= set(_3dhotspots.keys())


def processalterationevents(eventfile, outfile, defaultCancerType, cancerTypeMap, retainonlycuratedgenes):
    if os.path.isfile(outfile):
        cacheannotated(outfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+', 1000)
    with open(eventfile, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)
        outf.write(headers['^-$'])
        # outf.write("\tmutation_effect")
        outf.write("\toncogenic")
        for l in levels:
            outf.write('\t' + l)
        outf.write("\tHighest_level")
        outf.write("is-a-hotspot")
        outf.write("is-a-3d-hotspot")
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

            if start is None and iproteinpos >= 0 and row[iproteinpos] != "" and row[iproteinpos] != ".":
                poss = row[iproteinpos].split('/')[0].split('-')
                try:
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

            oncokblevels = pulloncokb(hugo, hgvs, None, consequence, start, end, cancertype)
            row.append(oncokblevels)

            hotspot = pullsinglehotspots(hugo, hgvs, None, consequence, start, end, cancertype)
            row.append(hotspot)

            _3dhotspot = pull3dhotspots(hugo, hgvs, None, consequence, start, end, cancertype)
            row.append(_3dhotspot)

            outf.write('\t'.join(row) + "\n")

    outf.close()


def processsv(svdata, outfile, defaultCancerType, cancerTypeMap, retainonlycuratedgenes):
    if os.path.isfile(outfile):
        cacheannotated(outfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+')
    with open(svdata, 'r') as infile:
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
        isample = geIndexOfHeader(headers, ['SAMPLE_ID', 'TUMOR_SAMPLE_BARCODE'])
        icancertype = geIndexOfHeader(headers, ['ONCOTREE_CODE', 'CANCER_TYPE'])

        i = 0
        for row in reader:
            i = i + 1
            if i % 100 == 0:
                print i

            sample = getsampleid(row[isample])
            gene1 = row[igene1]
            gene2 = row[igene2]
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
                # oncokblevels = pulloncokb(gene1+'-'+gene2, '', 'fusion', None, None, None, cancertype)

                oncokblevels = pulloncokb(gene1, gene1 + '-' + gene2 + " fusion", None, None, None, None, cancertype)
                if oncokblevels == unknownvaraint:
                    oncokblevels = pulloncokb(gene1, gene2 + '-' + gene1 + " fusion", None, None, None, None,
                                              cancertype)
                if oncokblevels == unknownvaraint:
                    oncokblevels = pulloncokb(gene2, gene1 + '-' + gene2 + " fusion", None, None, None, None,
                                              cancertype)
                if oncokblevels == unknownvaraint:
                    oncokblevels = pulloncokb(gene2, gene2 + '-' + gene1 + " fusion", None, None, None, None,
                                              cancertype)
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


def processcnagisticdata(cnafile, outfile, defaultCancerType, cancerTypeMap, retainonlycuratedgenes):
    if os.path.isfile(outfile):
        cacheannotated(outfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+', 1000)
    with open(cnafile, 'r') as infile:
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
                if rawsample in headers:
                    cna = row[headers[rawsample]]
                    if cna in cnaEventMap:
                        alteration = cnaEventMap[cna]
                        if alteration == "Amplification" or alteration == "Deletion":
                            cancertype = defaultCancerType
                            sample = getsampleid(rawsample)
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

                if oncogenic == "oncogenic" or oncogenic == "likely oncogenic":
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

    with open(clinicalfile, 'r') as clinfile:
        reader = csv.reader(clinfile, delimiter='\t')
        headers = readheaders(reader)
        outf.write(headers['^-$'])
        for l in levels:
            outf.write('\t' + l)
        outf.write('\tHIGHEST_LEVEL\toncogenic_mutations\t#actionable_mutations\t#oncogenic_mutations\n')
        isample = headers['SAMPLE_ID']

        for row in reader:
            # print row
            outf.write('\t'.join(row))
            sample = row[isample]

            for l in levels:
                outf.write('\t')
                if sample in samplelevels and l in samplelevels[sample]:
                    outf.write(";".join(samplelevels[sample][l]))

            highestlevel = ''
            if sample in sampleleveltreatments:
                highestlevel = gethighestsensitivitylevel(sampleleveltreatments[sample])
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


def processmutationdata(mutfile, outfile, clinicaldata):
    outf = open(outfile, 'w+')
    with open(mutfile, 'r') as infile:
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
    with open(annotatedfile, 'r') as infile:
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
        if hugo in indelsinglehotspots and (consequence == "inframe_insertion" or consequence == "inframe_insertion"):
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
    if hugo not in curatedgenes and alterationtype and alterationtype.lower() != 'fusion':
        return ""

    key = '-'.join([hugo, proteinchange, cancertype])
    if key not in oncokbcache:
        # url = 'http://dashi-dev.cbio.mskcc.org:8080/oncokb/api/indicator.json?source=cbioportal'
        # url = 'http://localhost:8080/oncokb/api/indicator.json?source=cbioportal'
        # url = 'http://dashi.cbio.mskcc.org:38080/internal/legacy-api/indicator.json?source=cbioportal'
        url = baseurl+'/legacy-api/indicator.json?source=cbioportal'
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
            if not evidences['geneExist'] or (not evidences['variantExist'] and not evidences['alleleExist']):
                return ''

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
                    continue
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
        if l.startswith("LEVEL_R") or l not in oncokbdata:
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
    with open(clinicalFile, 'r') as infile:
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
