#!/usr/bin/python
import json
import sys
import csv
import requests
import os.path
import logging
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import date
import ctypes as ct

logging.basicConfig(level=logging.INFO)
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

log = logging.getLogger('AnnotatorCore')

csv.field_size_limit(int(ct.c_ulong(-1).value // 2)) # Deal with overflow problem on Windows, https://stackoverflow.com/questions/15063936/csv-error-field-larger-than-field-limit-131072
sizeLimit = csv.field_size_limit()
csv.field_size_limit(sizeLimit) # for reading large files

oncokbapiurl = "https://www.oncokb.org/api/v1"
oncokbapibearertoken = ""

def setoncokbbaseurl(u):
    global oncokbapiurl
    oncokbapiurl = u.rstrip('/') + '/api/v1'

def setoncokbapitoken(t):
    global oncokbapibearertoken
    oncokbapibearertoken = t.strip()

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
    log.info(len(sampleidsfilter))


GENE_IN_ONCOKB_HEADER = 'GENE_IN_ONCOKB'
VARIANT_IN_ONCOKB_HEADER = 'VARIANT_IN_ONCOKB'

GENE_IN_ONCOKB_DEFAULT = 'False'
VARIANT_IN_ONCOKB_DEFAULT = 'False'

levels = [
    'LEVEL_1',
    'LEVEL_2',
    'LEVEL_3A',
    'LEVEL_3B',
    'LEVEL_4',
    'LEVEL_R1',
    'LEVEL_R2',
    'LEVEL_R3'
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

# column headers
HUGO_HEADERS = ['HUGO_SYMBOL', 'HUGO_GENE_SYMBOL', 'GENE']
CONSEQUENCE_HEADERS = ['VARIANT_CLASSIFICATION', 'MUTATION_TYPE']
HGVS_HEADERS = ['ALTERATION', 'HGVSP_SHORT', 'HGVSP', 'AMINO_ACID_CHANGE', 'FUSION']
SAMPLE_HEADERS = ['SAMPLE_ID', 'TUMOR_SAMPLE_BARCODE']
START_HEADERS = ['PROTEIN_START']
END_HEADERS = ['PROTEIN_END']
PROTEIN_POSITION_HEADERS = ['PROTEIN_POSITION']
CANCER_TYPE_HEADERS = ['ONCOTREE_CODE', 'CANCER_TYPE']
FUSION_HEADERS = ['FUSION']

POST_QUERIES_THRESHOLD = 1000

def getsampleid(rawsampleid):
    if rawsampleid.startswith("TCGA"):
        return rawsampleid[:15]
    return rawsampleid


def getOncokbInfo():
    ret = ['Files annotated on ' + date.today().strftime('%m/%d/%Y') + "\nOncoKB API URL: "+oncokbapiurl]
    try:
        info = requests.get(oncokbapiurl + "/info").json()
        ret.append('\nOncoKB data version: ' + info['dataVersion']['version']+', released on ' + info['dataVersion']['date'])
    except:
        log.error("error when fetch OncoKB info")
    return ''.join(ret)


def generateReadme(outfile):
    outf = open(outfile, 'w+', 1000)
    outf.write(getOncokbInfo())
    outf.close()

def gethotspots(url, type):
    hotspots = {}
    response = requests.get(url)
    if response.status_code == 200:
        hotspotsjson = response.json()

        for hs in hotspotsjson:
            gene = hs['hugoSymbol']
            start = hs['aminoAcidPosition']['start']
            end = hs['aminoAcidPosition']['end']
            if type is None or hs['type'] == type:
                if gene not in hotspots:
                    hotspots[gene] = set()
                for i in range(start, end + 1):
                    hotspots[gene].add(i)
    else:
        log.error("error when processing %s \n" % url +
                  "reason: %s" % response.reason)
    return hotspots


def makeoncokbpostrequest(url, body):
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer %s' % oncokbapibearertoken
    }
    return requests.post(url, headers=headers, data=json.dumps(body, default=lambda o: o.__dict__))


def makeoncokbgetrequest(url):
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer %s' % oncokbapibearertoken
    }
    return requests.get(url, headers=headers)


def getcuratedgenes():
    global curatedgenes
    url = oncokbapiurl + "/utils/allCuratedGenes.json"
    response = makeoncokbgetrequest(url)
    if response.status_code == 200:
        curatedgenesjson = response.json()

        for curatedgene in curatedgenesjson:
            if curatedgene['hugoSymbol'] is not None:
                curatedgenes.append(curatedgene['hugoSymbol'])
    else:
        log.error("error when processing %s \n" % url +
                  "reason: %s" % response.reason)

missensesinglehotspots = None
indelsinglehotspots = None
_3dhotspots = None
curatedgenes = []

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


conversiondict = {'Ala': 'A',
                  'Asx': 'B',
                  'Cys': 'C',
                  'Asp': 'D',
                  'Glu': 'E',
                  'Phe': 'F',
                  'Gly': 'G',
                  'His': 'H',
                  'Ile': 'I',
                  'Lys': 'K',
                  'Leu': 'L',
                  'Met': 'M',
                  'Asn': 'N',
                  'Pro': 'P',
                  'Gln': 'Q',
                  'Arg': 'R',
                  'Ser': 'S',
                  'Thr': 'T',
                  'Val': 'V',
                  'Trp': 'W',
                  'Tyr': 'Y',
                  'Glx': 'Z'
                  }
conversionlist = conversiondict.keys()
def conversion(hgvs):
    threecharactersearch = re.findall('[a-zA-Z]{3}\d+', hgvs, flags=re.IGNORECASE)
    if threecharactersearch:
        if any(letters.lower() in hgvs.lower() for letters in conversionlist):
            return replace_all(hgvs)
    return hgvs

def replace_all(hgvs):
    # Author: Thomas Glaessle
    pattern = re.compile('|'.join(conversionlist), re.IGNORECASE)
    return pattern.sub(lambda m: conversiondict[m.group().capitalize()], hgvs)


def append_annotation_to_file(outf, rows, annotations):
    if len(rows) != len(annotations):
        log.error('The length of the rows and annotations do not match')

    for index, annotation in enumerate(annotations):
        row = rows[index]
        if annotation is not None:
            row.append(annotation)
        outf.write('\t'.join(row) + "\n")

def processalterationevents(eventfile, outfile, previousoutfile, defaultCancerType, cancerTypeMap,
                            retainonlycuratedgenes, annotatehotspots):
    if annotatehotspots:
        inithotspots()
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+', 1000)
    with open(eventfile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)

        ncols = headers["length"]

        outf.write(headers['^-$'])

        if annotatehotspots:
            outf.write("\tis-a-hotspot")
            outf.write("\tis-a-3d-hotspot")

        outf.write("\t" + GENE_IN_ONCOKB_HEADER)
        outf.write("\t" + VARIANT_IN_ONCOKB_HEADER)

        outf.write("\tMUTATION_EFFECT")
        outf.write("\tONCOGENIC")

        for l in levels:
            outf.write('\t' + l)

        outf.write("\tHIGHEST_LEVEL")

        outf.write("\tCITATIONS")

        outf.write("\n")

        ihugo = geIndexOfHeader(headers, HUGO_HEADERS)
        iconsequence = geIndexOfHeader(headers, CONSEQUENCE_HEADERS)
        ihgvs = geIndexOfHeader(headers, HGVS_HEADERS)
        isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
        istart = geIndexOfHeader(headers, START_HEADERS)
        iend = geIndexOfHeader(headers, END_HEADERS)
        iproteinpos = geIndexOfHeader(headers, PROTEIN_POSITION_HEADERS)
        icancertype = geIndexOfHeader(headers, CANCER_TYPE_HEADERS)

        posp = re.compile('[0-9]+')

        i = 0
        queries = []
        rows = []
        for row in reader:
            i = i + 1

            if i % POST_QUERIES_THRESHOLD == 0:
                log.info(i)

            row = padrow(row, ncols)

            sample = getsampleid(row[isample])
            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            hugo = row[ihugo]

            consequence = None
            if iconsequence >= 0 and row[iconsequence] != 'NULL':
                consequence = row[iconsequence]
            if consequence in mutationtypeconsequencemap:
                consequence = '%2B'.join(mutationtypeconsequencemap[consequence])

            hgvs = row[ihgvs]
            if hgvs.startswith('p.'):
                hgvs = hgvs[2:]
            if hugo=='TERT' and (row[iconsequence]=='5\'Flank' or row[iconsequence]=='5\'UTR'):
                hgvs = "Promoter Mutation"

            cancertype = defaultCancerType
            if icancertype >= 0:
                cancertype = row[icancertype]
            if sample in cancerTypeMap:
                cancertype = cancerTypeMap[sample]
            if cancertype == "":
                log.info("Cancer type for the sample should be defined for a more accurate result\nline %s: %s\n" % (i, row))
                # continue

            hgvs = conversion(hgvs)

            start = None
            if istart >= 0 and row[istart] != 'NULL' and row[istart] != '':
                start = row[istart]

            end = None
            if iend >= 0 and row[iend] != 'NULL' and row[iend] != '':
                end = row[iend]

            if start is None and iproteinpos >= 0 and row[iproteinpos] != "" and row[iproteinpos] != "." and row[iproteinpos] != "-":
                poss = row[iproteinpos].split('/')[0].split('-')
                try:
                    if len(poss) > 0:
                        start = int(poss[0])
                    if len(poss) == 2:
                        end = int(poss[1])
                except ValueError:
                    log.info("position wrong at line %s: %s" % (str(i), row[iproteinpos]))

            if start is None and consequence == "missense_variant":
                m = posp.search(hgvs)
                if m:
                    start = m.group()

            if start is not None and end is None:
                end = start

            if annotatehotspots:
                hotspot = pullsinglehotspots(hugo, hgvs, None, consequence, start, end, cancertype)
                row.append(hotspot)

                _3dhotspot = pull3dhotspots(hugo, hgvs, None, consequence, start, end, cancertype)
                row.append(_3dhotspot)

            if not retainonlycuratedgenes or hugo in curatedgenes:
                query = ProteinChangeQuery(hugo, hgvs, cancertype, consequence, start, end)
                queries.append(query)
                rows.append(row)
            else:
                # Include Gene in OncoKB and Variant in OncoKB
                row.append(GENE_IN_ONCOKB_DEFAULT + '\t' + VARIANT_IN_ONCOKB_DEFAULT)
                outf.write('\t'.join(row) + "\n")

            if len(queries) == POST_QUERIES_THRESHOLD:
                annotations = pull_mutation_info(queries)
                append_annotation_to_file(outf, rows, annotations)
                queries = []
                rows = []

        if len(queries) > 0:
            annotations = pull_mutation_info(queries)
            append_annotation_to_file(outf, rows, annotations)


    outf.close()


def getgenesfromfusion(fusion, nameregex=None):
    GENES_REGEX = "([A-Za-z\d]+-[A-Za-z\d]+)" if nameregex is None else nameregex
    searchresult = re.search(GENES_REGEX, fusion, flags=re.IGNORECASE)
    gene1=None
    gene2=None
    if searchresult:
        parts = searchresult.group(1).split("-")
        gene1 = parts[0]
        gene2 = gene1
        if len(parts) > 1 and parts[1] != "intragenic":
            gene2 = parts[1]
    else:
        gene1=gene2=fusion
    return gene1, gene2

def processsv(svdata, outfile, previousoutfile, defaultCancerType, cancerTypeMap, retainonlycuratedgenes, nameregex):
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+')
    with open(svdata, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)

        ncols = headers["length"]

        outf.write(headers['^-$'])
        outf.write("\t" + GENE_IN_ONCOKB_HEADER)
        outf.write("\t" + VARIANT_IN_ONCOKB_HEADER)
        outf.write("\tMUTATION_EFFECT")
        outf.write("\tONCOGENIC")
        for l in levels:
            outf.write('\t' + l)
        outf.write("\tHIGHEST_LEVEL")
        outf.write("\tCITATIONS\n")

        igene1 = geIndexOfHeader(headers, ['GENE1'])
        igene2 = geIndexOfHeader(headers, ['GENE2'])
        ifusion = geIndexOfHeader(headers, FUSION_HEADERS)
        isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
        icancertype = geIndexOfHeader(headers, CANCER_TYPE_HEADERS)

        i = 0
        queries = []
        rows = []
        for row in reader:
            i = i + 1
            if i % POST_QUERIES_THRESHOLD == 0:
                log.info(i)

            row = padrow(row, ncols)

            sample = getsampleid(row[isample])

            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            gene1 = None
            gene2 = None
            if igene1 >= 0:
                gene1 = row[igene1]
            if igene2 >= 0:
                gene2 = row[igene2]
            if igene1 < 0 and igene2 < 0 and ifusion >= 0:
                fusion = row[ifusion]
                gene1, gene2 = getgenesfromfusion(fusion, nameregex)

            cancertype = defaultCancerType
            if icancertype >= 0:
                cancertype = row[icancertype]
            if sample in cancerTypeMap:
                cancertype = cancerTypeMap[sample]
            if cancertype == "":
                log.info("Cancer type for the sample should be defined for a more accurate result\nline %s: %s\n" % (i, row))
                # continueor


            if not retainonlycuratedgenes or gene1 in curatedgenes or gene2 in curatedgenes:
                queries.append(StructuralVariantQuery(gene1, gene2, 'FUSION', cancertype))
                rows.append(row)

                if len(queries) == POST_QUERIES_THRESHOLD:
                    annotations = pull_structural_variant_info(queries)
                    append_annotation_to_file(outf, rows, annotations)
                    queries = []
                    rows = []
            else:
                # Include default Gene in OncoKB and Variant in OncoKB
                row.append(GENE_IN_ONCOKB_DEFAULT + '\t' + VARIANT_IN_ONCOKB_DEFAULT)
                outf.write('\t'.join(row) + "\n")

        if len(queries) > 0:
            annotations = pull_structural_variant_info(queries)
            append_annotation_to_file(outf, rows, annotations)
    outf.close()


def processcnagisticdata(cnafile, outfile, previousoutfile, defaultCancerType, cancerTypeMap, retainonlycuratedgenes, annotate_gain_loss=False):
    CNA_AMPLIFICATION_TXT = 'Amplification'
    CNA_DELETION_TXT = 'Deletion'
    CNA_LOSS_TXT = 'Loss'
    CNA_GAIN_TXT = 'Gain'

    cnaEventMap = {
        "-2": CNA_DELETION_TXT,
        "-1.5": CNA_DELETION_TXT,
        "2": CNA_AMPLIFICATION_TXT
    }

    if annotate_gain_loss:
        cnaEventMap.update({
            "-1": CNA_LOSS_TXT,
            "1": CNA_GAIN_TXT
        })

    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+', 1000)
    with open(cnafile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')
        headers = readheaders(reader)
        startofsamples = getfirstcolumnofsampleingisticdata(headers['^-$'].split('\t'))
        # Get the header again without upper case the sample id
        infile.seek(0)
        headers = readheaders(reader, startofsamples)
        rawsamples = headers['^-$'].split('\t')[startofsamples:]
        samples = []
        for rs in rawsamples:
            samples.append(getsampleid(rs))

        if defaultCancerType == '' and not set(cancerTypeMap.keys()).issuperset(set(samples)):
            log.info(
                "Cancer type for all samples should be defined for a more accurate result\nsamples in cna file: %s\n" % (
                    samples))

        outf.write('SAMPLE_ID\tCANCER_TYPE\tHUGO_SYMBOL\tALTERATION')
        outf.write("\t"+GENE_IN_ONCOKB_HEADER)
        outf.write("\t"+VARIANT_IN_ONCOKB_HEADER)
        outf.write("\tMUTATION_EFFECT")
        outf.write("\tONCOGENIC")
        for l in levels:
            outf.write('\t' + l)
        outf.write("\tHIGHEST_LEVEL")
        outf.write("\tCITATIONS\n")

        i = 0
        rows = []
        queries = []
        for row in reader:
            i = i + 1
            if i % POST_QUERIES_THRESHOLD == 0:
                log.info(i)

            hugo = row[0]

            for rawsample in rawsamples:
                if rawsample in headers:
                    cna = row[headers[rawsample]]
                    if cna in cnaEventMap:
                        cna_type = cnaEventMap[cna]
                        if cna_type is not None:
                            cancertype = defaultCancerType
                            sample = getsampleid(rawsample)

                            if sampleidsfilter and sample not in sampleidsfilter:
                                continue

                            if sample in cancerTypeMap:
                                cancertype = cancerTypeMap[sample]

                            if not retainonlycuratedgenes or hugo in curatedgenes:
                                rows.append([sample, cancertype, hugo, cna_type])
                                queries.append(CNAQuery(hugo, cna_type, cancertype))

                                if len(queries) == POST_QUERIES_THRESHOLD:
                                    annotations = pull_cna_info(queries)
                                    append_annotation_to_file(outf, rows, annotations)
                                    rows = []
                                    queries = []
                            else:
                                # Include Gene in OncoKB and Variant in OncoKB
                                append_annotation_to_file(outf, [[sample, cancertype, hugo, cna_type]],
                                                          [GENE_IN_ONCOKB_DEFAULT + '\t' + VARIANT_IN_ONCOKB_DEFAULT])

        if len(queries) > 0:
            annotations = pull_cna_info(queries)
            append_annotation_to_file(outf, rows, annotations)

    outf.close()

def getfirstcolumnofsampleingisticdata(headers):
    header0 = headers[0].lower()
    if header0 != "hugo_symbol" and header0 != "gene symbol":
        log.info("Gistic data should start with Hugo_Symbol")
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
        with open(annotatedmutfile, 'rU') as mutfile:
            reader = csv.reader(mutfile, delimiter='\t')
            headers = readheaders(reader)

            ncols = headers["length"]

            igene1 = geIndexOfHeader(headers, ['GENE1'] + HUGO_HEADERS)  # fusion
            igene2 = geIndexOfHeader(headers, ['GENE2'] + HUGO_HEADERS)  # fusion

            ihugo = geIndexOfHeader(headers, HUGO_HEADERS)
            iconsequence = geIndexOfHeader(headers, CONSEQUENCE_HEADERS)
            ihgvs = geIndexOfHeader(headers, HGVS_HEADERS)
            isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
            istart = geIndexOfHeader(headers, START_HEADERS)
            iend = geIndexOfHeader(headers, END_HEADERS)
            icancertype = geIndexOfHeader(headers, CANCER_TYPE_HEADERS)
            # imutationeffect = headers['MUTATION_EFFECT']
            ioncogenic = headers['ONCOGENIC']

            isfusion = igene1 != -1 & igene2 != -1
            ismutorcna = ihugo != -1 & ihgvs != -1

            if not isfusion and not ismutorcna:
                log.error("missing proper header")
                exit()

            for row in reader:

                row = padrow(row, ncols)

                sample = getsampleid(row[isample])

                oncogenic = ""
                if ioncogenic < len(row):
                    oncogenic = row[ioncogenic].lower()
                if sample not in samplelevels:
                    samplelevels[sample] = {}
                    sampleleveltreatments[sample] = {}
                    sampledrivers[sample] = []
                    sampleactionablecount[sample] = {}

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
                            sampleactionablecount[sample][variant] = True

    outf = open(outfile, 'w+')

    # export to anntoated file
    with open(clinicalfile, 'rU') as clinfile:
        reader = csv.reader(clinfile, delimiter='\t')
        headers = readheaders(reader)
        outf.write(headers['^-$'])
        for l in levels:
            outf.write('\t' + l)
        outf.write('\tHIGHEST_LEVEL\tONCOGENIC_MUTATIONS\t#ACTIONABLE_MUTATIONS\t#ONCOGENIC_MUTATIONS\n')
        isample = headers['SAMPLE_ID']

        for row in reader:
            sample = row[isample]

            if sampleidsfilter and sample not in sampleidsfilter:
                continue

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
                actionablecount = len(sampleactionablecount[sample].keys())

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
        isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
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
        'LEVEL_2': '#1F78B4',
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
        'LEVEL_2': 'Level 2',
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
            accumlevelcancerperc = list(map(sum, zip(accumlevelcancerperc,levelcancerperc)))

        ax = plt.gca()
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
        plt.gcf().text(0.90, 0.1, "Generated by OncoKB\n[Chakravarty et al., JCO PO 2017]", fontsize=6,
                       horizontalalignment='right', verticalalignment='bottom')

        # plt.show()
        f.savefig(outfile, bbox_inches='tight')

def processmutationdata(mutfile, outfile, clinicaldata):
    outf = open(outfile, 'w+')
    with open(mutfile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')
        headers = readheaders(reader)

        ihugo = geIndexOfHeader(headers, HUGO_HEADERS)
        iconsequence = geIndexOfHeader(headers, CONSEQUENCE_HEADERS)
        ihgvs = geIndexOfHeader(headers, HGVS_HEADERS)
        isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
        istart = geIndexOfHeader(headers, START_HEADERS)
        iend = geIndexOfHeader(headers, END_HEADERS)

        i = 0
        for row in reader:
            if i % POST_QUERIES_THRESHOLD == 0:
                log.info(i)
            i = i + 1

            sample = row[isample]
            hugo = row[ihugo]
            consequence = row[iconsequence]
            if consequence == 'NULL' or consequence == '':
                consequence = None
            if consequence in mutationtypeconsequencemap:
                consequence = '%2B'.join(mutationtypeconsequencemap[consequence])

            hgvs = row[ihgvs]
            if hgvs.startswith('p.'):
                hgvs = hgvs[2:]
            cancertype = None

            start = row[istart]
            if start == 'NULL' or start == '':
                start = None

            end = row[iend]
            if end == 'NULL' or end == '':
                end = None

            if sample in clinicaldata:
                cancertype = clinicaldata[sample]
            oncokbevidences = pull_mutation_info(hugo, hgvs, consequence, start, end, cancertype)
            annotatedrow = [hugo, consequence, start, end, hgvs, sample, cancertype, oncokbevidences]
            outf.write('\t'.join(annotatedrow) + "\n")

    outf.close()


oncokbcache = {}


def cacheannotated(annotatedfile, defaultCancerType, cancerTypeMap):
    with open(annotatedfile, 'rU') as infile:
        try:
            reader = csv.reader(infile, delimiter='\t')
            headers = readheaders(reader)

            ihugo = geIndexOfHeader(headers, HUGO_HEADERS)
            iconsequence = geIndexOfHeader(headers, CONSEQUENCE_HEADERS)
            ihgvs = geIndexOfHeader(headers, HGVS_HEADERS)
            isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
            istart = geIndexOfHeader(headers, START_HEADERS)
            iend = geIndexOfHeader(headers, END_HEADERS)
            icancertype = geIndexOfHeader(headers, CANCER_TYPE_HEADERS)
            imutationeffect = headers['MUTATION_EFFECT']
            icitations = headers['CITATIONS']
            ioncogenic = headers['ONCOGENIC']
            igeneannotated = headers[GENE_IN_ONCOKB_HEADER]
            ivariantannotated = headers[VARIANT_IN_ONCOKB_HEADER]

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
                    oncokbcache[key][GENE_IN_ONCOKB_HEADER] = row[igeneannotated]
                    oncokbcache[key][VARIANT_IN_ONCOKB_HEADER] = row[ivariantannotated]
                    oncokbcache[key]['mutation_effect'] = row[imutationeffect]
                    oncokbcache[key]['citations'] = row[icitations]
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
        log.error("%s: %s-%s" % (hugo, str(start), str(end)))
    return ""


def pull3dhotspots(hugo, proteinchange, alterationtype, consequence, start, end, cancertype):
    try:
        if hugo in _3dhotspots and consequence == "missense_variant":
            for i in range(int(start), int(end) + 1):
                if i in _3dhotspots[hugo]:
                    return "Y"
    except TypeError:
        log.error("%s: %s-%s" % (hugo, str(start), str(end)))
    return ""

def appendoncokbcitations(citations, pmids, abstracts):
    if citations is None:
        citations = []

    if pmids is not None:
        for pmid in pmids:
            if pmid not in citations:
                citations.append(pmid)

    if abstracts is not None:
        for abstract in abstracts:
            abstractStr = abstract['abstract'] + '(' + abstract['link'] + ')'
            if abstractStr not in citations:
                citations.append(abstractStr)

    return citations


class Gene:
    def __init__(self, hugo):
        self.hugoSymbol = hugo


class ProteinChangeQuery:
    def __init__(self, hugo, hgvs, cancertype, consequence=None, start=None, end=None):
        self.gene = Gene(hugo)
        self.alteration = hgvs
        self.consequence = consequence
        self.proteinStart = start
        self.proteinEnd = end
        self.tumorType = cancertype

class CNAQuery:
    def __init__(self, hugo, cnatype, cancertype):
        self.gene = Gene(hugo)
        self.copyNameAlterationType = cnatype.upper()
        self.tumorType = cancertype

class StructuralVariantQuery:
    def __init__(self, hugoA, hugoB, structural_variant_type, cancertype):

        # Assume all structural variants in the file are functional fusions
        is_functional_fusion = True
        if hugoA == hugoB:
            is_functional_fusion = False
            structural_variant_type = 'DELETION'

        self.geneA = Gene(hugoA)
        self.geneB = Gene(hugoB)
        self.functionalFusion = is_functional_fusion
        self.structuralVariantType = structural_variant_type.upper()
        self.tumorType = cancertype


def pull_mutation_info(queries):
    url = 'https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange'
    response = makeoncokbpostrequest(url, queries)
    annotation = []
    if response.status_code == 200:
        annotation = response.json()
    else:
        for query in queries:
            geturl = url + '?'
            geturl += 'hugoSymbol=' + query.gene.hugoSymbol
            geturl += '&alteration=' + query.alteration
            geturl += '&tumorType=' + query.tumorType
            if query.consequence:
                geturl += '&consequence=' + query.consequence
            if query.proteinStart and query.proteinStart != '\\N' and query.proteinStart != 'NULL' and query.proteinStart != '':
                geturl += '&proteinStart=' + str(query.proteinStart)
            if query.proteinEnd and query.proteinEnd != '\\N' and query.proteinEnd != 'NULL' and query.proteinEnd != '':
                geturl += '&proteinEnd=' + str(query.proteinEnd)
            getresponse = makeoncokbgetrequest(geturl)
            if getresponse.status_code == 200:
                annotation.append(getresponse.json())
            else:
                # if the api call fails, we should still push a None into the list
                # to keep the same length of the queries
                annotation.append(None)

    processed_annotation = []
    for query_annotation in annotation:
        processed_annotation.append(process_oncokb_annotation(query_annotation))
    return processed_annotation


def pull_cna_info(queries):
    url = oncokbapiurl + '/annotate/copyNumberAlterations?'

    response = makeoncokbpostrequest(url, queries)
    annotation = []
    if response.status_code == 200:
        annotation = response.json()
    else:
        for query in queries:
            geturl = url + '?'
            geturl += 'hugoSymbol=' + query.gene.hugoSymbol
            geturl += '&copyNameAlterationType=' + query.copyNameAlterationType
            geturl += '&tumorType=' + query.tumorType
            getresponse = makeoncokbgetrequest(geturl)
            if getresponse.status_code == 200:
                annotation.append(getresponse.json())
            else:
                # if the api call fails, we should still push a None into the list
                # to keep the same length of the queries
                annotation.append(None)

    processed_annotation = []
    for query_annotation in annotation:
        processed_annotation.append(process_oncokb_annotation(query_annotation))
    return processed_annotation



def pull_structural_variant_info(queries):
    url = oncokbapiurl + '/annotate/structuralVariants'

    response = makeoncokbpostrequest(url, queries)
    annotation = []
    if response.status_code == 200:
        annotation = response.json()
    else:
        for query in queries:
            geturl = url + '?'
            geturl += 'hugoSymbolA=' + query.geneA.hugoSymbol
            geturl += '&hugoSymbolB=' + query.geneB.hugoSymbol
            geturl += '&structuralVariantType=' + query.structuralVariantType
            geturl += '&isFunctionalFusion=' + str(query.functionalFusion).upper() if type(query.functionalFusion) is bool else query.functionalFusion
            geturl += '&tumorType=' + query.tumorType

            getresponse = makeoncokbgetrequest(geturl)
            if getresponse.status_code == 200:
                annotation.append(getresponse.json())
            else:
                # if the api call fails, we should still push a None into the list
                # to keep the same length of the queries
                annotation.append(None)

    processed_annotation = []
    for query_annotation in annotation:
        processed_annotation.append(process_oncokb_annotation(query_annotation))
    return processed_annotation



def process_oncokb_annotation(annotation):
    if annotation is None:
        return None

    oncokbdata = {}
    for l in levels:
        oncokbdata[l] = []

    oncokbdata[GENE_IN_ONCOKB_HEADER] = GENE_IN_ONCOKB_DEFAULT
    oncokbdata[VARIANT_IN_ONCOKB_HEADER] = VARIANT_IN_ONCOKB_DEFAULT
    oncokbdata['mutation_effect'] = ""
    oncokbdata['citations'] = []
    oncokbdata['oncogenic'] = ""

    try:
        # oncogenic
        oncokbdata[GENE_IN_ONCOKB_HEADER] = GENE_IN_ONCOKB_DEFAULT if annotation['geneExist'] is None else str(annotation['geneExist'])
        oncokbdata[VARIANT_IN_ONCOKB_HEADER] = VARIANT_IN_ONCOKB_DEFAULT if annotation['variantExist'] is None else str(annotation['variantExist'])

        # oncogenic
        oncokbdata['oncogenic'] = annotation['oncogenic']

        # if not evidences['geneExist'] or (not evidences['variantExist'] and not evidences['alleleExist']):
        #     return ''

        # mutation effect
        if (annotation['mutationEffect'] is not None):
            oncokbdata['mutation_effect'] = annotation['mutationEffect']['knownEffect']
            oncokbdata['citations'] = appendoncokbcitations(oncokbdata['citations'],
                                                            annotation['mutationEffect']['citations']['pmids'],
                                                            annotation['mutationEffect']['citations']['abstracts'])

        # oncogenic
        oncokbdata['oncogenic'] = annotation['oncogenic']

        # get treatment
        for treatment in annotation['treatments']:
            level = treatment['level']

            if level not in levels:
                log.info("%s is ignored" % level)
                # oncokbdata[level].append('')
            else:
                drugs = treatment['drugs']

                oncokbdata['citations'] = appendoncokbcitations(oncokbdata['citations'], treatment['pmids'],
                                                                treatment['abstracts'])

                if len(drugs) == 0:
                    oncokbdata[level].append('[NOT SPECIFIED]')
                else:
                    drugnames = []
                    for drug in drugs:
                        drugnames.append(drug['drugName'])
                    oncokbdata[level].append('+'.join(drugnames))
    except:
        log.error("error when processing %s " % annotation)
        # sys.exit()


    ret = []
    ret.append(oncokbdata[GENE_IN_ONCOKB_HEADER])
    ret.append(oncokbdata[VARIANT_IN_ONCOKB_HEADER])
    ret.append(oncokbdata['mutation_effect'])
    ret.append(oncokbdata['oncogenic'])
    for l in levels:
        ret.append(','.join(oncokbdata[l]))
    ret.append(gethighestsensitivitylevel(oncokbdata))
    ret.append(';'.join(oncokbdata['citations']))

    ret = "\t".join(ret)
    ret = ret.encode('ascii', 'ignore').decode('ascii')  # ignore unicode
    return ret


def gethighestsensitivitylevel(oncokbdata):
    r1 = set()
    if "LEVEL_R1" in oncokbdata:
        r1 = set(oncokbdata["LEVEL_R1"])
    for l in levels:
        if l.startswith("LEVEL_R") or l not in oncokbdata or oncokbdata[l] == '':
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


def readheaders(reader, caseinsensitiveIndex=None):
    headers = {}
    for row in reader:
        if not row[0].startswith("#"):
            headers["^-$"] = '\t'.join(row)  # the whole line
            headers["length"] = len(row)
            i = 0
            for h in row:
                header = h if caseinsensitiveIndex is not None and i >= caseinsensitiveIndex else h.upper()
                headers[header] = i
                i = i + 1
            break
    return headers

def padrow(row, n):
    nr = len(row)
    if nr == n:
        return row

    if nr < n:
        return row + [""] * (n - len(row))

    else:  # nr<n
        return row[0:n]
