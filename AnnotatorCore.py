#!/usr/bin/python
import json
import sys
import csv
from enum import Enum

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

# API timeout is set to two minutes
REQUEST_TIMEOUT = 240

csv.field_size_limit(int(ct.c_ulong(-1).value // 2)) # Deal with overflow problem on Windows, https://stackoverflow.co/120m/questions/15063936/csv-error-field-larger-than-field-limit-131072
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

dxLevels = [
    'LEVEL_Dx1',
    'LEVEL_Dx2',
    'LEVEL_Dx3'
]

pxLevels = [
    'LEVEL_Px1',
    'LEVEL_Px2',
    'LEVEL_Px3'
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
ALTERATION_HEADER = 'ALTERATION'
HGVSP_SHORT_HEADER = 'HGVSP_SHORT'
HGVSP_HEADER = 'HGVSP'
HGVSG_HEADER = 'HGVSG'
HGVS_HEADERS = [ALTERATION_HEADER, HGVSP_SHORT_HEADER, HGVSP_HEADER, HGVSG_HEADER, 'AMINO_ACID_CHANGE', 'FUSION']
SAMPLE_HEADERS = ['SAMPLE_ID', 'TUMOR_SAMPLE_BARCODE']
PROTEIN_START_HEADERS = ['PROTEIN_START']
PROTEIN_END_HEADERS = ['PROTEIN_END']
PROTEIN_POSITION_HEADERS = ['PROTEIN_POSITION']
CANCER_TYPE_HEADERS = ['ONCOTREE_CODE', 'CANCER_TYPE']
FUSION_HEADERS = ['FUSION']
REFERENCE_GENOME_HEADERS = ['NCBI_BUILD', 'REFERENCE_GENOME']

# columns for genomic change annotation
GC_CHROMOSOME_HEADER = 'CHROMOSOME'
GC_START_POSITION_HEADER = 'START_POSITION'
GC_END_POSITION_HEADER = 'END_POSITION'
GC_REF_ALLELE_HEADER = 'REFERENCE_ALLELE'
GC_VAR_ALLELE_1_HEADER = 'TUMOR_SEQ_ALLELE1'
GC_VAR_ALLELE_2_HEADER = 'TUMOR_SEQ_ALLELE2'
GENOMIC_CHANGE_HEADERS = [GC_CHROMOSOME_HEADER, GC_START_POSITION_HEADER, GC_END_POSITION_HEADER, GC_REF_ALLELE_HEADER, GC_VAR_ALLELE_1_HEADER, GC_VAR_ALLELE_2_HEADER]


class QueryType(Enum):
    HGVSP_SHORT = 'HGVSP_SHORT'
    HGVSP = 'HGVSP'
    HGVSG = 'HGVSG'
    GENOMIC_CHANGE = 'GENOMIC_CHANGE'


class ReferenceGenome(Enum):
    GRCH37 = 'GRCh37'
    GRCH38 = 'GRCh38'


REQUIRED_QUERY_TYPE_COLUMNS = {
    QueryType.HGVSP_SHORT: [HGVSP_SHORT_HEADER],
    QueryType.HGVSP: [HGVSP_HEADER],
    QueryType.HGVSG: [HGVSG_HEADER],
    QueryType.GENOMIC_CHANGE: GENOMIC_CHANGE_HEADERS
}

POST_QUERIES_THRESHOLD = 200
POST_QUERIES_THRESHOLD_GC_HGVSG = 100

def getOncokbInfo():
    ret = ['Files annotated on ' + date.today().strftime('%m/%d/%Y') + "\nOncoKB API URL: "+oncokbapiurl]
    try:
        info = requests.get(oncokbapiurl + "/info", timeout=REQUEST_TIMEOUT).json()
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
    response = requests.get(url, timeout=REQUEST_TIMEOUT)
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
    return requests.post(url, headers=headers, data=json.dumps(body, default=lambda o: o.__dict__),
                         timeout=REQUEST_TIMEOUT)


def makeoncokbgetrequest(url):
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer %s' % oncokbapibearertoken
    }
    return requests.get(url, headers=headers, timeout=REQUEST_TIMEOUT)


_3dhotspots = None

def init_3d_hotspots():
    global _3dhotspots
    _3dhotspots = gethotspots(_3dhotspotsbaseurl+"/api/hotspots/3d", None)


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


def append_annotation_to_file(outf, ncols, rows, annotations):
    if len(rows) != len(annotations):
        log.error('The length of the rows and annotations do not match')

    for index, annotation in enumerate(annotations):
        row = rows[index]
        if annotation is not None:
            row = row + annotation

        row = padrow(row, ncols)
        rowstr = '\t'.join(row)
        rowstr = rowstr.encode('ascii', 'ignore').decode('ascii')
        outf.write(rowstr + "\n")


def get_tumor_type_from_row(row, row_index, defaultCancerType, icancertype, cancerTypeMap, sample):
    cancertype = defaultCancerType
    if icancertype >= 0:
        row_cancer_type = get_cell_content(row, icancertype)
        if row_cancer_type is not None:
            cancertype = row_cancer_type
    if sample in cancerTypeMap:
        cancertype = cancerTypeMap[sample]
    if cancertype == "":
        log.info("Cancer type for the sample should be defined for a more accurate result\nline %s: %s\n" % (row_index, row))
        # continue
    return cancertype

def has_desired_headers(desired_headers, file_headers):
    has_required_headers = True
    for header in desired_headers:
        if header not in file_headers:
            has_required_headers = False
            break

    return has_required_headers


def resolve_query_type(user_input_query_type, headers):
    selected_query_type = None
    if isinstance(user_input_query_type, QueryType):
        selected_query_type = user_input_query_type

    if selected_query_type is None and HGVSP_SHORT_HEADER in headers:
        selected_query_type = QueryType.HGVSP_SHORT
    if selected_query_type is None and HGVSP_HEADER in headers:
        selected_query_type = QueryType.HGVSP
    if selected_query_type is None and HGVSG_HEADER in headers:
        selected_query_type = QueryType.HGVSG

    if selected_query_type is None and has_desired_headers(REQUIRED_QUERY_TYPE_COLUMNS[QueryType.GENOMIC_CHANGE], headers):
        selected_query_type = QueryType.GENOMIC_CHANGE

    # default to HGVSp_Short
    if selected_query_type is None:
        selected_query_type = QueryType.HGVSP_SHORT

    # check the file has required columns
    if has_desired_headers(REQUIRED_QUERY_TYPE_COLUMNS[selected_query_type], headers) == False:
        # when it is False, it will never be GENOMIC_CHANGE. For other types, we need to check whether ALTERATION column is available
        if ALTERATION_HEADER not in headers:
            raise Exception("The file does not have required columns "
                            + ', '.join(REQUIRED_QUERY_TYPE_COLUMNS[user_input_query_type])
                            + " for the query type: " + user_input_query_type.value)

    return selected_query_type


def get_reference_genome_from_row(row_reference_genome, default_reference_genome):
    reference_genome = default_reference_genome
    if row_reference_genome is not None and row_reference_genome != '':
        try:
            reference_genome = ReferenceGenome[row_reference_genome.upper()]
        except KeyError:
            log.warning('Unexpected reference genome, only GRCh37 and GRCh38 are supported.' + (
                ' Use default.' if default_reference_genome is not None else ' Skipping.'))
    return reference_genome


def processalterationevents(eventfile, outfile, previousoutfile, defaultCancerType, cancerTypeMap,
                            annotatehotspots, user_input_query_type, default_reference_genome):
    if annotatehotspots:
        init_3d_hotspots()
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+', 1000)
    with open(eventfile, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)

        ncols = headers["length"]
        if ncols == 0:
            return
        newncols = 0

        outf.write(headers['^-$'])

        if annotatehotspots:
            outf.write("\tIS-A-HOTSPOT")
            outf.write("\tIS-A-3D-HOTSPOT")
            newncols += 2

        outf.write("\t" + GENE_IN_ONCOKB_HEADER)
        outf.write("\t" + VARIANT_IN_ONCOKB_HEADER)

        outf.write("\tMUTATION_EFFECT")
        outf.write("\tMUTATION_EFFECT_CITATIONS")
        outf.write("\tONCOGENIC")

        newncols += 5

        for l in levels:
            outf.write('\t' + l)
        newncols += len(levels)

        outf.write("\tHIGHEST_LEVEL")
        outf.write("\tTX_CITATIONS")
        newncols += 2

        for l in dxLevels:
            outf.write('\t' + l)
        newncols += len(dxLevels)

        outf.write("\tHIGHEST_DX_LEVEL")
        outf.write("\tDX_CITATIONS")
        newncols += 2

        for l in pxLevels:
            outf.write('\t' + l)
        newncols += len(pxLevels)

        outf.write("\tHIGHEST_PX_LEVEL")
        outf.write("\tPX_CITATIONS")
        newncols += 2

        outf.write("\n")

        query_type = resolve_query_type(user_input_query_type, headers)
        if (query_type == QueryType.HGVSP_SHORT):
            process_alteration(reader, outf, headers, [HGVSP_SHORT_HEADER, ALTERATION_HEADER], ncols, newncols,
                               defaultCancerType,
                               cancerTypeMap, annotatehotspots, default_reference_genome)

        if (query_type == QueryType.HGVSP):
            process_alteration(reader, outf, headers, [HGVSP_HEADER, ALTERATION_HEADER], ncols, newncols, defaultCancerType,
                               cancerTypeMap, annotatehotspots, default_reference_genome)

        if (query_type == QueryType.HGVSG):
            process_hvsg(reader, outf, headers, [HGVSG_HEADER, ALTERATION_HEADER], ncols, newncols, defaultCancerType,
                         cancerTypeMap, annotatehotspots, default_reference_genome)

        if (query_type == QueryType.GENOMIC_CHANGE):
            process_genomic_change(reader, outf, headers, ncols, newncols, defaultCancerType, cancerTypeMap, annotatehotspots, default_reference_genome)

    outf.close()


def get_cell_content(row, index, return_empty_string=False):
    if index >= 0 and row[index] != 'NULL' and row[index] != '':
        return row[index]
    elif return_empty_string:
        return ''
    else:
        return None

def process_alteration(maffilereader, outf, maf_headers, alteration_column_names, ncols, nannotationcols, defaultCancerType, cancerTypeMap,
                       annotatehotspots, default_reference_genome):
    ihugo = geIndexOfHeader(maf_headers, HUGO_HEADERS)
    iconsequence = geIndexOfHeader(maf_headers, CONSEQUENCE_HEADERS)
    ihgvs = geIndexOfHeader(maf_headers, alteration_column_names)
    isample = geIndexOfHeader(maf_headers, SAMPLE_HEADERS)
    istart = geIndexOfHeader(maf_headers, PROTEIN_START_HEADERS)
    iend = geIndexOfHeader(maf_headers, PROTEIN_END_HEADERS)
    iproteinpos = geIndexOfHeader(maf_headers, PROTEIN_POSITION_HEADERS)
    icancertype = geIndexOfHeader(maf_headers, CANCER_TYPE_HEADERS)
    ireferencegenome= geIndexOfHeader(maf_headers, REFERENCE_GENOME_HEADERS)

    posp = re.compile('[0-9]+')

    i = 0
    queries = []
    rows = []
    for row in maffilereader:
        i = i + 1

        if i % POST_QUERIES_THRESHOLD == 0:
            log.info(i)

        row = padrow(row, ncols)

        sample = row[isample]
        if sampleidsfilter and sample not in sampleidsfilter:
            continue

        hugo = row[ihugo]

        consequence = get_cell_content(row, iconsequence)
        if consequence in mutationtypeconsequencemap:
            consequence = '%2B'.join(mutationtypeconsequencemap[consequence])

        hgvs = row[ihgvs]
        if hgvs.startswith('p.'):
            hgvs = hgvs[2:]

        cancertype = get_tumor_type_from_row(row, i, defaultCancerType, icancertype, cancerTypeMap, sample)
        reference_genome = get_reference_genome_from_row(get_cell_content(row, ireferencegenome), default_reference_genome)

        hgvs = conversion(hgvs)

        start = get_cell_content(row, istart)

        end = get_cell_content(row, iend)

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

        query = ProteinChangeQuery(hugo, hgvs, cancertype, reference_genome, consequence, start, end)
        queries.append(query)
        rows.append(row)

        if len(queries) == POST_QUERIES_THRESHOLD:
            annotations = pull_protein_change_info(queries,annotatehotspots)
            append_annotation_to_file(outf, ncols + nannotationcols, rows, annotations)
            queries = []
            rows = []

    if len(queries) > 0:
        annotations = pull_protein_change_info(queries,annotatehotspots)
        append_annotation_to_file(outf, ncols + nannotationcols, rows, annotations)

# this method is from genome-nexus annotation-tools
# https://github.com/genome-nexus/annotation-tools/blob/53ff7f7fe673e961282f871ebc78d2ecc0831919/standardize_mutation_data.py
def get_var_allele(ref_allele, tumor_seq_allele1, tumor_seq_allele2):
    # set the general tumor_seq_allele as the first non-ref allele encountered
    # this will be used to resolve the variant classification and variant type
    # if there are no tumor alleles that do not match the ref allele then use empty string
    # in the event that this happens then there might be something wrong with the data itself
    # if both alleles are different, use allele2. Stick with the logic of GenomeNexus
    try:
        tumor_seq_allele = ""
        if ref_allele != tumor_seq_allele2:
            tumor_seq_allele = tumor_seq_allele2
        elif ref_allele != tumor_seq_allele1:
            tumor_seq_allele = tumor_seq_allele1
    except:
        tumor_seq_allele = ""

    return tumor_seq_allele

def process_genomic_change(maffilereader, outf, maf_headers, ncols, nannotationcols, defaultCancerType, cancerTypeMap, annotatehotspots, default_reference_genome):
    ichromosome = geIndexOfHeader(maf_headers, [GC_CHROMOSOME_HEADER])
    istart = geIndexOfHeader(maf_headers, [GC_START_POSITION_HEADER])
    iend = geIndexOfHeader(maf_headers, [GC_END_POSITION_HEADER])
    irefallele = geIndexOfHeader(maf_headers, [GC_REF_ALLELE_HEADER])
    ivarallele1 = geIndexOfHeader(maf_headers, [GC_VAR_ALLELE_1_HEADER])
    ivarallele2 = geIndexOfHeader(maf_headers, [GC_VAR_ALLELE_2_HEADER])

    isample = geIndexOfHeader(maf_headers, SAMPLE_HEADERS)
    icancertype = geIndexOfHeader(maf_headers, CANCER_TYPE_HEADERS)
    ireferencegenome= geIndexOfHeader(maf_headers, REFERENCE_GENOME_HEADERS)

    posp = re.compile('[0-9]+')

    i = 0
    queries = []
    rows = []
    for row in maffilereader:
        i = i + 1

        if i % POST_QUERIES_THRESHOLD_GC_HGVSG == 0:
            log.info(i)

        row = padrow(row, ncols)

        sample = row[isample]
        if sampleidsfilter and sample not in sampleidsfilter:
            continue

        cancertype = get_tumor_type_from_row(row, i, defaultCancerType, icancertype, cancerTypeMap, sample)
        reference_genome = get_reference_genome_from_row(get_cell_content(row, ireferencegenome), default_reference_genome)

        chromosome = get_cell_content(row, ichromosome, True)
        start = get_cell_content(row, istart, True)
        end = get_cell_content(row, iend, True)
        ref_allele = get_cell_content(row, irefallele, True)
        var_allele_1 = get_cell_content(row, ivarallele1, True)
        var_allele_2 = get_cell_content(row, ivarallele2, True)
        var_allele = get_var_allele(ref_allele, var_allele_1, var_allele_2)

        query = GenomicChangeQuery(chromosome, start, end, ref_allele, var_allele, cancertype, reference_genome)
        queries.append(query)
        rows.append(row)

        if len(queries) == POST_QUERIES_THRESHOLD_GC_HGVSG:
            annotations = pull_genomic_change_info(queries,annotatehotspots)
            append_annotation_to_file(outf, ncols+nannotationcols, rows, annotations)
            queries = []
            rows = []

    if len(queries) > 0:
        annotations = pull_genomic_change_info(queries,annotatehotspots)
        append_annotation_to_file(outf, ncols+nannotationcols, rows, annotations)

def process_hvsg(maffilereader, outf, maf_headers, alteration_column_names, ncols, nannotationcols, defaultCancerType, cancerTypeMap, annotatehotspots, default_reference_genome):
    ihgvsg = geIndexOfHeader(maf_headers, alteration_column_names)
    isample = geIndexOfHeader(maf_headers, SAMPLE_HEADERS)
    icancertype = geIndexOfHeader(maf_headers, CANCER_TYPE_HEADERS)
    ireferencegenome= geIndexOfHeader(maf_headers, REFERENCE_GENOME_HEADERS)

    i = 0
    queries = []
    rows = []
    for row in maffilereader:
        i = i + 1

        if i % POST_QUERIES_THRESHOLD_GC_HGVSG == 0:
            log.info(i)

        row = padrow(row, ncols)

        sample = row[isample]
        if sampleidsfilter and sample not in sampleidsfilter:
            continue

        hgvsg = get_cell_content(row, ihgvsg)

        cancertype = get_tumor_type_from_row(row, i, defaultCancerType, icancertype, cancerTypeMap, sample)
        reference_genome = get_reference_genome_from_row(get_cell_content(row, ireferencegenome), default_reference_genome)

        if hgvsg is None:
            if annotatehotspots:
                default_cols = [['', '', GENE_IN_ONCOKB_DEFAULT, VARIANT_IN_ONCOKB_DEFAULT]]
            else:
                default_cols = [[GENE_IN_ONCOKB_DEFAULT, VARIANT_IN_ONCOKB_DEFAULT]]
            append_annotation_to_file(outf, ncols + nannotationcols, [row],
                                      default_cols)
        else:
            query = HGVSgQuery(hgvsg, cancertype, reference_genome)
            queries.append(query)
            rows.append(row)

        if len(queries) == POST_QUERIES_THRESHOLD_GC_HGVSG:
            annotations = pull_hgvsg_info(queries, annotatehotspots)
            append_annotation_to_file(outf, ncols+nannotationcols, rows, annotations)
            queries = []
            rows = []

    if len(queries) > 0:
        annotations = pull_hgvsg_info(queries,annotatehotspots)
        append_annotation_to_file(outf, ncols+nannotationcols, rows, annotations)


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

def processsv(svdata, outfile, previousoutfile, defaultCancerType, cancerTypeMap, nameregex):
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+')
    with open(svdata, 'rU') as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)

        ncols = headers["length"]

        if ncols == 0:
            return

        outf.write(headers['^-$'])
        outf.write("\t" + GENE_IN_ONCOKB_HEADER)
        outf.write("\t" + VARIANT_IN_ONCOKB_HEADER)
        outf.write("\tMUTATION_EFFECT")
        outf.write("\tMUTATION_EFFECT_CITATIONS")
        outf.write("\tONCOGENIC")
        for l in levels:
            outf.write('\t' + l)
        outf.write("\tHIGHEST_LEVEL")
        outf.write("\tTX_CITATIONS")

        for l in dxLevels:
            outf.write('\t' + l)
        outf.write("\tHIGHEST_DX_LEVEL")
        outf.write("\tDX_CITATIONS")

        for l in pxLevels:
            outf.write('\t' + l)
        outf.write("\tHIGHEST_PX_LEVEL")
        outf.write("\tPX_CITATIONS")
        outf.write("\n")

        newcols = ncols + 11 + len(levels) + len(dxLevels) + len(pxLevels)

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

            sample = row[isample]

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

            cancertype = get_tumor_type_from_row(row, i, defaultCancerType, icancertype, cancerTypeMap, sample)


            queries.append(StructuralVariantQuery(gene1, gene2, 'FUSION', cancertype))
            rows.append(row)

            if len(queries) == POST_QUERIES_THRESHOLD:
                annotations = pull_structural_variant_info(queries)
                append_annotation_to_file(outf, newcols, rows, annotations)
                queries = []
                rows = []

        if len(queries) > 0:
            annotations = pull_structural_variant_info(queries)
            append_annotation_to_file(outf, newcols, rows, annotations)
    outf.close()


def processcnagisticdata(cnafile, outfile, previousoutfile, defaultCancerType, cancerTypeMap, annotate_gain_loss=False):
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
        samples = []
        rawsamples = []
        if headers["length"] != 0:
            startofsamples = getfirstcolumnofsampleingisticdata(headers['^-$'].split('\t'))
            rawsamples = headers['^-$'].split('\t')[startofsamples:]
        for rs in rawsamples:
            samples.append(rs)

        if defaultCancerType == '' and not set(cancerTypeMap.keys()).issuperset(set(samples)):
            log.info(
                "Cancer type for all samples should be defined for a more accurate result\nsamples in cna file: %s\n" % (
                    samples))

        outf.write('SAMPLE_ID\tCANCER_TYPE\tHUGO_SYMBOL\tALTERATION')
        outf.write("\t"+GENE_IN_ONCOKB_HEADER)
        outf.write("\t"+VARIANT_IN_ONCOKB_HEADER)
        outf.write("\tMUTATION_EFFECT")
        outf.write("\tMUTATION_EFFECT_CITATIONS")
        outf.write("\tONCOGENIC")
        for l in levels:
            outf.write('\t' + l)
        outf.write("\tHIGHEST_LEVEL")
        outf.write("\tTX_CITATIONS")

        for l in dxLevels:
            outf.write('\t' + l)
        outf.write("\tHIGHEST_DX_LEVEL")
        outf.write("\tDX_CITATIONS")

        for l in pxLevels:
            outf.write('\t' + l)
        outf.write("\tHIGHEST_PX_LEVEL")
        outf.write("\tPX_CITATIONS")
        outf.write("\n")

        ncols = 15 + len(levels) + len(dxLevels) + len(pxLevels)

        i = 0
        rows = []
        queries = []
        for row in reader:
            i = i + 1
            if i % POST_QUERIES_THRESHOLD == 0:
                log.info(i)

            hugo = row[0]
            if len(row) == 1:
                log.warning("No CNA specified for gene " + hugo)
                continue

            for rawsample in rawsamples:
                if rawsample in headers:
                    if len(row) <= headers[rawsample]:
                        log.warning('No CNA specified for ' + row[0] + ' ' + rawsample)
                        continue
                    cna = row[headers[rawsample]]
                    if cna in cnaEventMap:
                        cna_type = cnaEventMap[cna]
                        if cna_type is not None:
                            cancertype = defaultCancerType
                            sample = rawsample

                            if sampleidsfilter and sample not in sampleidsfilter:
                                continue

                            if sample in cancerTypeMap:
                                cancertype = cancerTypeMap[sample]

                            rows.append([sample, cancertype, hugo, cna_type])
                            queries.append(CNAQuery(hugo, cna_type, cancertype))

                            if len(queries) == POST_QUERIES_THRESHOLD:
                                annotations = pull_cna_info(queries)
                                append_annotation_to_file(outf, ncols, rows, annotations)
                                rows = []
                                queries = []

        if len(queries) > 0:
            annotations = pull_cna_info(queries)
            append_annotation_to_file(outf, ncols, rows, annotations)

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
    sampledxlevels = {}
    samplepxlevels = {}
    sampleleveltreatments = {}
    sampledrivers = {}
    samplemutationswithdiagnosis = {}
    samplemutationswithprognosis = {}
    sampleactionablecount = {}
    samplealterationcount = {}
    for annotatedmutfile in annotatedmutfiles:
        with open(annotatedmutfile, 'rU') as mutfile:
            reader = csv.reader(mutfile, delimiter='\t')
            headers = readheaders(reader)

            ncols = headers["length"]

            if ncols == 0:
                return

            igene1 = geIndexOfHeader(headers, ['GENE1'] + HUGO_HEADERS)  # fusion
            igene2 = geIndexOfHeader(headers, ['GENE2'] + HUGO_HEADERS)  # fusion
            ifusion = geIndexOfHeader(headers, ['FUSION'])

            ihugo = geIndexOfHeader(headers, HUGO_HEADERS)
            iconsequence = geIndexOfHeader(headers, CONSEQUENCE_HEADERS)
            ihgvs = geIndexOfHeader(headers, HGVS_HEADERS)
            isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
            istart = geIndexOfHeader(headers, PROTEIN_START_HEADERS)
            iend = geIndexOfHeader(headers, PROTEIN_END_HEADERS)
            icancertype = geIndexOfHeader(headers, CANCER_TYPE_HEADERS)
            # imutationeffect = headers['MUTATION_EFFECT']
            ioncogenic = headers['ONCOGENIC']

            isfusion = (igene1 != -1 & igene2 != -1) or ifusion != -1
            ismutorcna = ihugo != -1 & ihgvs != -1

            if not isfusion and not ismutorcna:
                log.error("missing proper header")
                exit()

            for row in reader:

                row = padrow(row, ncols)

                sample = row[isample]

                oncogenic = ""
                if ioncogenic < len(row):
                    oncogenic = row[ioncogenic].lower()
                if sample not in samplelevels:
                    samplelevels[sample] = {}
                    sampledxlevels[sample] = []
                    samplepxlevels[sample] = []
                    sampleleveltreatments[sample] = {}
                    sampledrivers[sample] = []
                    sampleactionablecount[sample] = {}

                if sample not in samplemutationswithdiagnosis:
                    samplemutationswithdiagnosis[sample] = []

                if sample not in samplemutationswithprognosis:
                    samplemutationswithprognosis[sample] = []

                if sample not in samplealterationcount:
                    samplealterationcount[sample] = 1
                else:
                    samplealterationcount[sample] += 1

                hugo = row[ihugo]
                alteration = row[ihgvs]
                gene1 = row[igene1]
                gene2 = row[igene2]

                variant = "NA"
                if ismutorcna:
                    variant = hugo + " " + alteration
                elif isfusion:
                    if ifusion != -1:
                        variant = row[ifusion]
                    else:
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

                for l in dxLevels:
                    il = headers[l]
                    if il < len(row) and row[il] != '':
                        if l not in samplelevels[sample]:
                            samplelevels[sample][l] = []
                        samplelevels[sample][l].append(row[il] + "(" + variant + ")")

                for l in pxLevels:
                    il = headers[l]
                    if il < len(row) and row[il] != '':
                        if l not in samplelevels[sample]:
                            samplelevels[sample][l] = []
                        samplelevels[sample][l].append(row[il] + "(" + variant + ")")

                ihighestdxlevel = geIndexOfHeader(headers, ['HIGHEST_DX_LEVEL'])
                if ihighestdxlevel != -1:
                    if row[ihighestdxlevel] != '':
                        samplemutationswithdiagnosis[sample].append(variant)
                        sampledxlevels[sample].append(row[ihighestdxlevel])

                ihighestpxlevel = geIndexOfHeader(headers, ['HIGHEST_PX_LEVEL'])
                if ihighestpxlevel != -1:
                    if row[ihighestpxlevel] != '':
                        samplemutationswithprognosis[sample].append(variant)
                        samplepxlevels[sample].append(row[ihighestpxlevel])

    outf = open(outfile, 'w+')

    # export to anntoated file
    with open(clinicalfile, 'rU') as clinfile:
        reader = csv.reader(clinfile, delimiter='\t')
        headers = readheaders(reader)
        outf.write(headers['^-$'])
        for l in levels:
            outf.write('\t' + l)
        outf.write('\tHIGHEST_LEVEL')
        for l in dxLevels:
            outf.write('\t' + l)
        outf.write('\tHIGHEST_DX_LEVEL')
        for l in pxLevels:
            outf.write('\t' + l)
        outf.write('\tHIGHEST_PX_LEVEL')
        outf.write('\tONCOGENIC_MUTATIONS\t#ONCOGENIC_MUTATIONS\t#MUTATIONS_WITH_THERAPEUTIC_IMPLICATIONS\t#MUTATIONS_WITH_DIAGNOSTIC_IMPLICATIONS\t#MUTATIONS_WITH_PROGNOSTIC_IMPLICATIONS\t#MUTATIONS\n')
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
            highestdxlevel = ''
            highestpxlevel = ''
            if sample in sampleleveltreatments:
                highestlevel = gethighestsensitivitylevel(sampleleveltreatments[sample])
            if sample in sampledxlevels:
                highestdxlevel = gethighestDxPxlevel(dxLevels, sampledxlevels[sample])
            if sample in samplepxlevels:
                highestpxlevel = gethighestDxPxlevel(pxLevels, samplepxlevels[sample])
            # if highestlevel == '':
            #     if sample in sampledrivers and len(sampledrivers[sample])>0:
            #         highestlevel = 'Oncogenic, no level'
            #     else:
            #         highestlevel = "VUS"
            outf.write('\t' + highestlevel)

            for l in dxLevels:
                outf.write('\t')
                if sample in samplelevels and l in samplelevels[sample]:
                    outf.write(";".join(samplelevels[sample][l]))

            outf.write('\t' + highestdxlevel)

            for l in pxLevels:
                outf.write('\t')
                if sample in samplelevels and l in samplelevels[sample]:
                    outf.write(";".join(samplelevels[sample][l]))
            outf.write('\t' + highestpxlevel)


            actionablecount = 0
            if sample in sampleactionablecount:
                actionablecount = len(sampleactionablecount[sample].keys())

            alterationcount = 0
            if sample in samplealterationcount:
                alterationcount = samplealterationcount[sample]

            drivercount = 0
            diagnosiscount = 0
            prognosiscount = 0
            drivermutations = ""
            if sample in sampledrivers:
                drivercount = len(sampledrivers[sample])
                drivermutations = ";".join(sampledrivers[sample])
            if sample in samplemutationswithdiagnosis:
                diagnosiscount = len(samplemutationswithdiagnosis[sample])
            if sample in samplemutationswithprognosis:
                prognosiscount = len(samplemutationswithprognosis[sample])

            outf.write('\t' + drivermutations)
            outf.write('\t' + str(drivercount))
            outf.write('\t' + str(actionablecount))
            outf.write('\t' + str(diagnosiscount))
            outf.write('\t' + str(prognosiscount))
            outf.write('\t' + str(alterationcount))

            outf.write('\n')

    outf.close()

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
    order = reversed(sorted(range(ncat),key=lambda x:(catactionabilityarray[x],catoncogenicarray[x])))
    drawplot(ax, 'OncoKB Actionability', extlevels, levelcatsamplecount, catarray, catsamplecount, order, parameters["thresholdcat"])

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
    catarray = [] # cancer types
    catactionabilityarray = [] # actionabiligy percentages per cancer type
    for cat in catsamplecount:
        if catsamplecount[cat] >= parameters["thresholdcat"]:
            catarray.append(cat)
            catactionabilityarray.append(catactionablesamplecount[cat] * 100.0 / catsamplecount[cat])

    ncat = len(catarray)
    order = reversed(sorted(range(ncat),key=lambda x:(catactionabilityarray[x])))
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
        'LEVEL_R3': '#FCD6D3',

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
        'LEVEL_R3': 'Level R3',

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
            accumlevelcancerperc = list(map(sum, zip(accumlevelcancerperc,levelcancerperc)))

        ax = plt.gca()
        ax.set_axisbelow(True)
        ax.set_aspect(0.1)

        ax.tick_params(axis='y', which='major', labelsize=6)
        ax.set_ylabel('% of samples', fontsize=6)
        ax.set_title(title, fontsize=8)
        ax.set_xticks([i+0.5 for i in ind])
        ax.set_xticklabels(catarray, rotation=60, ha="right", fontsize=4)
        # plt.yticks(np.arange(0, 81, 10))
        ax.legend(plts, legends, fontsize=6, bbox_to_anchor=(1.01, 1), loc="upper left")


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
            istart = geIndexOfHeader(headers, PROTEIN_START_HEADERS)
            iend = geIndexOfHeader(headers, PROTEIN_END_HEADERS)
            icancertype = geIndexOfHeader(headers, CANCER_TYPE_HEADERS)
            imutationeffect = headers['MUTATION_EFFECT']
            icitations = headers['CITATIONS']
            ioncogenic = headers['ONCOGENIC']
            igeneannotated = headers[GENE_IN_ONCOKB_HEADER]
            ivariantannotated = headers[VARIANT_IN_ONCOKB_HEADER]

            for row in reader:
                try:
                    hugo = row[ihugo]

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


def pull3dhotspots(hugo, consequence, start, end):
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
    def __init__(self, hugo, hgvs, cancertype, reference_genome=None, consequence=None, start=None, end=None):
        self.gene = Gene(hugo)
        self.alteration = hgvs
        if consequence is not None:
            self.consequence = consequence
        if start is not None:
            self.proteinStart = start
        if end is not None:
            self.proteinEnd = end
        self.tumorType = cancertype
        if reference_genome is not None:
            self.referenceGenome = reference_genome.value


class HGVSgQuery:
    def __init__(self, hgvsg, cancertype, reference_genome=None):
        self.hgvsg = hgvsg
        self.tumorType = cancertype
        if reference_genome is not None:
            self.referenceGenome = reference_genome.value


def gettumortypename(tumortype):
    if 'code' in tumortype and tumortype['code'] is not None and tumortype['code'] != '':
        return tumortype['code']
    elif 'name' in tumortype and tumortype['name'] is not None and tumortype['name'] != '':
        return tumortype['name']
    else:
        return tumortype['mainType']['name']


def getimplications(oncokbdata, implication_type, levels, implications):
    citation_column_key = implication_type + '_citations'
    for implication in implications:
        oncokbdata[citation_column_key] = appendoncokbcitations(oncokbdata[citation_column_key], implication['pmids'],
                                                                implication['abstracts'])
        level = implication['levelOfEvidence']

        if level is not None:
            if level not in levels:
                log.info(level + " is ignored")
            else:
                if 'tumorType' in implication:
                    tumortypename = gettumortypename(implication['tumorType'])
                    if tumortypename not in oncokbdata[level]:
                        oncokbdata[level].append(tumortypename)


class GenomicChangeQuery:
    def __init__(self, chromosome, start, end, ref_allele, var_allele, cancertype, reference_genome=None):
        self.genomicLocation = ','.join([chromosome, start, end, ref_allele, var_allele])
        self.tumorType = cancertype
        if reference_genome is not None:
            self.referenceGenome = reference_genome.value

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


def pull_protein_change_info(queries, annotate_hotspot):
    url = oncokbapiurl + '/annotate/mutations/byProteinChange'
    response = makeoncokbpostrequest(url, queries)
    if response.status_code == 401:
        raise Exception('unauthorized')
    annotation = []
    if response.status_code == 200:
        annotation = response.json()
    else:
        for query in queries:
            geturl = url + '?'
            geturl += 'hugoSymbol=' + query.gene.hugoSymbol
            geturl += '&alteration=' + query.alteration
            geturl += '&tumorType=' + query.tumorType
            if hasattr(query, 'consequence') and query.consequence:
                geturl += '&consequence=' + query.consequence
            if hasattr(query, 'proteinStart') and query.proteinStart and query.proteinStart != '\\N' and query.proteinStart != 'NULL' and query.proteinStart != '':
                geturl += '&proteinStart=' + str(query.proteinStart)
            if hasattr(query, 'proteinEnd') and query.proteinEnd and query.proteinEnd != '\\N' and query.proteinEnd != 'NULL' and query.proteinEnd != '':
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
        processed_annotation.append(process_oncokb_annotation(query_annotation, annotate_hotspot))
    return processed_annotation


def pull_hgvsg_info(queries, annotate_hotspot):
    url = oncokbapiurl + '/annotate/mutations/byHGVSg'
    response = makeoncokbpostrequest(url, queries)
    if response.status_code == 401:
        raise Exception('unauthorized')
    annotation = []
    if response.status_code == 200:
        annotation = response.json()
    else:
        for query in queries:
            geturl = url + '?'
            geturl += 'hgvsg=' + query.hgvsg
            geturl += '&tumorType=' + query.tumorType
            getresponse = makeoncokbgetrequest(geturl)
            if getresponse.status_code == 200:
                annotation.append(getresponse.json())
            else:
                # if the api call fails, we should still push a None into the list
                # to keep the same length of the queries
                print('Error on annotating the url ' + geturl)
                annotation.append(None)

    processed_annotation = []
    for query_annotation in annotation:
        processed_annotation.append(process_oncokb_annotation(query_annotation, annotate_hotspot))
    return processed_annotation

def pull_genomic_change_info(queries, annotate_hotspot):
    url = oncokbapiurl + '/annotate/mutations/byGenomicChange'
    response = makeoncokbpostrequest(url, queries)
    if response.status_code == 401:
        raise Exception('unauthorized')
    annotation = []
    if response.status_code == 200:
        annotation = response.json()
    else:
        for query in queries:
            geturl = url + '?'
            geturl += 'genomicLocation=' + query.genomicLocation
            geturl += '&tumorType=' + query.tumorType
            getresponse = makeoncokbgetrequest(geturl)
            if getresponse.status_code == 200:
                annotation.append(getresponse.json())
            else:
                # if the api call fails, we should still push a None into the list
                # to keep the same length of the queries
                print('Error on annotating the url ' + geturl)
                annotation.append(None)

    processed_annotation = []
    for query_annotation in annotation:
        processed_annotation.append(process_oncokb_annotation(query_annotation, annotate_hotspot))
    return processed_annotation


def pull_cna_info(queries):
    url = oncokbapiurl + '/annotate/copyNumberAlterations'

    response = makeoncokbpostrequest(url, queries)
    if response.status_code == 401:
        raise Exception('unauthorized')
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
                print('Error on annotating the url ' + geturl)
                annotation.append(None)

    processed_annotation = []
    for query_annotation in annotation:
        processed_annotation.append(process_oncokb_annotation(query_annotation, annotate_hotspot=False))
    return processed_annotation



def pull_structural_variant_info(queries):
    url = oncokbapiurl + '/annotate/structuralVariants'

    response = makeoncokbpostrequest(url, queries)
    if response.status_code == 401:
        raise Exception('unauthorized')
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
                print('Error on annotating the url ' + geturl)
                annotation.append(None)

    processed_annotation = []
    for query_annotation in annotation:
        processed_annotation.append(process_oncokb_annotation(query_annotation, annotate_hotspot=False))
    return processed_annotation



def process_oncokb_annotation(annotation, annotate_hotspot):
    if annotation is None:
        return None

    oncokbdata = {}
    for l in levels:
        oncokbdata[l] = []
    for l in dxLevels:
        oncokbdata[l] = []
    for l in pxLevels:
        oncokbdata[l] = []

    oncokbdata[GENE_IN_ONCOKB_HEADER] = GENE_IN_ONCOKB_DEFAULT
    oncokbdata[VARIANT_IN_ONCOKB_HEADER] = VARIANT_IN_ONCOKB_DEFAULT
    oncokbdata['mutation_effect'] = ""
    oncokbdata['mutation_effect_citations'] = []
    oncokbdata['citations'] = []
    oncokbdata['oncogenic'] = ""
    oncokbdata['tx_citations'] = []
    oncokbdata['dx_citations'] = []
    oncokbdata['px_citations'] = []

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
            oncokbdata['mutation_effect_citations'] = appendoncokbcitations(oncokbdata['mutation_effect_citations'],
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

                oncokbdata['tx_citations'] = appendoncokbcitations(oncokbdata['tx_citations'], treatment['pmids'],
                                                                treatment['abstracts'])

                if len(drugs) == 0:
                    oncokbdata[level].append('[NOT SPECIFIED]')
                else:
                    drugnames = []
                    for drug in drugs:
                        drugnames.append(drug['drugName'])
                    treatmentname = '+'.join(drugnames)
                    if treatmentname not in oncokbdata[level]:
                        oncokbdata[level].append('+'.join(drugnames))
        if annotation['diagnosticImplications'] is not None:
            getimplications(oncokbdata, 'dx', dxLevels, annotation['diagnosticImplications'])

        if annotation['prognosticImplications'] is not None:
            getimplications(oncokbdata, 'px', pxLevels, annotation['prognosticImplications'])

        oncokbdata['highestDiagnosticImplicationLevel'] = annotation['highestDiagnosticImplicationLevel']
        oncokbdata['highestPrognosticImplicationLevel'] = annotation['highestPrognosticImplicationLevel']
    except:
        log.error("error when processing %s " % annotation)
        # sys.exit()


    ret = []
    if annotate_hotspot:
        if annotation['hotspot']:
            ret.append('Y')
        else:
            ret.append('')

        _3dhotspot = pull3dhotspots(annotation['query']['hugoSymbol'], annotation['query']['consequence'], annotation['query']['proteinStart'], annotation['query']['proteinEnd'])
        ret.append(_3dhotspot)

    ret.append(oncokbdata[GENE_IN_ONCOKB_HEADER])
    ret.append(oncokbdata[VARIANT_IN_ONCOKB_HEADER])
    ret.append(oncokbdata['mutation_effect'])
    ret.append(';'.join(oncokbdata['mutation_effect_citations']))
    ret.append(oncokbdata['oncogenic'])
    for l in levels:
        ret.append(','.join(oncokbdata[l]))
    ret.append(gethighestsensitivitylevel(oncokbdata))
    ret.append(';'.join(oncokbdata['tx_citations']))

    for l in dxLevels:
        ret.append(','.join(oncokbdata[l]))
    ret.append(gethighestDxPxlevel(dxLevels, [oncokbdata['highestDiagnosticImplicationLevel']]))
    ret.append(';'.join(oncokbdata['dx_citations']))

    for l in pxLevels:
        ret.append(','.join(oncokbdata[l]))
    ret.append(gethighestDxPxlevel(pxLevels, [oncokbdata['highestPrognosticImplicationLevel']]))
    ret.append(';'.join(oncokbdata['px_citations']))

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

def gethighestDxPxlevel(levels, oncokbdata):
    for l in levels:
        if l not in oncokbdata:
            continue
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
    headers["length"] = 0
    for row in reader:
        if not row[0].startswith("#"):
            headers["^-$"] = '\t'.join(row)  # the whole line
            headers["length"] = len(row)
            i = 0
            for h in row:
                h=h.strip()
                headers[h.upper()] = i
                headers[h] = i
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
