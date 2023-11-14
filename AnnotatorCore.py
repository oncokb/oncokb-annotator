#!/usr/bin/python
import datetime
import json
import csv
import sys

import requests
import os.path
import logging
import re
import ctypes as ct

from enum import Enum
from requests.adapters import HTTPAdapter
from urllib3 import Retry
from datetime import date

logging.basicConfig(level=logging.INFO)
logging.getLogger("requests").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

log = logging.getLogger('AnnotatorCore')

# API timeout is set to two minutes
REQUEST_TIMEOUT = 240

API_REQUEST_RETRY_STATUS_FORCELIST = [429, 500, 502, 503, 504]

csv.field_size_limit(int(ct.c_ulong(
    -1).value // 2))  # Deal with overflow problem on Windows, https://stackoverflow.co/120m/questions/15063936/csv-error-field-larger-than-field-limit-131072
sizeLimit = csv.field_size_limit()
csv.field_size_limit(sizeLimit)  # for reading large files

DEFAULT_ONCOKB_URL = "https://www.oncokb.org"
oncokb_api_url = DEFAULT_ONCOKB_URL + "/api"
oncokb_annotation_api_url = oncokb_api_url + "/v1"

oncokb_api_bearer_token = ""

# 'U', which no longer has any effect in python 3. 3.11 completely removed the option.
# It will throw error if we don't do condition check.
# https://stackoverflow.com/questions/56791545/what-is-the-non-deprecated-version-of-open-u-mode
DEFAULT_READ_FILE_MODE = 'r' if sys.version_info.major > 2 else 'rU'


def setoncokbbaseurl(u):
    if u and u is not None:
        global oncokb_api_url
        global oncokb_annotation_api_url
        oncokb_api_url = u.rstrip('/') + '/api'
        oncokb_annotation_api_url = oncokb_api_url + '/v1'


def setoncokbapitoken(t):
    global oncokb_api_bearer_token
    oncokb_api_bearer_token = t.strip()


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


ANNOTATED_HEADER = 'ANNOTATED'
GENE_IN_ONCOKB_HEADER = 'GENE_IN_ONCOKB'
VARIANT_IN_ONCOKB_HEADER = 'VARIANT_IN_ONCOKB'

GENE_IN_ONCOKB_DEFAULT = 'False'
VARIANT_IN_ONCOKB_DEFAULT = 'False'

levels = [
    'LEVEL_R1',
    'LEVEL_1',
    'LEVEL_2',
    'LEVEL_3A',
    'LEVEL_3B',
    'LEVEL_4',
    'LEVEL_R2'
]
sensitive_levels = [
    'LEVEL_1',
    'LEVEL_2',
    'LEVEL_3A',
    'LEVEL_3B',
    'LEVEL_4',
]
resistance_levels = [
    'LEVEL_R1',
    'LEVEL_R2'
]
TX_TYPE_SENSITIVE = 'sensitive'
TX_TYPE_RESISTANCE = 'resistance'

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

CNA_AMPLIFICATION_TXT = 'Amplification'
CNA_DELETION_TXT = 'Deletion'
CNA_LOSS_TXT = 'Loss'
CNA_GAIN_TXT = 'Gain'

CNAS = [
    CNA_DELETION_TXT,
    CNA_LOSS_TXT,
    CNA_GAIN_TXT,
    CNA_AMPLIFICATION_TXT,
]

GISTIC_CNA_MAP = {
    "-2": CNA_DELETION_TXT,
    "-1.5": CNA_DELETION_TXT,
    "-1": CNA_LOSS_TXT,
    "1": CNA_GAIN_TXT,
    "2": CNA_AMPLIFICATION_TXT
}

CNA_FILE_FORMAT_GISTIC = 'gistic'
CNA_FILE_FORMAT_INDIVIDUAL = 'individual'
CND_FILE_FORMAT = [CNA_FILE_FORMAT_GISTIC, CNA_FILE_FORMAT_INDIVIDUAL]

# column headers
HUGO_HEADERS = ['HUGO_SYMBOL', 'HUGO_GENE_SYMBOL', 'GENE']
CONSEQUENCE_HEADERS = ['VARIANT_CLASSIFICATION', 'MUTATION_TYPE']
ALTERATION_HEADER = 'ALTERATION'
HGVSP_SHORT_HEADER = 'HGVSP_SHORT'
HGVSP_HEADER = 'HGVSP'
HGVSG_HEADER = 'HGVSG'
# columns for copy number alteration
CNA_HEADERS = [ALTERATION_HEADER, 'COPY_NUMBER_ALTERATION', 'CNA', 'GISTIC']
HGVS_HEADERS = [ALTERATION_HEADER, HGVSP_SHORT_HEADER, HGVSP_HEADER, HGVSG_HEADER, 'AMINO_ACID_CHANGE',
                'FUSION'] + CNA_HEADERS
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
GENOMIC_CHANGE_HEADERS = [GC_CHROMOSOME_HEADER, GC_START_POSITION_HEADER, GC_END_POSITION_HEADER, GC_REF_ALLELE_HEADER,
                          GC_VAR_ALLELE_2_HEADER]

# columns for structural variant annotation
SV_GENEA_HEADER = ['SITE1_GENE', 'GENEA', 'GENE1', 'SITE1_HUGO_SYMBOL']
SV_GENEB_HEADER = ['SITE2_GENE', 'GENEB', 'GENE2', 'SITE2_HUGO_SYMBOL']
SV_TYPE_HEADER = ['SV_CLASS_NAME', 'SV_TYPE', 'CLASS']
SV_TYPES = ['DELETION', 'TRANSLOCATION', 'DUPLICATION', 'INSERTION', 'INVERSION', 'FUSION', 'UNKNOWN']

DESCRIPTION_HEADERS = ['GENE_SUMMARY', 'VARIANT_SUMMARY', 'TUMOR_TYPE_SUMMARY', 'DIAGNOSTIC_SUMMARY',
                       'PROGNOSTIC_SUMMARY', 'MUTATION_EFFECT_DESCRIPTION']

ONCOKB_ANNOTATION_HEADERS_GC = ["ONCOKB_HUGO_SYMBOL", "ONCOKB_PROTEIN_CHANGE", "ONCOKB_CONSEQUENCE"]

UNKNOWN = 'UNKNOWN'


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
    ret = ['Files annotated on ' + date.today().strftime('%m/%d/%Y') + "\nOncoKB API URL: " + oncokb_annotation_api_url]
    try:
        info = requests.get(oncokb_annotation_api_url + "/info", timeout=REQUEST_TIMEOUT).json()
        ret.append(
            '\nOncoKB data version: ' + info['dataVersion']['version'] + ', released on ' + info['dataVersion']['date'])
    except Exception:
        log.error("error when fetch OncoKB info")
    return ''.join(ret)


def validate_oncokb_token():
    if not oncokb_annotation_api_url.startswith(DEFAULT_ONCOKB_URL):
        log.warning(
            "OncoKB base url has been specified by the user that is different from the default www.oncokb.org. The token validation is skipped.")
        return None

    if oncokb_api_bearer_token is None or not oncokb_api_bearer_token:
        log.error("Please specify your OncoKB token")
        exit()

    response = requests.get(oncokb_api_url + "/tokens/" + oncokb_api_bearer_token, timeout=REQUEST_TIMEOUT)
    if response.status_code == 200:
        token = response.json()
        time_stamp = datetime.datetime.strptime(token['expiration'], "%Y-%m-%dT%H:%M:%SZ")
        days_from_expiration = time_stamp - datetime.datetime.now()
        if (days_from_expiration.days < 0):
            log.error(
                "Your OncoKB API token already expired. Please reach out to us to renew your token.")
            exit()
        elif (days_from_expiration.days < 7):
            log.warning(
                "Your OncoKB API token will expire soon, please be on the lookout for an OncoKB email to renew your token. Expire on " + str(
                    time_stamp) + ' UTC')
        else:
            log.info("Your OncoKB API token is valid and will expire on " + str(time_stamp) + ' UTC')
    else:
        try:
            response_json = response.json()
            reason = response_json["title"]
            if response_json["detail"]:
                reason = response_json["detail"]
        except Exception:
            reason = response.reason

        log.error("Error when validating token, " + "reason: %s" % reason)
        exit()


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
        log.error("error when processing %s \n" % url + "reason: %s" % response.reason)
    return hotspots


def requests_retry_session(
        retries=3,
        backoff_factor=0.3,
        status_forcelist=API_REQUEST_RETRY_STATUS_FORCELIST,
        allowed_methods=('GET', 'HEAD'),
        session=None,
):
    session = session or requests.Session()
    retry = Retry(
        total=retries,
        read=retries,
        connect=retries,
        backoff_factor=backoff_factor,
        status_forcelist=status_forcelist,
        allowed_methods=allowed_methods,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount('http://', adapter)
    session.mount('https://', adapter)
    return session


def makeoncokbpostrequest(url, body):
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer %s' % oncokb_api_bearer_token
    }
    return requests_retry_session(allowed_methods=["POST"]).post(url, headers=headers,
                                                                 data=json.dumps(body, default=lambda o: o.__dict__),
                                                                 timeout=REQUEST_TIMEOUT)


def makeoncokbgetrequest(url):
    headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer %s' % oncokb_api_bearer_token
    }
    return requests_retry_session(allowed_methods=["HEAD", "GET"]).get(url, headers=headers, timeout=REQUEST_TIMEOUT)


_3dhotspots = None


def init_3d_hotspots():
    global _3dhotspots
    _3dhotspots = gethotspots(_3dhotspotsbaseurl + "/api/hotspots/3d", None)


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
    threecharactersearch = re.findall(r'[a-zA-Z]{3}\d+', hgvs, flags=re.IGNORECASE)
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
        log.info(
            "Cancer type for the sample should be defined for a more accurate result. \tLine %s" % (row_index))
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

    if selected_query_type is None and has_desired_headers(REQUIRED_QUERY_TYPE_COLUMNS[QueryType.GENOMIC_CHANGE],
                                                           headers):
        selected_query_type = QueryType.GENOMIC_CHANGE

    # default to HGVSp_Short
    if selected_query_type is None:
        selected_query_type = QueryType.HGVSP_SHORT

    # check the file has required columns
    if has_desired_headers(REQUIRED_QUERY_TYPE_COLUMNS[selected_query_type], headers) is False:
        # when it is False, it will never be GENOMIC_CHANGE. For other types, we need to check whether ALTERATION column is available
        if ALTERATION_HEADER not in headers:
            raise Exception(
                "The file does not have required columns "
                + ', '.join(REQUIRED_QUERY_TYPE_COLUMNS[user_input_query_type])
                + " for the query type: "
                + user_input_query_type.value
            )

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


def append_headers(outf, newncols, include_descriptions, genomic_change_annotation):
    oncokb_annotation_headers = get_oncokb_annotation_column_headers(include_descriptions, genomic_change_annotation)
    outf.write("\t".join(oncokb_annotation_headers))
    newncols += len(oncokb_annotation_headers)

    outf.write("\n")
    return newncols


def processalterationevents(eventfile, outfile, previousoutfile, defaultCancerType, cancerTypeMap,
                            annotatehotspots, user_input_query_type, default_reference_genome, include_descriptions):
    if annotatehotspots:
        init_3d_hotspots()
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+', 1000)
    with open(eventfile, DEFAULT_READ_FILE_MODE) as infile:
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

        outf.write("\t")

        query_type = resolve_query_type(user_input_query_type, headers)
        if (query_type == QueryType.HGVSP_SHORT):
            newncols = append_headers(outf, newncols, include_descriptions, False)
            process_alteration(reader, outf, headers, [HGVSP_SHORT_HEADER, ALTERATION_HEADER], ncols, newncols,
                               defaultCancerType,
                               cancerTypeMap, annotatehotspots, default_reference_genome, include_descriptions)

        if (query_type == QueryType.HGVSP):
            newncols = append_headers(outf, newncols, include_descriptions, False)
            process_alteration(reader, outf, headers, [HGVSP_HEADER, ALTERATION_HEADER], ncols, newncols,
                               defaultCancerType,
                               cancerTypeMap, annotatehotspots, default_reference_genome, include_descriptions)

        if (query_type == QueryType.HGVSG):
            newncols = append_headers(outf, newncols, include_descriptions, True)
            process_hvsg(reader, outf, headers, [HGVSG_HEADER, ALTERATION_HEADER], ncols, newncols, defaultCancerType,
                         cancerTypeMap, annotatehotspots, default_reference_genome, include_descriptions)

        if (query_type == QueryType.GENOMIC_CHANGE):
            newncols = append_headers(outf, newncols, include_descriptions, True)
            process_genomic_change(reader, outf, headers, ncols, newncols, defaultCancerType, cancerTypeMap,
                                   annotatehotspots, default_reference_genome, include_descriptions)

    outf.close()


def get_cell_content(row, index, return_empty_string=False):
    if index >= 0 and row[index] != 'NULL' and row[index] != '':
        return row[index]
    elif return_empty_string:
        return ''
    else:
        return None


def get_oncokb_annotation_column_headers(include_descriptions, genomic_change_annotation):
    headers = [ANNOTATED_HEADER]
    if genomic_change_annotation:
        headers.extend(ONCOKB_ANNOTATION_HEADERS_GC)

    headers.extend([GENE_IN_ONCOKB_HEADER,
                    VARIANT_IN_ONCOKB_HEADER,
                    "MUTATION_EFFECT",
                    "MUTATION_EFFECT_CITATIONS",
                    "ONCOGENIC"])

    for level in sorted(levels):
        headers.append(level)
    headers.append("HIGHEST_LEVEL")
    headers.append("HIGHEST_SENSITIVE_LEVEL")
    headers.append("HIGHEST_RESISTANCE_LEVEL")
    headers.append("TX_CITATIONS")

    for dx_level in dxLevels:
        headers.append(dx_level)
    headers.append("HIGHEST_DX_LEVEL")
    headers.append("DX_CITATIONS")

    for px_level in pxLevels:
        headers.append(px_level)
    headers.append("HIGHEST_PX_LEVEL")
    headers.append("PX_CITATIONS")

    if include_descriptions:
        headers.extend(DESCRIPTION_HEADERS)

    return headers


def process_alteration(maffilereader, outf, maf_headers, alteration_column_names, ncols, nannotationcols,
                       defaultCancerType, cancerTypeMap,
                       annotatehotspots, default_reference_genome, include_descriptions):
    ihugo = geIndexOfHeader(maf_headers, HUGO_HEADERS)
    iconsequence = geIndexOfHeader(maf_headers, CONSEQUENCE_HEADERS)
    ihgvs = geIndexOfHeader(maf_headers, alteration_column_names)
    isample = geIndexOfHeader(maf_headers, SAMPLE_HEADERS)
    istart = geIndexOfHeader(maf_headers, PROTEIN_START_HEADERS)
    iend = geIndexOfHeader(maf_headers, PROTEIN_END_HEADERS)
    iproteinpos = geIndexOfHeader(maf_headers, PROTEIN_POSITION_HEADERS)
    icancertype = geIndexOfHeader(maf_headers, CANCER_TYPE_HEADERS)
    ireferencegenome = geIndexOfHeader(maf_headers, REFERENCE_GENOME_HEADERS)

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
        reference_genome = get_reference_genome_from_row(get_cell_content(row, ireferencegenome),
                                                         default_reference_genome)

        hgvs = conversion(hgvs)

        start = get_cell_content(row, istart)

        end = get_cell_content(row, iend)

        if start is None and iproteinpos >= 0 and row[iproteinpos] != "" and row[iproteinpos] != "." and \
                row[iproteinpos] != "-":
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
            annotations = pull_protein_change_info(queries, include_descriptions, annotatehotspots)
            append_annotation_to_file(outf, ncols + nannotationcols, rows, annotations)
            queries = []
            rows = []

    if len(queries) > 0:
        annotations = pull_protein_change_info(queries, include_descriptions, annotatehotspots)
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
    except Exception:
        tumor_seq_allele = ""

    return tumor_seq_allele


def process_genomic_change(maffilereader, outf, maf_headers, ncols, nannotationcols, defaultCancerType, cancerTypeMap,
                           annotatehotspots, default_reference_genome, include_descriptions):
    ichromosome = geIndexOfHeader(maf_headers, [GC_CHROMOSOME_HEADER])
    istart = geIndexOfHeader(maf_headers, [GC_START_POSITION_HEADER])
    iend = geIndexOfHeader(maf_headers, [GC_END_POSITION_HEADER])
    irefallele = geIndexOfHeader(maf_headers, [GC_REF_ALLELE_HEADER])
    ivarallele1 = geIndexOfHeader(maf_headers, [GC_VAR_ALLELE_1_HEADER])
    ivarallele2 = geIndexOfHeader(maf_headers, [GC_VAR_ALLELE_2_HEADER])

    isample = geIndexOfHeader(maf_headers, SAMPLE_HEADERS)
    icancertype = geIndexOfHeader(maf_headers, CANCER_TYPE_HEADERS)
    ireferencegenome = geIndexOfHeader(maf_headers, REFERENCE_GENOME_HEADERS)

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
        reference_genome = get_reference_genome_from_row(get_cell_content(row, ireferencegenome),
                                                         default_reference_genome)

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
            annotations = pull_genomic_change_info(queries, include_descriptions, annotatehotspots)
            append_annotation_to_file(outf, ncols + nannotationcols, rows, annotations)
            queries = []
            rows = []

    if len(queries) > 0:
        annotations = pull_genomic_change_info(queries, include_descriptions, annotatehotspots)
        append_annotation_to_file(outf, ncols + nannotationcols, rows, annotations)


def process_hvsg(maffilereader, outf, maf_headers, alteration_column_names, ncols, nannotationcols, defaultCancerType,
                 cancerTypeMap, annotatehotspots, default_reference_genome, include_descriptions):
    ihgvsg = geIndexOfHeader(maf_headers, alteration_column_names)
    isample = geIndexOfHeader(maf_headers, SAMPLE_HEADERS)
    icancertype = geIndexOfHeader(maf_headers, CANCER_TYPE_HEADERS)
    ireferencegenome = geIndexOfHeader(maf_headers, REFERENCE_GENOME_HEADERS)

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
        reference_genome = get_reference_genome_from_row(get_cell_content(row, ireferencegenome),
                                                         default_reference_genome)

        if hgvsg is None:
            if annotatehotspots:
                default_cols = [['', '', 'False']]
            else:
                default_cols = [['False']]
            append_annotation_to_file(outf, ncols + nannotationcols, [row],
                                      default_cols)
        else:
            query = HGVSgQuery(hgvsg, cancertype, reference_genome)
            queries.append(query)
            rows.append(row)

        if len(queries) == POST_QUERIES_THRESHOLD_GC_HGVSG:
            annotations = pull_hgvsg_info(queries, include_descriptions, annotatehotspots)
            append_annotation_to_file(outf, ncols + nannotationcols, rows, annotations)
            queries = []
            rows = []

    if len(queries) > 0:
        annotations = pull_hgvsg_info(queries, include_descriptions, annotatehotspots)
        append_annotation_to_file(outf, ncols + nannotationcols, rows, annotations)


def getgenesfromfusion(fusion, nameregex=None):
    GENES_REGEX = r"([A-Za-z\d]+-[A-Za-z\d]+)" if nameregex is None else nameregex
    searchresult = re.search(GENES_REGEX, fusion, flags=re.IGNORECASE)
    geneA = None
    geneB = None
    if searchresult:
        parts = searchresult.group(1).split("-")
        geneA = parts[0]
        geneB = geneA
        if len(parts) > 1 and parts[1] != "intragenic":
            geneB = parts[1]
    else:
        geneA = geneB = fusion
    return geneA, geneB


def process_fusion(svdata, outfile, previousoutfile, defaultCancerType, cancerTypeMap, nameregex, include_descriptions):
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+')
    with open(svdata, DEFAULT_READ_FILE_MODE) as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)

        ncols = headers["length"]

        if ncols == 0:
            return

        outf.write(headers['^-$'])
        oncokb_annotation_headers = get_oncokb_annotation_column_headers(include_descriptions, False)
        outf.write("\t")
        outf.write("\t".join(oncokb_annotation_headers))
        outf.write("\n")

        newcols = ncols + len(oncokb_annotation_headers)

        igeneA = geIndexOfHeader(headers, SV_GENEA_HEADER)
        igeneB = geIndexOfHeader(headers, SV_GENEB_HEADER)
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

            geneA = None
            geneB = None
            if igeneA >= 0:
                geneA = row[igeneA]
            if igeneB >= 0:
                geneB = row[igeneB]
            if igeneA < 0 and igeneB < 0 and ifusion >= 0:
                fusion = row[ifusion]
                geneA, geneB = getgenesfromfusion(fusion, nameregex)

            cancertype = get_tumor_type_from_row(row, i, defaultCancerType, icancertype, cancerTypeMap, sample)

            queries.append(StructuralVariantQuery(geneA, geneB, 'FUSION', cancertype))
            rows.append(row)

            if len(queries) == POST_QUERIES_THRESHOLD:
                annotations = pull_structural_variant_info(queries, include_descriptions)
                append_annotation_to_file(outf, newcols, rows, annotations)
                queries = []
                rows = []

        if len(queries) > 0:
            annotations = pull_structural_variant_info(queries, include_descriptions)
            append_annotation_to_file(outf, newcols, rows, annotations)
    outf.close()


def process_sv(svdata, outfile, previousoutfile, defaultCancerType, cancerTypeMap, include_descriptions):
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)
    outf = open(outfile, 'w+')
    with open(svdata, DEFAULT_READ_FILE_MODE) as infile:
        reader = csv.reader(infile, delimiter='\t')

        headers = readheaders(reader)

        ncols = headers["length"]

        if ncols == 0:
            return

        outf.write(headers['^-$'])
        oncokb_annotation_headers = get_oncokb_annotation_column_headers(include_descriptions, False)
        outf.write("\t")
        outf.write("\t".join(oncokb_annotation_headers))
        outf.write("\n")

        newcols = ncols + len(oncokb_annotation_headers)

        igeneA = geIndexOfHeader(headers, SV_GENEA_HEADER)
        igeneB = geIndexOfHeader(headers, SV_GENEB_HEADER)
        isvtype = geIndexOfHeader(headers, SV_TYPE_HEADER)
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

            if igeneA < 0 or igeneB < 0:
                log.warning("Please specify two genes")
                continue

            svtype = None
            if isvtype >= 0:
                svtype = row[isvtype].upper()
                if svtype not in SV_TYPES:
                    svtype = None
            if svtype is None:
                svtype = UNKNOWN

            cancertype = get_tumor_type_from_row(row, i, defaultCancerType, icancertype, cancerTypeMap, sample)

            sv_query = StructuralVariantQuery(row[igeneA], row[igeneB], svtype, cancertype)
            queries.append(sv_query)
            rows.append(row)

            if len(queries) == POST_QUERIES_THRESHOLD:
                annotations = pull_structural_variant_info(queries, include_descriptions)
                append_annotation_to_file(outf, newcols, rows, annotations)
                queries = []
                rows = []

        if len(queries) > 0:
            annotations = pull_structural_variant_info(queries, include_descriptions)
            append_annotation_to_file(outf, newcols, rows, annotations)
    outf.close()


def get_cna(cell_value, annotate_gain_loss=False):
    cna = None
    if cell_value is not None and cell_value != '':
        if cell_value in GISTIC_CNA_MAP:
            cna = GISTIC_CNA_MAP[cell_value]
        else:
            for default_cna in CNAS:
                if cell_value.upper() == default_cna.upper():
                    cna = default_cna
    if not annotate_gain_loss and cna is not None and cna.upper() in [CNA_GAIN_TXT.upper(), CNA_LOSS_TXT.upper()]:
        cna = None
    return cna


def process_gistic_data(outf, gistic_data_file, defaultCancerType, cancerTypeMap, annotate_gain_loss,
                        include_descriptions):
    with open(gistic_data_file, DEFAULT_READ_FILE_MODE) as infile:
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
                    cna_type = get_cna(row[headers[rawsample]], annotate_gain_loss)
                    if cna_type is not None:
                        cancer_type = defaultCancerType
                        sample = rawsample

                        if sampleidsfilter and sample not in sampleidsfilter:
                            continue

                        if sample in cancerTypeMap:
                            cancer_type = cancerTypeMap[sample]

                        rows.append([sample, cancer_type, hugo, cna_type])
                        queries.append(CNAQuery(hugo, cna_type, cancer_type))

        headers = ['SAMPLE_ID', 'CANCER_TYPE', 'HUGO_SYMBOL', 'ALTERATION'] + get_oncokb_annotation_column_headers(
            include_descriptions, False)
        outf.write('\t'.join(headers))
        outf.write('\n')
        return headers, rows, queries


def process_individual_cna_file(outf, cna_data_file, defaultCancerType, cancerTypeMap, annotate_gain_loss,
                                include_descriptions):
    with open(cna_data_file, DEFAULT_READ_FILE_MODE) as infile:
        reader = csv.reader(infile, delimiter='\t')
        headers = readheaders(reader)
        row_headers = headers['^-$'].split('\t') + get_oncokb_annotation_column_headers(include_descriptions, False)

        i = 0
        rows = []
        queries = []

        outf.write('\t'.join(row_headers))
        outf.write('\n')

        for row in reader:
            i = i + 1
            isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
            ihugo = geIndexOfHeader(headers, HUGO_HEADERS)
            icancertype = geIndexOfHeader(headers, CANCER_TYPE_HEADERS)
            icna = geIndexOfHeader(headers, CNA_HEADERS)

            hugo = row[ihugo] if ihugo >= 0 else None
            cna_type = get_cna(row[icna], annotate_gain_loss)
            sample = row[isample] if isample >= 0 else None
            cancer_type = get_tumor_type_from_row(row, i, defaultCancerType, icancertype, cancerTypeMap, sample)

            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            if hugo and cna_type:
                rows.append(row)
                queries.append(CNAQuery(hugo, cna_type, cancer_type))
            else:
                outf.write('\t'.join(row))
                outf.write('\n')
                if not hugo:
                    log.warning("Gene is not specified for row " + str(row))
                if not cna_type:
                    log.warning("CNA is not specified for row " + str(row))
        return row_headers, rows, queries


def process_cna_data(cnafile, outfile, previousoutfile, defaultCancerType, cancerTypeMap, include_descriptions,
                     annotate_gain_loss=False,
                     cna_format=CNA_FILE_FORMAT_GISTIC):
    if os.path.isfile(previousoutfile):
        cacheannotated(previousoutfile, defaultCancerType, cancerTypeMap)

    if not cna_format or cna_format not in CND_FILE_FORMAT:
        log.error('The CNA file format is not supported, only gistic or individual can be used ')
        return

    outf = open(outfile, 'w+', 1000)

    headers = []
    rows = []
    queries = []
    if cna_format == CNA_FILE_FORMAT_GISTIC:
        headers, rows, queries = process_gistic_data(outf, cnafile, defaultCancerType, cancerTypeMap,
                                                     annotate_gain_loss, include_descriptions)
    else:
        headers, rows, queries = process_individual_cna_file(outf, cnafile, defaultCancerType, cancerTypeMap,
                                                             annotate_gain_loss, include_descriptions)

    ncols = len(headers)

    i = 0
    while len(rows) > 0:
        i += POST_QUERIES_THRESHOLD
        log.info(i)
        rows_sec, rows = rows[:POST_QUERIES_THRESHOLD], rows[POST_QUERIES_THRESHOLD:]
        queries_sec, queries = queries[:POST_QUERIES_THRESHOLD], queries[POST_QUERIES_THRESHOLD:]
        annotations = pull_cna_info(queries_sec, include_descriptions)
        append_annotation_to_file(outf, ncols, rows_sec, annotations)

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


def process_clinical_data(annotatedmutfiles, clinicalfile, outfile):
    samplelevels = {}
    sampledxlevels = {}
    samplepxlevels = {}
    sampleleveltreatments = {}
    sampledrivers = {}
    sample_resistance = {}
    samplemutationswithdiagnosis = {}
    samplemutationswithprognosis = {}
    sample_tx_sensitive_count = {}
    sample_tx_resistance_count = {}
    samplealterationcount = {}
    for annotatedmutfile in annotatedmutfiles:
        with open(annotatedmutfile, DEFAULT_READ_FILE_MODE) as mutfile:
            reader = csv.reader(mutfile, delimiter='\t')
            headers = readheaders(reader)

            ncols = headers["length"]

            if ncols == 0:
                return

            igeneA = geIndexOfHeader(headers, SV_GENEA_HEADER)  # fusion
            igeneB = geIndexOfHeader(headers, SV_GENEB_HEADER)  # fusion
            ifusion = geIndexOfHeader(headers, ['FUSION'])

            ihugo = geIndexOfHeader(headers, HUGO_HEADERS)
            ihgvs = geIndexOfHeader(headers, HGVS_HEADERS)
            isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
            ioncogenic = headers['ONCOGENIC']

            isfusion = (igeneA != -1 and igeneB != -1) or ifusion != -1
            ismutorcna = ihugo != -1 and ihgvs != -1

            if not isfusion and not ismutorcna:
                log.error("file " + annotatedmutfile + " missing proper header")
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
                    sample_resistance[sample] = []
                    sample_tx_sensitive_count[sample] = {}
                    sample_tx_resistance_count[sample] = {}

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
                geneA = row[igeneA]
                geneB = row[igeneB]

                variant = "NA"
                if ismutorcna:
                    variant = hugo + " " + alteration
                elif isfusion:
                    if ifusion != -1:
                        variant = row[ifusion]
                    else:
                        if geneA == geneB:
                            variant = geneA + " intragenic deletion"
                        else:
                            variant = geneA + "-" + geneB + " fusion"

                if oncogenic == "oncogenic" or oncogenic == "likely oncogenic" or oncogenic == "predicted oncogenic":
                    sampledrivers[sample].append(variant)
                if oncogenic == "resistance":
                    sample_resistance[sample].append(variant)

                for level in levels:
                    il = geIndexOfHeader(headers, [level])
                    if 0 <= il < len(row) and row[il] != '':
                        if level not in samplelevels[sample]:
                            samplelevels[sample][level] = []
                            sampleleveltreatments[sample][level] = []
                        samplelevels[sample][level].append(row[il] + "(" + variant + ")")
                        sampleleveltreatments[sample][level].extend(row[il].split(","))

                        if level.startswith('LEVEL_R'):
                            sample_tx_resistance_count[sample][variant] = True
                        else:
                            sample_tx_sensitive_count[sample][variant] = True

                for dx_level in dxLevels:
                    il = geIndexOfHeader(headers, [dx_level])
                    if 0 <= il < len(row) and row[il] != '':
                        if dx_level not in samplelevels[sample]:
                            samplelevels[sample][dx_level] = []
                        samplelevels[sample][dx_level].append(row[il] + "(" + variant + ")")

                for px_level in pxLevels:
                    il = geIndexOfHeader(headers, [px_level])
                    if 0 <= il < len(row) and row[il] != '':
                        if px_level not in samplelevels[sample]:
                            samplelevels[sample][px_level] = []
                        samplelevels[sample][px_level].append(row[il] + "(" + variant + ")")

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

    # export to annotated file
    with open(clinicalfile, DEFAULT_READ_FILE_MODE) as clinfile:
        reader = csv.reader(clinfile, delimiter='\t')
        headers = readheaders(reader)
        outf.write(headers['^-$'])
        for level in sorted(levels):
            outf.write('\t' + level)
        outf.write('\tHIGHEST_LEVEL')
        outf.write('\tHIGHEST_SENSITIVE_LEVEL')
        outf.write('\tHIGHEST_RESISTANCE_LEVEL')
        for dx_level in dxLevels:
            outf.write('\t' + dx_level)
        outf.write('\tHIGHEST_DX_LEVEL')
        for px_level in pxLevels:
            outf.write('\t' + px_level)
        outf.write('\tHIGHEST_PX_LEVEL')
        outf.write(
            '\tONCOGENIC_MUTATIONS\t#ONCOGENIC_MUTATIONS\tRESISTANCE_MUTATIONS\t#RESISTANCE_MUTATIONS\t#MUTATIONS_WITH_SENSITIVE_THERAPEUTIC_IMPLICATIONS\t#MUTATIONS_WITH_RESISTANCE_THERAPEUTIC_IMPLICATIONS\t#MUTATIONS_WITH_DIAGNOSTIC_IMPLICATIONS\t#MUTATIONS_WITH_PROGNOSTIC_IMPLICATIONS\t#MUTATIONS\n')
        isample = headers['SAMPLE_ID']

        for row in reader:
            sample = row[isample]

            if sampleidsfilter and sample not in sampleidsfilter:
                continue

            outf.write('\t'.join(row))

            for level in sorted(levels):
                outf.write('\t')
                if sample in samplelevels and level in samplelevels[sample]:
                    outf.write(";".join(samplelevels[sample][level]))

            highestlevel = ''
            highest_sensitive_level = ''
            highest_resistance_level = ''
            highestdxlevel = ''
            highestpxlevel = ''
            if sample in sampleleveltreatments:
                highestlevel = get_highest_tx_level(sampleleveltreatments[sample])
                highest_sensitive_level = get_highest_tx_level(sampleleveltreatments[sample], TX_TYPE_SENSITIVE)
                highest_resistance_level = get_highest_tx_level(sampleleveltreatments[sample], TX_TYPE_RESISTANCE)
            if sample in sampledxlevels:
                highestdxlevel = get_highest_dxpx_level(dxLevels, sampledxlevels[sample])
            if sample in samplepxlevels:
                highestpxlevel = get_highest_dxpx_level(pxLevels, samplepxlevels[sample])
            # if highestlevel == '':
            #     if sample in sampledrivers and len(sampledrivers[sample])>0:
            #         highestlevel = 'Oncogenic, no level'
            #     else:
            #         highestlevel = "VUS"
            outf.write('\t' + highestlevel)
            outf.write('\t' + highest_sensitive_level)
            outf.write('\t' + highest_resistance_level)

            for dx_level in dxLevels:
                outf.write('\t')
                if sample in samplelevels and dx_level in samplelevels[sample]:
                    outf.write(";".join(samplelevels[sample][dx_level]))

            outf.write('\t' + highestdxlevel)

            for px_level in pxLevels:
                outf.write('\t')
                if sample in samplelevels and px_level in samplelevels[sample]:
                    outf.write(";".join(samplelevels[sample][px_level]))
            outf.write('\t' + highestpxlevel)

            tx_sensitive_count = 0
            tx_resistance_count = 0
            if sample in sample_tx_sensitive_count:
                tx_sensitive_count = len(sample_tx_sensitive_count[sample].keys())
            if sample in sample_tx_resistance_count:
                tx_resistance_count = len(sample_tx_resistance_count[sample].keys())

            alterationcount = 0
            if sample in samplealterationcount:
                alterationcount = samplealterationcount[sample]

            drivercount = 0
            resistance_count = 0
            diagnosiscount = 0
            prognosiscount = 0
            drivermutations = ""
            resistance_mutations = ""
            if sample in sampledrivers:
                drivercount = len(sampledrivers[sample])
                drivermutations = ";".join(sampledrivers[sample])
            if sample in sample_resistance:
                resistance_count = len(sample_resistance[sample])
                resistance_mutations = ";".join(sample_resistance[sample])
            if sample in samplemutationswithdiagnosis:
                diagnosiscount = len(samplemutationswithdiagnosis[sample])
            if sample in samplemutationswithprognosis:
                prognosiscount = len(samplemutationswithprognosis[sample])

            outf.write('\t' + drivermutations)
            outf.write('\t' + str(drivercount))
            outf.write('\t' + str(resistance_mutations))
            outf.write('\t' + str(resistance_count))
            outf.write('\t' + str(tx_sensitive_count))
            outf.write('\t' + str(tx_resistance_count))
            outf.write('\t' + str(diagnosiscount))
            outf.write('\t' + str(prognosiscount))
            outf.write('\t' + str(alterationcount))

            outf.write('\n')

    outf.close()


oncokbcache = {}


def cacheannotated(annotatedfile, defaultCancerType, cancerTypeMap):
    with open(annotatedfile, DEFAULT_READ_FILE_MODE) as infile:
        try:
            reader = csv.reader(infile, delimiter='\t')
            headers = readheaders(reader)

            ihugo = geIndexOfHeader(headers, HUGO_HEADERS)
            ihgvs = geIndexOfHeader(headers, HGVS_HEADERS)
            isample = geIndexOfHeader(headers, SAMPLE_HEADERS)
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
                    for level in levels:
                        il = headers[level]
                        if il < len(row):
                            oncokbcache[key][level] = row[il].split(',')
                        else:
                            oncokbcache[key][level] = []
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

    def __str__(self):
        return self.hugoSymbol


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

    def __repr__(self):
        return ",".join([self.gene.hugoSymbol, self.alteration, self.tumorType, self.consequence, self.proteinStart,
                         self.proteinEnd, self.referenceGenome])


class HGVSgQuery:
    def __init__(self, hgvsg, cancertype, reference_genome=None):
        self.hgvsg = hgvsg
        self.tumorType = cancertype
        if reference_genome is not None:
            self.referenceGenome = reference_genome.value

    def __repr__(self):
        return ",".join([self.hgvsg, self.tumorType, self.referenceGenome])


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
        if chromosome is not None:
            chromosome = chromosome.strip()
            if chromosome.startswith('chr'):
                chromosome = chromosome[3:]
        self.genomicLocation = ','.join([chromosome, start, end, ref_allele, var_allele])
        self.tumorType = cancertype
        if reference_genome is not None:
            self.referenceGenome = reference_genome.value

    def __repr__(self):
        return " ".join([self.genomicLocation, self.tumorType])


class CNAQuery:
    def __init__(self, hugo, cnatype, cancertype):
        self.gene = Gene(hugo)
        self.copyNameAlterationType = cnatype.upper()
        self.tumorType = cancertype

    def __repr__(self):
        return ",".join([self.gene.hugoSymbol, self.copyNameAlterationType, self.tumorType])


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

    def __repr__(self):
        return ",".join(
            [self.geneA.hugoSymbol, self.geneB.hugoSymbol, str(self.functionalFusion), self.structuralVariantType,
             self.tumorType])


def pull_protein_change_info(queries, include_descriptions, annotate_hotspot):
    url = oncokb_annotation_api_url + '/annotate/mutations/byProteinChange'
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
            if hasattr(query,
                       'proteinStart') and query.proteinStart and query.proteinStart != '\\N' and query.proteinStart != 'NULL' and query.proteinStart != '':
                geturl += '&proteinStart=' + str(query.proteinStart)
            if hasattr(query,
                       'proteinEnd') and query.proteinEnd and query.proteinEnd != '\\N' and query.proteinEnd != 'NULL' and query.proteinEnd != '':
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
        processed_annotation.append(
            process_oncokb_annotation(query_annotation, include_descriptions, False, annotate_hotspot))
    return processed_annotation


def pull_hgvsg_info(queries, include_descriptions, annotate_hotspot):
    url = oncokb_annotation_api_url + '/annotate/mutations/byHGVSg'
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
        processed_annotation.append(
            process_oncokb_annotation(query_annotation, include_descriptions, True, annotate_hotspot))
    return processed_annotation


def pull_genomic_change_info(queries, include_descriptions, annotate_hotspot):
    url = oncokb_annotation_api_url + '/annotate/mutations/byGenomicChange'
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
        processed_annotation.append(
            process_oncokb_annotation(query_annotation, include_descriptions, True, annotate_hotspot))
    return processed_annotation


def pull_cna_info(queries, include_descriptions):
    url = oncokb_annotation_api_url + '/annotate/copyNumberAlterations'

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
        processed_annotation.append(
            process_oncokb_annotation(query_annotation, include_descriptions, False, annotate_hotspot=False))
    return processed_annotation


def pull_structural_variant_info(queries, include_descriptions):
    url = oncokb_annotation_api_url + '/annotate/structuralVariants'

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
            geturl += '&isFunctionalFusion=' + str(query.functionalFusion).upper() if type(
                query.functionalFusion) is bool else query.functionalFusion
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
        processed_annotation.append(
            process_oncokb_annotation(query_annotation, include_descriptions, False, annotate_hotspot=False))
    return processed_annotation


def process_oncokb_annotation(annotation, include_descriptions, genomic_change_annotation, annotate_hotspot):
    if annotation is None:
        return ['False']

    oncokbdata = {}
    for level in levels:
        oncokbdata[level] = []
    for dx_level in dxLevels:
        oncokbdata[dx_level] = []
    for px_level in pxLevels:
        oncokbdata[px_level] = []

    oncokbdata[GENE_IN_ONCOKB_HEADER] = GENE_IN_ONCOKB_DEFAULT
    oncokbdata[VARIANT_IN_ONCOKB_HEADER] = VARIANT_IN_ONCOKB_DEFAULT
    oncokbdata['mutation_effect'] = ""
    oncokbdata['mutation_effect_citations'] = []
    oncokbdata['mutation_effect_description'] = ""
    oncokbdata['citations'] = []
    oncokbdata['oncogenic'] = ""
    oncokbdata['tx_citations'] = []
    oncokbdata['dx_citations'] = []
    oncokbdata['px_citations'] = []

    try:
        # oncogenic
        oncokbdata[GENE_IN_ONCOKB_HEADER] = GENE_IN_ONCOKB_DEFAULT if annotation['geneExist'] is None else str(
            annotation['geneExist'])
        oncokbdata[VARIANT_IN_ONCOKB_HEADER] = VARIANT_IN_ONCOKB_DEFAULT if annotation['variantExist'] is None else str(
            annotation['variantExist'])

        # oncogenic
        oncokbdata['oncogenic'] = annotation['oncogenic']

        # if not evidences['geneExist'] or (not evidences['variantExist'] and not evidences['alleleExist']):
        #     return ''

        # mutation effect
        if (annotation['mutationEffect'] is not None):
            oncokbdata['mutation_effect'] = annotation['mutationEffect']['knownEffect']
            oncokbdata['mutation_effect_description'] = annotation['mutationEffect']['description']
            oncokbdata['mutation_effect_citations'] = appendoncokbcitations(oncokbdata['mutation_effect_citations'],
                                                                            annotation['mutationEffect']['citations'][
                                                                                'pmids'],
                                                                            annotation['mutationEffect']['citations'][
                                                                                'abstracts'])

        # oncogenic
        oncokbdata['oncogenic'] = annotation['oncogenic']

        # get treatment
        for treatment in annotation['treatments']:
            level = treatment['level']

            if level not in levels:
                log.info("%s is ignored" % level)
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
    except Exception:
        log.error("error when processing %s " % annotation)

    ret = []
    if annotate_hotspot:
        if annotation['hotspot']:
            ret.append('Y')
        else:
            ret.append('')

        _3dhotspot = pull3dhotspots(annotation['query']['hugoSymbol'], annotation['query']['consequence'],
                                    annotation['query']['proteinStart'], annotation['query']['proteinEnd'])
        ret.append(_3dhotspot)

    ret.append('True')

    if genomic_change_annotation:
        query_hugo_symbol = annotation['query']['hugoSymbol']
        ret.append('' if query_hugo_symbol is None else query_hugo_symbol)

        query_alteration = annotation['query']['alteration']
        ret.append('' if query_alteration is None else query_alteration)

        query_consequence = annotation['query']['consequence']
        ret.append('' if query_consequence is None else query_consequence)

    ret.append(oncokbdata[GENE_IN_ONCOKB_HEADER])
    ret.append(oncokbdata[VARIANT_IN_ONCOKB_HEADER])
    ret.append(oncokbdata['mutation_effect'])
    ret.append(';'.join(oncokbdata['mutation_effect_citations']))
    ret.append(oncokbdata['oncogenic'])
    for level in sorted(levels):
        ret.append(','.join(oncokbdata[level]))
    ret.append(get_highest_tx_level(oncokbdata))
    ret.append(get_highest_tx_level(oncokbdata, TX_TYPE_SENSITIVE))
    ret.append(get_highest_tx_level(oncokbdata, TX_TYPE_RESISTANCE))
    ret.append(';'.join(oncokbdata['tx_citations']))

    for dx_level in dxLevels:
        ret.append(','.join(oncokbdata[dx_level]))
    ret.append(get_highest_dxpx_level(dxLevels, [oncokbdata['highestDiagnosticImplicationLevel']]))
    ret.append(';'.join(oncokbdata['dx_citations']))

    for px_level in pxLevels:
        ret.append(','.join(oncokbdata[px_level]))
    ret.append(get_highest_dxpx_level(pxLevels, [oncokbdata['highestPrognosticImplicationLevel']]))
    ret.append(';'.join(oncokbdata['px_citations']))

    if include_descriptions:
        ret.append(annotation['geneSummary'])
        ret.append(annotation['variantSummary'])
        ret.append(annotation['tumorTypeSummary'])
        ret.append(annotation['diagnosticSummary'])
        ret.append(annotation['prognosticSummary'])
        ret.append(oncokbdata['mutation_effect_description'])

    return ret


def get_highest_tx_level(oncokb_data, tx_type=None):
    target_levels = levels
    if tx_type is not None and tx_type:
        if tx_type.lower() == TX_TYPE_SENSITIVE:
            target_levels = sensitive_levels
        elif tx_type.lower() == TX_TYPE_RESISTANCE:
            target_levels = resistance_levels
    for level in target_levels:
        if level in oncokb_data and oncokb_data[level] is not None and len(oncokb_data[level]) > 0:
            return level
    return ""


def get_highest_dxpx_level(dxpx_levels, oncokbdata):
    for level in dxpx_levels:
        if level not in oncokbdata:
            continue
        return level
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
    with open(clinicalFile, DEFAULT_READ_FILE_MODE) as infile:
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
                h = h.strip()
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
