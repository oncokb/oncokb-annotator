#!/usr/bin/python

import argparse
from AnnotatorCore import *
import logging

logging.basicConfig(level=logging.INFO)
log = logging.getLogger('MafAnnotator')


def main(argv):
    if argv.help:
        log.info(
            '\n'
            'MafAnnotator.py -i <input MAF file> -o <output MAF file> [-p previous results] [-c <input clinical file>] '
            '[-s sample list filter] [-t <default tumor type>] [-u oncokb-base-url] [-b oncokb api bear token] [-a] [-q query type]\n'
            'Essential MAF columns (case insensitive):\n'
            '    HUGO_SYMBOL: Hugo gene symbol\n'
            '    VARIANT_CLASSIFICATION: Translational effect of variant allele\n'
            '    TUMOR_SAMPLE_BARCODE: sample ID\n'
            '    HGVSP_SHORT: protein change in HGVSP format\n'
            '    PROTEIN_START: protein start\n'
            '    PROTEIN_END: protein end\n'
            '    PROTEIN_POSITION: can be used instead of PROTEIN_START and PROTEIN_END (in the output of vcf2map)\n'
            'Essential clinical columns:\n'
            '    SAMPLE_ID: sample ID\n'
            '    ONCOTREE_CODE: tumor type code from oncotree (oncotree.mskcc.org)\n'
            'Cancer type will be assigned based on the following priority:\n'
            '    1) ONCOTREE_CODE in clinical data file\n'
            '    2) ONCOTREE_CODE exist in MAF\n'
            '    3) default tumor type (-t)\n'
            'Query type only allows the following values (case-insensitive):\n'
            '    - HGVSp_Short \n'
            '      It reads from column HGVSp_Short or Alteration\n'
            '    - HGVSp\n'
            '      It reads from column HGVSp or Alteration\n'
            '    - HGVSg\n'
            '      It reads from column HGVSg or Alteration\n'
            '    - Genomic_Change\n'
            '      It reads from columns Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele1 and Tumor_Seq_Allele2  \n'
            'Default OncoKB base url is https://www.oncokb.org.\n'
            'Use -a to annotate mutational hotspots\n'
        )
        sys.exit()
    if argv.input_file == '' or argv.output_file == '' or argv.oncokb_api_bearer_token == '':
        log.info('For help: python MafAnnotator.py -h')
        sys.exit(2)

    if argv.sample_ids_filter:
        setsampleidsfileterfile(argv.sample_ids_filter)
    if argv.cancer_hotspots_base_url:
        setcancerhotspotsbaseurl(argv.cancer_hotspots_base_url)
    if argv.oncokb_api_url:
        setoncokbbaseurl(argv.oncokb_api_url)
    setoncokbapitoken(argv.oncokb_api_bearer_token)
    getcuratedgenes()

    cancertypemap = {}
    if argv.input_clinical_file:
        readCancerTypes(argv.input_clinical_file, cancertypemap)

    log.info('annotating %s ...' % argv.input_file)

    user_input_query_type = None
    if argv.query_type is not None:
        try:
            user_input_query_type = QueryType[argv.query_type.upper()]
        except KeyError:
            # if not isinstance(argv.query_type.upper(), QueryType):
            print(
                'Query type is not acceptable. Only the following allows(case insensitive): HGVSp_Short, HGVSp, HGVSg, Genomic_Change')
            raise

    processalterationevents(argv.input_file, argv.output_file, argv.previous_result_file, argv.default_cancer_type,
                            cancertypemap, True, argv.annotate_hotspots, user_input_query_type)

    log.info('done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-h', dest='help', action="store_true", default=False)
    parser.add_argument('-i', dest='input_file', default='', type=str)
    parser.add_argument('-o', dest='output_file', default='', type=str)
    parser.add_argument('-p', dest='previous_result_file', default='', type=str)
    parser.add_argument('-c', dest='input_clinical_file', default='', type=str)
    parser.add_argument('-s', dest='sample_ids_filter', default='', type=str)
    parser.add_argument('-t', dest='default_cancer_type', default='', type=str)
    parser.add_argument('-u', dest='oncokb_api_url', default='', type=str)
    parser.add_argument('-a', dest='annotate_hotspots', action="store_true", default=False)
    parser.add_argument('-v', dest='cancer_hotspots_base_url', default='', type=str)
    parser.add_argument('-b', dest='oncokb_api_bearer_token', default='', type=str)
    parser.add_argument('-q', dest='query_type', default=None, type=str)
    parser.set_defaults(func=main)

    args = parser.parse_args()
    args.func(args)
