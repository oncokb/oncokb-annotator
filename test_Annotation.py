#!/usr/bin/python

from AnnotatorCore import *
import os

setoncokbapitoken(os.environ("ONCOKB_API_TOKEN"))

log.info('test-----------', os.getenv("ONCOKB_API_TOKEN"), '------')

VARIANT_EXISTS_INDEX = 1
MUTATION_EFFECT_INDEX = 2
ONCOGENIC_INDEX = 3
HIGHEST_LEVEL_INDEX = 12
UNKNOWN = 'Unknown'


def get_annotation_from_string(content_str):
    return [] if content_str is None else content_str.split('\t')


def test_check_atypical_alts():
    queries = [
        ProteinChangeQuery('Other Biomarkers', 'MSI-H', 'Colorectal Cancer'),
        ProteinChangeQuery('Other Biomarkers', 'MSI-H', 'Leukemia')
    ]

    annotations = pull_mutation_info(queries)
    assert len(annotations) == 2

    annotation = get_annotation_from_string(annotations[0])
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = get_annotation_from_string(annotations[1])
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''


def test_check_fusions():
    queries = [
        StructuralVariantQuery('ALK', 'EML4', 'FUSION', 'NSCLC'),
        StructuralVariantQuery('ALK', 'EML4', 'FUSION', 'Melanoma'),
    ]

    annotations = pull_structural_variant_info(queries)
    assert len(annotations) == 2

    annotation = get_annotation_from_string(annotations[0])
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = get_annotation_from_string(annotations[1])
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_3B'


def test_cna():
    queries = [
        CNAQuery('BRCA2', 'DELETION', 'Ovarian Cancer'),
        CNAQuery('ERBB2', 'Amplification', 'Breast Cancer'),
        CNAQuery('ERBB2', 'Amplification', 'Colorectal Cancer'),
    ]

    annotations = pull_cna_info(queries)
    assert len(annotations) == 3

    annotation = get_annotation_from_string(annotations[0])
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Loss-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = get_annotation_from_string(annotations[1])
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = get_annotation_from_string(annotations[2])
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_2'
