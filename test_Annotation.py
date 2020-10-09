#!/usr/bin/python
import pytest

from AnnotatorCore import *
import os

ONCOKB_API_TOKEN = os.environ["ONCOKB_API_TOKEN"]
setoncokbapitoken(ONCOKB_API_TOKEN)

log.info('test-----------', os.environ["ONCOKB_API_TOKEN"], '------')

VARIANT_EXISTS_INDEX = 1
MUTATION_EFFECT_INDEX = 2
ONCOGENIC_INDEX = 3
HIGHEST_LEVEL_INDEX = 12
UNKNOWN = 'Unknown'

def fake_gene_one_query_suite(annotations):
    assert len(annotations) == 1

    annotation = annotations[0]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == ''
    assert annotation[HIGHEST_LEVEL_INDEX] == ''


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_protein_change():
    queries = [
        ProteinChangeQuery('BRAF', 'V600E', 'Colorectal Cancer')
    ]

    annotations = pull_protein_change_info(queries)
    assert len(annotations) == 1

    annotation = annotations[0]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_fake_gene_protein_change():
    queries = [
        ProteinChangeQuery('test1', 'V600E', 'Ovarian Cancer')
    ]

    annotations = pull_protein_change_info(queries)
    fake_gene_one_query_suite(annotations)


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_atypical_alts():
    queries = [
        ProteinChangeQuery('Other Biomarkers', 'MSI-H', 'Colorectal Cancer'),
        ProteinChangeQuery('Other Biomarkers', 'MSI-H', 'Leukemia'),
        ProteinChangeQuery('TERT', 'Promoter Mutation', 'Bladder Cancer'),
        ProteinChangeQuery('TERT', 'Promoter Mutation', 'Bladder Cancer', '5\'Flank')
    ]

    annotations = pull_protein_change_info(queries)
    assert len(annotations) == 4

    annotation = annotations[0]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

    annotation = annotations[2]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

    annotation_dup = annotations[3]
    assert len(annotation_dup) == 14
    assert annotation == annotation_dup


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_hgvsg():
    queries = [
        HGVSgQuery('12:g.25398285C>A', 'LUAD'),
        HGVSgQuery('12:g.25398285_25398286delinsAG', 'LUAD'),
        HGVSgQuery('5:g.1295228G>A', 'LUAD'),
    ]

    annotations = pull_hgvsg_info(queries)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_3A'

    annotation = annotations[1]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_3A'

    annotation = annotations[2]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_genomic_change():
    queries = [
        GenomicChangeQuery('12', '25398285', '25398285', 'C', 'A', 'LUAD'),
        GenomicChangeQuery('12', '25398285', '25398286', 'CA', 'AG', 'LUAD'),
        GenomicChangeQuery('5', '1295228', '1295228', 'G', 'A', 'LUAD'),
    ]

    annotations = pull_genomic_change_info(queries)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_3A'

    annotation = annotations[1]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_3A'

    annotation = annotations[2]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_fusions():
    queries = [
        StructuralVariantQuery('ALK', 'EML4', 'FUSION', 'NSCLC'),
        StructuralVariantQuery('ALK', 'EML4', 'FUSION', 'Melanoma'),
    ]

    annotations = pull_structural_variant_info(queries)
    assert len(annotations) == 2

    annotation = annotations[0]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_3B'


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_fake_fusion_gene():
    queries = [
        StructuralVariantQuery('test1', 'test2', 'FUSION', 'NSCLC'),
    ]

    annotations = pull_structural_variant_info(queries)
    fake_gene_one_query_suite(annotations)


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_cna():
    queries = [
        CNAQuery('BRCA2', 'DELETION', 'Ovarian Cancer'),
        CNAQuery('ERBB2', 'Amplification', 'Breast Cancer'),
        CNAQuery('ERBB2', 'Amplification', 'Colorectal Cancer'),
    ]

    annotations = pull_cna_info(queries)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Loss-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[2]
    assert len(annotation) == 14
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_2'


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_fake_cna():
    queries = [
        CNAQuery('test1', 'Amplification', 'Breast Cancer'),
    ]

    annotations = pull_cna_info(queries)
    fake_gene_one_query_suite(annotations)
