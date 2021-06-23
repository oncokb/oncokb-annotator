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
HIGHEST_DX_LEVEL_INDEX = 17
HIGHEST_PX_LEVEL_INDEX = 21
UNKNOWN = 'Unknown'

def fake_gene_one_query_suite(annotations):
    assert len(annotations) == 1

    annotation = annotations[0]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == ''
    assert annotation[HIGHEST_LEVEL_INDEX] == ''


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_protein_change():
    queries = [
        ProteinChangeQuery('BRAF', 'V600E', 'Colorectal Cancer'),
        ProteinChangeQuery('ABL1', 'BCR-ABL1 Fusion', 'Acute Leukemias of Ambiguous Lineage'),
    ]

    annotations = pull_protein_change_info(queries, False)
    assert len(annotations) == 2

    annotation = annotations[0]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''
    assert annotation[HIGHEST_DX_LEVEL_INDEX] == 'LEVEL_Dx1'
    assert annotation[HIGHEST_PX_LEVEL_INDEX] == 'LEVEL_Px1'


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_reference_genome():
    queries = [
        GenomicChangeQuery('7', '140453136', '140453136', 'A', 'T', 'LUAD', ReferenceGenome.GRCH37),
        GenomicChangeQuery('7', '140753336', '140753336', 'A', 'T', 'LUAD', ReferenceGenome.GRCH38)
    ]

    annotations = pull_protein_change_info(queries, False)
    assert len(annotations) == 2

    annotation37 = annotations[0]
    annotation38 = annotations[1]
    assert annotation37 == annotation38

    queries = [
        ProteinChangeQuery('MYD88', 'M232T', 'Ovarian Cancer', ReferenceGenome.GRCH37),
        ProteinChangeQuery('MYD88', 'M219T', 'Ovarian Cancer', ReferenceGenome.GRCH38)
    ]

    annotations = pull_protein_change_info(queries, False)
    assert len(annotations) == 2

    annotation37 = annotations[0]
    annotation38 = annotations[1]
    assert annotation37 == annotation38


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_fake_gene_protein_change():
    queries = [
        ProteinChangeQuery('test1', 'V600E', 'Ovarian Cancer')
    ]

    annotations = pull_protein_change_info(queries, False)
    fake_gene_one_query_suite(annotations)


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_atypical_alts():
    queries = [
        ProteinChangeQuery('Other Biomarkers', 'MSI-H', 'Colorectal Cancer'),
        ProteinChangeQuery('Other Biomarkers', 'MSI-H', 'Leukemia'),
        ProteinChangeQuery('TERT', 'Promoter Mutation', 'Bladder Cancer'),
        ProteinChangeQuery('TERT', 'Promoter Mutation', 'Bladder Cancer', None, '5\'Flank')
    ]

    annotations = pull_protein_change_info(queries, False)
    assert len(annotations) == 4

    annotation = annotations[0]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

    annotation = annotations[2]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

    annotation_dup = annotations[3]
    assert len(annotation_dup) == 22
    assert annotation == annotation_dup


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_hgvsg():
    queries = [
        HGVSgQuery('12:g.25398285C>A', 'LUAD'),
        HGVSgQuery('12:g.25398285_25398286delinsAG', 'LUAD'),
        HGVSgQuery('5:g.1295228G>A', 'LUAD'),
    ]

    annotations = pull_hgvsg_info(queries, False)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[2]
    assert len(annotation) == 22
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

    annotations = pull_genomic_change_info(queries, False)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[2]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_fusions():
    queries = [
        StructuralVariantQuery('ALK', 'EML4', 'FUSION', 'NSCLC'),
        StructuralVariantQuery('ALK', 'EML4', 'FUSION', 'Melanoma'),
        StructuralVariantQuery('BCR', 'ABL1', 'FUSION', 'Acute Leukemias of Ambiguous Lineage'),
    ]

    annotations = pull_structural_variant_info(queries)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_3B'

    annotation = annotations[2]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''
    assert annotation[HIGHEST_DX_LEVEL_INDEX] == 'LEVEL_Dx1'
    assert annotation[HIGHEST_PX_LEVEL_INDEX] == 'LEVEL_Px1'


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
        CNAQuery('CDKN2A', 'Deletion', 'AML with BCR-ABL1'),
    ]

    annotations = pull_cna_info(queries)
    assert len(annotations) == 4

    annotation = annotations[0]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Loss-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[2]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_2'

    annotation = annotations[3]
    assert len(annotation) == 22
    assert annotation[MUTATION_EFFECT_INDEX] == 'Loss-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''
    assert annotation[HIGHEST_DX_LEVEL_INDEX] == 'LEVEL_Dx2'
    assert annotation[HIGHEST_PX_LEVEL_INDEX] == ''


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_fake_cna():
    queries = [
        CNAQuery('test1', 'Amplification', 'Breast Cancer'),
    ]

    annotations = pull_cna_info(queries)
    fake_gene_one_query_suite(annotations)
