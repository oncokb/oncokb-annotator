#!/usr/bin/python
import pytest
import os
import logging

from AnnotatorCore import pull_hgvsg_info, DESCRIPTION_HEADERS, ONCOKB_ANNOTATION_HEADERS_GC
from AnnotatorCore import pull_genomic_change_info
from AnnotatorCore import pull_protein_change_info
from AnnotatorCore import pull_structural_variant_info
from AnnotatorCore import pull_cna_info
from AnnotatorCore import setoncokbapitoken
from AnnotatorCore import ProteinChangeQuery
from AnnotatorCore import GenomicChangeQuery
from AnnotatorCore import StructuralVariantQuery
from AnnotatorCore import CNAQuery
from AnnotatorCore import HGVSgQuery
from AnnotatorCore import ReferenceGenome

ONCOKB_API_TOKEN = os.environ["ONCOKB_API_TOKEN"]
setoncokbapitoken(ONCOKB_API_TOKEN)

log = logging.getLogger('test_Annotation')
log.info('test-----------', os.environ["ONCOKB_API_TOKEN"], '------')

VARIANT_EXISTS_INDEX = 2
MUTATION_EFFECT_INDEX = VARIANT_EXISTS_INDEX + 1
ONCOGENIC_INDEX = MUTATION_EFFECT_INDEX + 2
LEVEL_1_INDEX = ONCOGENIC_INDEX + 1
LEVEL_2_INDEX = LEVEL_1_INDEX + 1
LEVEL_3A_INDEX = LEVEL_1_INDEX + 2
HIGHEST_LEVEL_INDEX = LEVEL_1_INDEX + 7
HIGHEST_DX_LEVEL_INDEX = HIGHEST_LEVEL_INDEX + 7
HIGHEST_PX_LEVEL_INDEX = HIGHEST_DX_LEVEL_INDEX + 5
UNKNOWN = 'Unknown'
NUMBER_OF_ANNOTATION_COLUMNS = 27
NUMBER_OF_DESCRIPTION_COLUMNS = len(DESCRIPTION_HEADERS)
NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS = len(ONCOKB_ANNOTATION_HEADERS_GC)
NUMBER_OF_ANNOTATION_COLUMNS_WITH_DESCRIPTIONS = NUMBER_OF_ANNOTATION_COLUMNS + NUMBER_OF_DESCRIPTION_COLUMNS
NUMBER_OF_GC_ANNOTATION_COLUMNS = NUMBER_OF_ANNOTATION_COLUMNS + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS
NUMBER_OF_GC_ANNOTATION_COLUMNS_WITH_DESCRIPTIONS = NUMBER_OF_GC_ANNOTATION_COLUMNS + NUMBER_OF_DESCRIPTION_COLUMNS


def fake_gene_one_query_suite(annotations, include_descriptions):
    assert len(annotations) == 1

    annotation = annotations[0]
    assert len(
        annotation) == NUMBER_OF_ANNOTATION_COLUMNS if include_descriptions is False else NUMBER_OF_ANNOTATION_COLUMNS_WITH_DESCRIPTIONS
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == UNKNOWN
    assert annotation[HIGHEST_LEVEL_INDEX] == ''


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_protein_change():
    queries = [
        ProteinChangeQuery('BRAF', 'V600E', 'Colorectal Cancer'),
        ProteinChangeQuery('ABL1', 'BCR-ABL1 Fusion', 'Acute Leukemias of Ambiguous Lineage'),
    ]

    annotations = pull_protein_change_info(queries, False, False)
    assert len(annotations) == 2

    annotation = annotations[0]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
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

    annotations = pull_genomic_change_info(queries, False, False)
    assert len(annotations) == 2

    annotation37 = annotations[0]
    annotation38 = annotations[1]
    assert annotation37 == annotation38

    queries = [
        ProteinChangeQuery('MYD88', 'M232T', 'Ovarian Cancer', ReferenceGenome.GRCH37),
        ProteinChangeQuery('MYD88', 'M219T', 'Ovarian Cancer', ReferenceGenome.GRCH38)
    ]

    annotations = pull_protein_change_info(queries, False, False)
    assert len(annotations) == 2

    annotation37 = annotations[0]
    annotation38 = annotations[1]
    assert annotation37 == annotation38


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_fake_gene_protein_change():
    queries = [
        ProteinChangeQuery('test1', 'V600E', 'Ovarian Cancer')
    ]

    annotations = pull_protein_change_info(queries, False, False)
    fake_gene_one_query_suite(annotations, False)

    annotations = pull_protein_change_info(queries, False, False)
    fake_gene_one_query_suite(annotations, True)


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_atypical_alts():
    queries = [
        ProteinChangeQuery('Other Biomarkers', 'MSI-H', 'Colorectal Cancer'),
        ProteinChangeQuery('Other Biomarkers', 'MSI-H', 'Leukemia'),
        ProteinChangeQuery('TERT', 'Promoter Mutation', 'Bladder Cancer'),
        ProteinChangeQuery('TERT', 'Promoter Mutation', 'Bladder Cancer', None, '5\'Flank')
    ]

    annotations = pull_protein_change_info(queries, False, False)
    assert len(annotations) == 4

    annotation = annotations[0]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == UNKNOWN
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

    annotation = annotations[2]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == 'Likely Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Likely Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == ''

    annotation_dup = annotations[3]
    assert len(annotation_dup) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation == annotation_dup


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_hgvsg():
    queries = [
        # KRAF G12C
        HGVSgQuery('12:g.25398285C>A', 'LUAD'),
        # KRAF G12C
        HGVSgQuery('12:g.25398285_25398286delinsAG', 'LUAD'),
        # TERT Promoter
        HGVSgQuery('5:g.1295167_1295168delinsAATG', 'LUAD'),
    ]

    annotations = pull_hgvsg_info(queries, False, False)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == NUMBER_OF_GC_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == NUMBER_OF_GC_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'LEVEL_1'

    annotation = annotations[2]
    assert len(annotation) == NUMBER_OF_GC_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Likely Gain-of-function'
    assert annotation[ONCOGENIC_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Likely Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == ''


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_genomic_change():
    queries = [
        # KRAF G12C
        GenomicChangeQuery('12', '25398285', '25398285', 'C', 'A', 'LUAD'),
        # KRAF G12C
        GenomicChangeQuery('12', '25398285', '25398286', 'CA', 'AG', 'LUAD'),
        # TERT Promoter
        GenomicChangeQuery('5', '1295167', '1295168', 'TC', 'AATG', 'LUAD'),
    ]

    annotations = pull_genomic_change_info(queries, False, False)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == NUMBER_OF_GC_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == NUMBER_OF_GC_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'LEVEL_1'

    annotation = annotations[2]
    assert len(annotation) == NUMBER_OF_GC_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Likely Gain-of-function'
    assert annotation[ONCOGENIC_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == 'Likely Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX + NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS] == ''


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_check_structural_variants():
    queries = [
        StructuralVariantQuery('ALK', 'EML4', 'FUSION', 'NSCLC'),
        StructuralVariantQuery('ALK', 'EML4', 'FUSION', 'Melanoma'),
        StructuralVariantQuery('BCR', 'ABL1', 'FUSION', 'Acute Leukemias of Ambiguous Lineage'),
    ]

    annotations = pull_structural_variant_info(queries, False)
    assert len(annotations) == 3

    annotation = annotations[0]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_3B'

    annotation = annotations[2]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
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

    annotations = pull_structural_variant_info(queries, False)
    fake_gene_one_query_suite(annotations, False)

    annotations = pull_structural_variant_info(queries, False)
    fake_gene_one_query_suite(annotations, True)


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_cna():
    queries = [
        CNAQuery('BRCA2', 'DELETION', 'Ovarian Cancer'),
        CNAQuery('ERBB2', 'Amplification', 'Breast Cancer'),
        CNAQuery('ERBB2', 'Amplification', 'Colorectal Cancer'),
        CNAQuery('CDKN2A', 'Deletion', 'AML with BCR-ABL1'),
    ]

    annotations = pull_cna_info(queries, False)
    assert len(annotations) == 4

    annotation = annotations[0]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == 'Loss-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[1]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[2]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[MUTATION_EFFECT_INDEX] == 'Gain-of-function'
    assert annotation[ONCOGENIC_INDEX] == 'Oncogenic'
    assert annotation[HIGHEST_LEVEL_INDEX] == 'LEVEL_1'

    annotation = annotations[3]
    assert len(annotation) == NUMBER_OF_ANNOTATION_COLUMNS
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

    annotations = pull_cna_info(queries, False)
    fake_gene_one_query_suite(annotations, False)

    annotations = pull_cna_info(queries, True)
    fake_gene_one_query_suite(annotations, True)


def check_brca2_s1882_without_cancertype(annotation, genomic_query=False):
    assert len(annotation) == NUMBER_OF_GC_ANNOTATION_COLUMNS if genomic_query else NUMBER_OF_ANNOTATION_COLUMNS
    assert annotation[(
            NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS + MUTATION_EFFECT_INDEX) if genomic_query else MUTATION_EFFECT_INDEX] == 'Likely Loss-of-function'
    assert annotation[(
            NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS + ONCOGENIC_INDEX) if genomic_query else ONCOGENIC_INDEX] == 'Likely Oncogenic'
    assert annotation[(
            NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS + HIGHEST_LEVEL_INDEX) if genomic_query else HIGHEST_LEVEL_INDEX] == 'LEVEL_1'
    assert annotation[(
            NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS + LEVEL_1_INDEX) if genomic_query else LEVEL_1_INDEX] == 'Olaparib,Olaparib+Bevacizumab,Rucaparib,Olaparib+Abiraterone+Prednisone,Niraparib,Olaparib+Abiraterone+Prednisolone,Talazoparib+Enzalutamide,Niraparib+Abiraterone Acetate+Prednisone'
    assert annotation[(
            NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS + LEVEL_2_INDEX) if genomic_query else LEVEL_2_INDEX] == 'Olaparib,Rucaparib,Niraparib'
    assert annotation[(
            NUMBER_OF_ONCOKB_ANNOTATION_GC_COLUMNS + LEVEL_3A_INDEX) if genomic_query else LEVEL_3A_INDEX] == 'Olaparib,Talazoparib'


@pytest.mark.skipif(ONCOKB_API_TOKEN in (None, ''), reason="oncokb api token required")
def test_duplicated_treatments():
    # there should not be any duplicated treatment listed when cancer type is not specified

    # test protein change query
    queries = [
        ProteinChangeQuery('BRCA2', 'S1882*', ''),
    ]
    annotations = pull_protein_change_info(queries, False, False)
    assert len(annotations) == 1

    check_brca2_s1882_without_cancertype(annotations[0])

    # test genomic change query
    queries = [
        GenomicChangeQuery('13', '32914137', '32914137', 'C', 'A', ''),
    ]
    annotations = pull_genomic_change_info(queries, False, False)
    assert len(annotations) == 1

    check_brca2_s1882_without_cancertype(annotations[0], True)
