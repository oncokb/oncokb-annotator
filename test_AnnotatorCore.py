#!/usr/bin/env python
import pytest
import csv

import AnnotatorCore

from AnnotatorCore import getgenesfromfusion
from AnnotatorCore import conversion
from AnnotatorCore import replace_all
from AnnotatorCore import resolve_query_type
from AnnotatorCore import resolve_query_type_for_row
from AnnotatorCore import get_row_mutation_status
from AnnotatorCore import get_highest_tx_level
from AnnotatorCore import get_cna
from AnnotatorCore import get_oncokb_annotation_column_headers
from AnnotatorCore import processalterationevents
from AnnotatorCore import StructuralVariantQuery
from AnnotatorCore import QueryType
from AnnotatorCore import ALTERATION_HEADER
from AnnotatorCore import HGVSP_HEADER
from AnnotatorCore import HGVSP_SHORT_HEADER
from AnnotatorCore import HGVSG_HEADER
from AnnotatorCore import HGVSC_HEADER
from AnnotatorCore import GC_REF_ALLELE_HEADER
from AnnotatorCore import GC_CHROMOSOME_HEADER
from AnnotatorCore import GC_START_POSITION_HEADER
from AnnotatorCore import GC_END_POSITION_HEADER
from AnnotatorCore import GC_VAR_ALLELE_1_HEADER
from AnnotatorCore import GC_VAR_ALLELE_2_HEADER
from AnnotatorCore import TX_TYPE_SENSITIVE
from AnnotatorCore import TX_TYPE_RESISTANCE
from AnnotatorCore import CNA_AMPLIFICATION_TXT
from AnnotatorCore import CNA_DELETION_TXT
from AnnotatorCore import CNA_GAIN_TXT
from AnnotatorCore import CNA_LOSS_TXT
from AnnotatorCore import MutationStatus


def test_getgenesfromfusion():
    AB_EXAMPLE = ('A', 'B')
    assert getgenesfromfusion('A-B') == AB_EXAMPLE
    assert getgenesfromfusion('A-B ') == AB_EXAMPLE
    assert getgenesfromfusion('a-b') == ('a', 'b')
    assert getgenesfromfusion('A') == ('A', 'A')
    assert getgenesfromfusion('A1-1B') == ('A1', '1B')

    # Test fusion case insensitive
    assert getgenesfromfusion('A-B fusion') == AB_EXAMPLE
    assert getgenesfromfusion('A-B Fusion') == AB_EXAMPLE

    # Test unnecessary characters will be trimmed off after fusion
    assert getgenesfromfusion('A-B fusion archer') == AB_EXAMPLE
    assert getgenesfromfusion('A-B fusion Archer') == AB_EXAMPLE
    assert getgenesfromfusion('A-B fusion -Archer') == AB_EXAMPLE
    assert getgenesfromfusion('A-B fusion -archer') == AB_EXAMPLE
    assert getgenesfromfusion('A-B fusion - archer') == AB_EXAMPLE
    assert getgenesfromfusion('A-B fusion - archer ') == AB_EXAMPLE

    assert getgenesfromfusion('A-B fusion test') == AB_EXAMPLE
    assert getgenesfromfusion('fusion A-B fusion') == AB_EXAMPLE

    # Test intragenic
    assert getgenesfromfusion('MLL2-intragenic') == ('MLL2', 'MLL2')


def test_conversion():
    # Test conversion case for case insensitivity
    assert conversion('tyr100') == 'Y100'
    assert conversion('tYr100') == 'Y100'
    assert conversion('Tyr100') == 'Y100'
    assert conversion('tyR100') == 'Y100'
    assert conversion('TyR100') == 'Y100'
    assert conversion('TYR100') == 'Y100'
    assert conversion('tYR100') == 'Y100'
    assert conversion('sEr100') == 'S100'

    # Test conversion only targets dict() keys
    assert conversion('hot100') == 'hot100'

    # Test conversion is not affected by empty string and whitespaces
    assert conversion('') == ''
    assert conversion(' sEr100') == ' S100'

    # Test conversion when the string contains three letter but not supposed to be converted
    assert conversion('Promoter') == 'Promoter'


def test_replace_all():
    # Test replace_all for case insensitivity
    assert replace_all('tyr') == 'Y'
    assert replace_all('tYr') == 'Y'
    assert replace_all('Tyr') == 'Y'
    assert replace_all('tyR') == 'Y'
    assert replace_all('TyR') == 'Y'
    assert replace_all('TYR') == 'Y'
    assert replace_all('tYR') == 'Y'
    assert replace_all('sEr') == 'S'

    # Test replace_all only targets the dict() keys
    assert replace_all('bubblegum juice cup dairy hot pot Tyr melon') == 'bubblegum juice cup dairy hot pot Y melon'
    assert replace_all('Ly Lys Pr Pro Gln Glad Ph PH Phe') == 'Ly K Pr P Q Glad Ph PH F'
    assert replace_all(
        'nOt can fat Tan Rat cat dog man Men FAn rot taR car fAr map TAP Zip poP') == 'nOt can fat Tan Rat cat dog man Men FAn rot taR car fAr map TAP Zip poP'

    # Test replace_all is not affected by numbers
    assert replace_all('Tyr600E Cys56734342342454562456') == 'Y600E C56734342342454562456'
    assert replace_all(
        '60 045 434 345 4 26 567 254 245 34 67567 8 56 8 364 56 6 345 7567 3455 6 8 99 89 7 3') == '60 045 434 345 4 26 567 254 245 34 67567 8 56 8 364 56 6 345 7567 3455 6 8 99 89 7 3'

    # Test replace_all is not affected by empty string and whitespaces
    assert replace_all('') == ''
    assert replace_all(' ') == ' '
    assert replace_all('Tyr Asn As n Ile Il e') == 'Y N As n I Il e'


def test_resolve_query_type():
    assert resolve_query_type(None, [HGVSG_HEADER]) == QueryType.HGVSG
    assert resolve_query_type(None, [HGVSP_HEADER]) == QueryType.HGVSP
    assert resolve_query_type(None, [HGVSP_SHORT_HEADER]) == QueryType.HGVSP_SHORT
    assert resolve_query_type(None, [HGVSG_HEADER, HGVSP_HEADER, HGVSP_SHORT_HEADER]) == QueryType.HGVSP_SHORT
    assert resolve_query_type(None, [GC_CHROMOSOME_HEADER, GC_START_POSITION_HEADER, GC_END_POSITION_HEADER,
                                     GC_REF_ALLELE_HEADER, GC_VAR_ALLELE_1_HEADER,
                                     GC_VAR_ALLELE_2_HEADER]) == QueryType.GENOMIC_CHANGE

    assert resolve_query_type(QueryType.HGVSG, [HGVSG_HEADER, HGVSP_HEADER, HGVSP_SHORT_HEADER]) == QueryType.HGVSG

    assert resolve_query_type(None, [HGVSC_HEADER], MutationStatus.GERMLINE) == QueryType.HGVSC
    assert resolve_query_type(QueryType.HGVSC, [HGVSC_HEADER], MutationStatus.GERMLINE) == QueryType.HGVSC
    assert resolve_query_type(None, [ALTERATION_HEADER], MutationStatus.GERMLINE) == QueryType.HGVSC
    assert resolve_query_type(None, [HGVSG_HEADER], MutationStatus.GERMLINE) == QueryType.HGVSG
    assert resolve_query_type(QueryType.HGVSG, [HGVSG_HEADER], MutationStatus.GERMLINE) == QueryType.HGVSG
    assert resolve_query_type(None, [GC_CHROMOSOME_HEADER, GC_START_POSITION_HEADER, GC_END_POSITION_HEADER,
                                     GC_REF_ALLELE_HEADER, GC_VAR_ALLELE_1_HEADER,
                                     GC_VAR_ALLELE_2_HEADER], MutationStatus.GERMLINE) == QueryType.GENOMIC_CHANGE

    # Test extreme cases
    with pytest.raises(Exception):
        assert resolve_query_type(None, [])
    assert resolve_query_type(None, [ALTERATION_HEADER]) == QueryType.HGVSP_SHORT

    # Raise exception when the file does not have asked header
    with pytest.raises(Exception):
        assert resolve_query_type(QueryType.HGVSG, [HGVSP_SHORT_HEADER])
    with pytest.raises(Exception):
        assert resolve_query_type(QueryType.GENOMIC_CHANGE, [GC_CHROMOSOME_HEADER, GC_START_POSITION_HEADER])

    with pytest.raises(Exception):
        assert resolve_query_type(QueryType.HGVSP_SHORT, [HGVSP_SHORT_HEADER], MutationStatus.GERMLINE)

    with pytest.raises(Exception):
        assert resolve_query_type(QueryType.HGVSC, [HGVSC_HEADER], MutationStatus.SOMATIC)


def test_get_row_mutation_status():
    assert get_row_mutation_status(['Germline'], 0) == MutationStatus.GERMLINE
    assert get_row_mutation_status(['germline'], 0) == MutationStatus.GERMLINE
    assert get_row_mutation_status(['Somatic'], 0) == MutationStatus.SOMATIC
    assert get_row_mutation_status(['LOH'], 0) == MutationStatus.SOMATIC
    assert get_row_mutation_status([''], 0) == MutationStatus.SOMATIC


def test_resolve_query_type_for_row():
    headers = {
        HGVSP_SHORT_HEADER: 0,
        HGVSP_HEADER: 1,
        HGVSG_HEADER: 2,
        HGVSC_HEADER: 3,
        ALTERATION_HEADER: 4,
    }
    germline_row = ['p.V600E', '', '', 'c.2369C>T', '']
    assert resolve_query_type_for_row(germline_row, headers, None, MutationStatus.GERMLINE) == QueryType.HGVSC

    germline_invalid_row = ['p.V600E', '', '', '', '']
    assert resolve_query_type_for_row(germline_invalid_row, headers, None, MutationStatus.GERMLINE) is None

    somatic_hgvsg_row = ['', '', '7:g.140453136A>T', '', '']
    assert resolve_query_type_for_row(somatic_hgvsg_row, headers, None, MutationStatus.SOMATIC) == QueryType.HGVSG

    somatic_alteration_row = ['', '', '', '', 'MSI-H']
    assert resolve_query_type_for_row(somatic_alteration_row, headers, None, MutationStatus.SOMATIC) == QueryType.HGVSP_SHORT


def test_combined_headers():
    headers = get_oncokb_annotation_column_headers(False, True, MutationStatus.SOMATIC, combined_mode=True)
    assert 'ONCOGENIC' in headers
    assert 'PATHOGENIC' in headers
    assert 'PENETRANCE' in headers
    assert 'GENOMIC_INDICATOR' in headers
    assert 'ONCOKB_HUGO_SYMBOL' in headers


def test_structural_variant_query_preserves_intragenic_type():
    query = StructuralVariantQuery('BCOR', 'BCOR', 'INVERSION', 'PAAD')

    assert query.geneA.hugoSymbol == 'BCOR'
    assert query.geneB.hugoSymbol == 'BCOR'
    assert query.functionalFusion is False
    assert query.structuralVariantType == 'INVERSION'


def test_processalterationevents_mixed(tmp_path, monkeypatch):
    input_file = tmp_path / "mixed.maf"
    output_file = tmp_path / "mixed.out.maf"
    input_file.write_text(
        "\t".join([
            "Hugo_Symbol", "Tumor_Sample_Barcode", "Mutation_Status", "HGVSp_Short", "HGVSc", "HGVSg"
        ]) + "\n"
        + "\t".join(["BRAF", "S1", "Somatic", "p.V600E", "", ""]) + "\n"
        + "\t".join(["BRCA2", "S2", "Germline", "", "c.1166del", ""]) + "\n"
        + "\t".join(["EGFR", "S3", "Germline", "p.T790M", "", ""]) + "\n"
    )

    calls = []

    def fake_pull_protein_change_info(queries, include_descriptions, annotate_hotspot, mutation_status, combined_mode):
        calls.append(("protein", [query.alteration for query in queries], mutation_status, combined_mode))
        return [['True'] for _ in queries]

    def fake_pull_hgvsc_info(queries, include_descriptions, annotate_hotspot, mutation_status, combined_mode):
        calls.append(("hgvsc", [query.hgvsc for query in queries], mutation_status, combined_mode))
        return [['True'] for _ in queries]

    monkeypatch.setattr(AnnotatorCore, 'pull_protein_change_info', fake_pull_protein_change_info)
    monkeypatch.setattr(AnnotatorCore, 'pull_hgvsc_info', fake_pull_hgvsc_info)

    processalterationevents(
        str(input_file), str(output_file), '', '', {}, False, None, None, False, None
    )

    assert ("protein", ["V600E"], MutationStatus.SOMATIC, True) in calls
    assert ("hgvsc", ["BRCA2:c.1166del"], MutationStatus.GERMLINE, True) in calls

    with output_file.open() as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        rows = list(reader)

    annotated_index = header.index('ANNOTATED')
    oncogenic_index = header.index('ONCOGENIC')
    pathogenic_index = header.index('PATHOGENIC')
    assert rows[0][annotated_index] == 'True'
    assert rows[1][annotated_index] == 'True'
    assert rows[1][oncogenic_index] == ''
    assert rows[2][annotated_index] == 'False'
    assert rows[2][pathogenic_index] == ''


def test_get_highest_tx_level():
    oncokb_data = {}
    assert get_highest_tx_level(oncokb_data) == ''
    assert get_highest_tx_level(oncokb_data, 'random') == ''
    assert get_highest_tx_level(oncokb_data, TX_TYPE_SENSITIVE) == ''
    assert get_highest_tx_level(oncokb_data, TX_TYPE_RESISTANCE) == ''

    oncokb_data = {'LEVEL_1': ['test'], 'LEVEL_R1': ['test'], 'LEVEL_R2': ['test']}
    assert get_highest_tx_level(oncokb_data) == 'LEVEL_R1'
    assert get_highest_tx_level(oncokb_data, 'random') == 'LEVEL_R1'
    assert get_highest_tx_level(oncokb_data, TX_TYPE_SENSITIVE) == 'LEVEL_1'
    assert get_highest_tx_level(oncokb_data, TX_TYPE_RESISTANCE) == 'LEVEL_R1'

    oncokb_data = {'LEVEL_1': ['test'], 'LEVEL_R2': ['test']}
    assert get_highest_tx_level(oncokb_data) == 'LEVEL_1'
    assert get_highest_tx_level(oncokb_data, 'random') == 'LEVEL_1'
    assert get_highest_tx_level(oncokb_data, TX_TYPE_SENSITIVE) == 'LEVEL_1'
    assert get_highest_tx_level(oncokb_data, TX_TYPE_RESISTANCE) == 'LEVEL_R2'


def test_cna():
    assert get_cna(None) is None
    assert get_cna('') is None
    assert get_cna('test') is None
    assert get_cna('Amplification') == CNA_AMPLIFICATION_TXT
    assert get_cna('Gain') is None
    assert get_cna('Deletion') == CNA_DELETION_TXT
    assert get_cna('Loss') is None
    assert get_cna('2') == CNA_AMPLIFICATION_TXT
    assert get_cna('1') is None
    assert get_cna('-2') == CNA_DELETION_TXT
    assert get_cna('-1.5') == CNA_DELETION_TXT
    assert get_cna('-1') is None
    assert get_cna('0') is None

    assert get_cna(None, False) is None
    assert get_cna('', False) is None
    assert get_cna('test', False) is None
    assert get_cna('Amplification', False) == CNA_AMPLIFICATION_TXT
    assert get_cna('Gain', False) is None
    assert get_cna('Deletion', False) == CNA_DELETION_TXT
    assert get_cna('Loss', False) is None

    assert get_cna(None, True) is None
    assert get_cna('', True) is None
    assert get_cna('test', True) is None
    assert get_cna('Amplification', True) == CNA_AMPLIFICATION_TXT
    assert get_cna('Gain', True) == CNA_GAIN_TXT
    assert get_cna('Deletion', True) == CNA_DELETION_TXT
    assert get_cna('Loss', True) == CNA_LOSS_TXT
    assert get_cna('2', True) == CNA_AMPLIFICATION_TXT
    assert get_cna('1', True) == CNA_GAIN_TXT
    assert get_cna('-2', True) == CNA_DELETION_TXT
    assert get_cna('-1.5', True) == CNA_DELETION_TXT
    assert get_cna('-1', True) == CNA_LOSS_TXT
    assert get_cna('0', True) is None
