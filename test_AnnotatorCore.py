#!/usr/bin/python
import pytest

from AnnotatorCore import *


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

    # Test extreme cases
    with pytest.raises(Exception):
        assert resolve_query_type(None, [])
    assert resolve_query_type(None, [ALTERATION_HEADER]) == QueryType.HGVSP_SHORT

    # Raise exception when the file does not have asked header
    with pytest.raises(Exception):
        assert resolve_query_type(QueryType.HGVSG, [HGVSP_SHORT_HEADER])
    with pytest.raises(Exception):
        assert resolve_query_type(QueryType.GENOMIC_CHANGE, [GC_CHROMOSOME_HEADER, GC_START_POSITION_HEADER])
