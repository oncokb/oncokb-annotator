#!/usr/bin/env python3

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
    assert conversion('tyr') == 'Y'
    assert conversion('tYr') == 'Y'
    assert conversion('Tyr') == 'Y'
    assert conversion('tyR') == 'Y'
    assert conversion('TyR') == 'Y'
    assert conversion('TYR') == 'Y'
    assert conversion('tYR') == 'Y'
    assert conversion('sEr') == 'S'

    # Test conversion only targets dict() keys
    assert conversion('hot cup mop lap tap zip Thr cheetos hand') == 'hot cup mop lap tap zip T cheetos hand'
    assert conversion('Ly s Lys Th r Thr X X As p Asp') == 'Ly s K Th r T X X As p D'

    # Test conversion is not affected by numbers
    assert conversion('8453 388 830 32 -2 -4 -23 2 -23 6 26 784 1 2 3 4 5 6 7 8 90 0') == '8453 388 830 32 -2 -4 -23 2 -23 6 26 784 1 2 3 4 5 6 7 8 90 0'
    assert conversion('Val500Thr Asn 788 88 His 5H8is') == 'V500T N 788 88 H 5H8is'

    # Test conversion is not affected by empty string and whitespaces
    assert conversion('') == ''
    assert conversion(' ') == ' '
    assert conversion('Tyr Asn As n Ile Il e') == 'Y N As n I Il e'

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
