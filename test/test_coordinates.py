#!/usr/bin/env python3
import sys
import os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')
from gencoor.exceptions import ChromosomeNotStrError, PositionNotIntegerError, NameNotStrError,\
                               StrandNotStrError, CoordinateFlipError
from gencoor.coordinates import GenCoor
import pytest

# @pytest.fixture(scope="module")


def test_GenCoor_ChromosomeNotStrError():
    with pytest.raises(ChromosomeNotStrError) as e_info:
        g1 = GenCoor(chrom=2, start=1, end=100, name="g1", strand=".")


def test_GenCoor_PositionNotIntegerError():
    with pytest.raises(PositionNotIntegerError) as e_info:
        g1 = GenCoor(chrom="chr1", start="1", end=100, name="g1", strand=".")
    with pytest.raises(PositionNotIntegerError) as e_info:
        g1 = GenCoor(chrom="chr1", start=1, end="100", name="g1", strand=".")


def test_GenCoor_NameNotStrError():
    with pytest.raises(NameNotStrError) as e_info:
        g1 = GenCoor(chrom="chr1", start=1, end=100, name=1000, strand=".")


def test_GenCoor_StrandNotStrError():
    with pytest.raises(StrandNotStrError) as e_info:
        g1 = GenCoor(chrom="chr1", start=1, end=100, name="g1", strand=100)


def test_GenCoor_CoordinateFlipError():
    with pytest.raises(CoordinateFlipError) as e_info:
        g1 = GenCoor(chrom="chr1", start=100, end=1, name="g1", strand=".")


def test_GenCoor_capital_name():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="capital", strand=".")
    g1.capital_name()
    assert g1.name == "CAPITAL"


def test_GenCoor_overlap1():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr2", start=1, end=100, name="test", strand=".")
    assert g1.overlap(a_gencoor=g2) is False


def test_GenCoor_overlap2():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=101, end=110, name="test", strand=".")
    assert g1.overlap(a_gencoor=g2) is False


def test_GenCoor_overlap3():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=100, end=110, name="test", strand=".")
    assert g1.overlap(a_gencoor=g2) is False


def test_GenCoor_overlap4():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=99, end=110, name="test", strand=".")
    assert g1.overlap(a_gencoor=g2) is True


def test_GenCoor_overlap5():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand="+")
    g2 = GenCoor(chrom="chr1", start=99, end=110, name="test", strand="-")
    assert g1.overlap(a_gencoor=g2, strandness=True) is False


def test_GenCoor_overlap6():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand="-")
    g2 = GenCoor(chrom="chr1", start=99, end=110, name="test", strand="-")
    assert g1.overlap(a_gencoor=g2, strandness=True) is True


def test_GenCoor_overlap7():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=50, end=50, name="test", strand=".")
    assert g1.overlap(a_gencoor=g2) is True


def test_GenCoor_overlap8():
    g1 = GenCoor(chrom="chr1", start=50, end=50, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=50, end=50, name="test", strand=".")
    assert g1.overlap(a_gencoor=g2) is True

