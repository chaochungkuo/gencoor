#!/usr/bin/env python3
import sys
import os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')
from gencoor.exceptions import ChromosomeNotStrError, PositionNotIntegerError, NameNotStrError,\
                               StrandNotStrError, CoordinateFlipError
from gencoor.coordinates import GenCoor, GenCoorSet
# from gencoor.coordinates_cy import GenCoor, GenCoorSet
import pytest

# @pytest.fixture(scope="module")

#################################################################
#### Tests on GenCoor

# def test_GenCoor_ChromosomeNotStrError():
#     with pytest.raises(ChromosomeNotStrError) as e_info:
#         g1 = GenCoor(chrom=2, start=1, end=100, name="g1", strand=".")
#
#
# def test_GenCoor_PositionNotIntegerError():
#     with pytest.raises(PositionNotIntegerError) as e_info:
#         g1 = GenCoor(chrom="chr1", start="1", end=100, name="g1", strand=".")
#     with pytest.raises(PositionNotIntegerError) as e_info:
#         g1 = GenCoor(chrom="chr1", start=1, end="100", name="g1", strand=".")
#
#
# def test_GenCoor_NameNotStrError():
#     with pytest.raises(NameNotStrError) as e_info:
#         g1 = GenCoor(chrom="chr1", start=1, end=100, name=1000, strand=".")
#
#
# def test_GenCoor_StrandNotStrError():
#     with pytest.raises(StrandNotStrError) as e_info:
#         g1 = GenCoor(chrom="chr1", start=1, end=100, name="g1", strand=100)
#
#
# def test_GenCoor_CoordinateFlipError():
#     with pytest.raises(CoordinateFlipError) as e_info:
#         g1 = GenCoor(chrom="chr1", start=100, end=1, name="g1", strand=".")


def test_GenCoor_capital_name():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="capital", strand=".")
    g1.capital_name()
    assert g1.name == "CAPITAL"


def test_GenCoor_overlap1():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr2", start=1, end=100, name="test", strand=".")
    assert g1.overlap(region=g2) is False

def test_GenCoor_overlap_pairing1():
    res = GenCoor.overlap_pairing("chr1", 1, 100, ".", "chr2", 1, 100, ".")
    assert res is False

def test_GenCoor_overlap():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=101, end=110, name="test", strand=".")
    assert g1.overlap(region=g2) is False

def test_GenCoor_overlap_pairing2():
    res = GenCoor.overlap_pairing("chr1", 1, 100, ".", "chr1", 101, 110, ".")
    assert res is False

def test_GenCoor_overlap3():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=100, end=110, name="test", strand=".")
    assert g1.overlap(region=g2) is False

def test_GenCoor_overlap_pairing3():
    res = GenCoor.overlap_pairing("chr1", 1, 100, ".", "chr1", 100, 110, ".")
    assert res is False

def test_GenCoor_overlap4():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=99, end=110, name="test", strand=".")
    assert g1.overlap(region=g2) is True

def test_GenCoor_overlap_pairing4():
    res = GenCoor.overlap_pairing("chr1", 1, 100, ".", "chr1", 99, 110, ".")
    assert res is True

def test_GenCoor_overlap5():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand="+")
    g2 = GenCoor(chrom="chr1", start=99, end=110, name="test", strand="-")
    assert g1.overlap(region=g2, strand_specific=True) is False


def test_GenCoor_overlap6():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand="-")
    g2 = GenCoor(chrom="chr1", start=99, end=110, name="test", strand="-")
    assert g1.overlap(region=g2, strand_specific=True) is True


def test_GenCoor_overlap7():
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=50, end=50, name="test", strand=".")
    assert g1.overlap(region=g2) is True


def test_GenCoor_overlap8():
    g1 = GenCoor(chrom="chr1", start=50, end=50, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=50, end=50, name="test", strand=".")
    assert g1.overlap(region=g2) is True


def test_GenCoor_distance1():
    g1 = GenCoor(chrom="chr1", start=10, end=50, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=60, end=70, name="test", strand=".")
    assert g1.distance(region=g2) == 10


def test_GenCoor_distance2():
    g1 = GenCoor(chrom="chr1", start=10, end=50, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=60, end=70, name="test", strand=".")
    assert g2.distance(region=g1) == 10


def test_GenCoor_distance3():
    g1 = GenCoor(chrom="chr1", start=10, end=50, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=60, end=70, name="test", strand=".")
    assert g2.distance(region=g1, sign=True) == -10


def test_GenCoor_distance4():
    g1 = GenCoor(chrom="chr1", start=10, end=65, name="test", strand=".")
    g2 = GenCoor(chrom="chr1", start=60, end=70, name="test", strand=".")
    assert g2.distance(region=g1, sign=True) == 0


def test_GenCoor_distance5():
    g1 = GenCoor(chrom="chr1", start=10, end=65, name="test", strand=".")
    g2 = GenCoor(chrom="chr2", start=60, end=70, name="test", strand=".")
    assert g2.distance(region=g1, sign=True) is None

def test_GenCoor_extract_bed12():
    g1 = GenCoor(chrom="chr22", start=1000, end=5000, name="cloneA", strand="+",
                 data="1000/t5000/t0/t2/t567, 488/t0, 3512", score=960)
    g1_extracts = g1.extract_bed12()
    assert len(g1_extracts) == 2
    assert len(g1_extracts[0]) == 567
    assert len(g1_extracts[1]) == 488

#################################################################
#### Tests on GenCoorSet




def test_GenCoorSet_add():
    genset = GenCoorSet(name="Test_set")
    genset.add(GenCoor(chrom="chr1", start=10, end=20, name="test", strand="."))
    genset.add(GenCoor(chrom="chr1", start=15, end=50, name="test", strand="."))
    genset.add(GenCoor(chrom="chr2", start=100, end=200, name="test", strand="."))
    assert len(genset.list) == 3


def test_GenCoorSet_len():
    genset = GenCoorSet(name="Test_set")
    genset.add(GenCoor(chrom="chr1", start=10, end=20, name="test", strand="."))
    genset.add(GenCoor(chrom="chr1", start=15, end=50, name="test", strand="."))
    genset.add(GenCoor(chrom="chr2", start=100, end=200, name="test", strand="."))
    assert len(genset) == 3

def test_GenCoorSet_load():
    genset = GenCoorSet(name="Test_set")
    genset.load(filename="/Users/jovesus/rgtdata/hg38/genes_Gencode_hg38.bed", filetype="BED")
    genset.save(filename="/Users/jovesus/projects/gencoor/test/test_res_bed.bed", filetype="BED")
    ### test BED12
    # genset.load(filename="/Users/jovesus/projects/gencoor/test/test_bed12.bed", filetype="BED12")
    # genset.save(filename="/Users/jovesus/projects/gencoor/test/test_res_bed12.bed", filetype="BED12")
def test_extend():
    genset = GenCoorSet(name="Test_set")
    genset.add(GenCoor(chrom="chr1", start=10, end=20, name="test", strand="+"))
    genset.add(GenCoor(chrom="chr1", start=15, end=50, name="test", strand="-"))
    genset.add(GenCoor(chrom="chr2", start=100, end=200, name="test", strand="."))
    ngcs = genset.extend(mode="left", length=5)
    assert ngcs.list[0].start == 5
    assert ngcs.list[0].end == 20
    assert ngcs.list[1].start == 10
    assert ngcs.list[1].end == 50
    ngcs = genset.extend(mode="right", length=5)
    assert ngcs.list[0].start == 10
    assert ngcs.list[0].end == 25
    assert ngcs.list[1].start == 15
    assert ngcs.list[1].end == 55
    ngcs = genset.extend(mode="5end", length=5)
    assert ngcs.list[0].start == 5
    assert ngcs.list[0].end == 20
    assert ngcs.list[1].start == 15
    assert ngcs.list[1].end == 55
    ngcs = genset.extend(mode="3end", length=5)
    assert ngcs.list[0].start == 10
    assert ngcs.list[0].end == 25
    assert ngcs.list[1].start == 10
    assert ngcs.list[1].end == 50
    ngcs = genset.extend(mode="both", length=5)
    assert ngcs.list[0].start == 5
    assert ngcs.list[0].end == 25
    assert ngcs.list[1].start == 10
    assert ngcs.list[1].end == 55
def test_split_by_strands():
    genset = GenCoorSet(name="Test_set")
    genset.load(filename="/Users/jovesus/rgtdata/hg38/genes_Gencode_hg38.bed", filetype="BED")
    genset_FWD, genset_REV = genset.split_by_strands()
    assert set([g.strand for g in genset_FWD]) == set(["+"])
    assert set([g.strand for g in genset_REV]) == set(["-"])
def test_merge():
    genset = GenCoorSet(name="Test_set")
    genset.add(GenCoor(chrom="chr1", start=10, end=20, name="test", strand="+"))
    genset.add(GenCoor(chrom="chr1", start=15, end=50, name="test", strand="-"))
    genset.add(GenCoor(chrom="chr2", start=100, end=200, name="test", strand="."))
    res = genset.merge(w_return=True)
    assert len(res) == 2
    assert res[0].start == 10
    assert res[0].end == 50
    assert res[0].strand == "."
    res = genset.merge(w_return=True, strand_specific=True)
    print(res[0])
    print(res[1])
    print(res[2])
    assert len(res) == 3
    assert res[0].start == 10
    assert res[0].end == 20
    assert res[0].strand == "+"

