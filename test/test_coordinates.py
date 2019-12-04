#!/usr/bin/env python3
import sys
import os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')
from gencoor.exceptions import ChromosomeNotStrError, PositionNotIntegerError, NameNotStrError,\
                               StrandNotStrError, CoordinateFlipError
from gencoor.coordinates import GenCoor, GenCoorSet
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
    assert g2.distance(region=g1, sign=True) == None

def test_GenCoor_extract_bed12():
    g1 = GenCoor(chrom="chr22", start=1000, end=5000, name="cloneA", strand="+",
                 data="1000/t5000/t0/t2/t567, 488/t0, 3512", score=960)
    g1_extracts = g1.extract_bed12()
    assert len(g1_extracts) == 2
    assert len(g1_extracts[0]) == 567
    assert len(g1_extracts[1]) == 488


        #
        # '3end as center'
        # '5end as 5end'
        # '3end as 3end'

def test_GenCoor_relocate1():
    g1 = GenCoor(chrom="chr1", start=10, end=20, name="test", strand=".")
    r = g1.relocate(mode='center as center', inplace=False)
    assert r.start == 10
    assert r.end == 20
    r = g1.relocate(mode='5end as center', inplace=False)
    assert r.start == 5
    assert r.end == 15
    r = g1.relocate(mode='3end as center', inplace=False)
    assert r.start == 15
    assert r.end == 25
    r = g1.relocate(mode='5end as 5end', inplace=False)
    assert r.start == 10
    assert r.end == 20
    r = g1.relocate(mode='3end as 3end', inplace=False)
    assert r.start == 10
    assert r.end == 20
def test_GenCoor_relocate2():
    g1 = GenCoor(chrom="chr1", start=10, end=20, name="test", strand="+")
    r = g1.relocate(mode='center as center', inplace=False)
    assert r.start == 10
    assert r.end == 20
    r = g1.relocate(mode='5end as center', inplace=False)
    assert r.start == 5
    assert r.end == 15
    r = g1.relocate(mode='3end as center', inplace=False)
    assert r.start == 15
    assert r.end == 25
    r = g1.relocate(mode='5end as 5end', inplace=False)
    assert r.start == 10
    assert r.end == 20
    r = g1.relocate(mode='3end as 3end', inplace=False)
    assert r.start == 10
    assert r.end == 20
def test_GenCoor_relocate3():
    g1 = GenCoor(chrom="chr1", start=10, end=20, name="test", strand="-")
    r = g1.relocate(mode='center as center', inplace=False)
    assert r.start == 10
    assert r.end == 20
    r = g1.relocate(mode='5end as center', inplace=False)
    assert r.start == 15
    assert r.end == 25
    r = g1.relocate(mode='3end as center', inplace=False)
    assert r.start == 5
    assert r.end == 15
    r = g1.relocate(mode='5end as 5end', inplace=False)
    assert r.start == 10
    assert r.end == 20
    r = g1.relocate(mode='3end as 3end', inplace=False)
    assert r.start == 10
    assert r.end == 20
def test_GenCoor_relocate4():
    g1 = GenCoor(chrom="chr1", start=10, end=20, name="test", strand="+")
    r = g1.relocate(mode='center as center', width=2, inplace=False)
    assert r.start == 14
    assert r.end == 16
    r = g1.relocate(mode='5end as center', width=2, inplace=False)
    assert r.start == 9
    assert r.end == 11
    r = g1.relocate(mode='3end as center', width=2, inplace=False)
    assert r.start == 19
    assert r.end == 21
    r = g1.relocate(mode='5end as 5end', width=2, inplace=False)
    assert r.start == 10
    assert r.end == 12
    r = g1.relocate(mode='3end as 3end', width=2, inplace=False)
    assert r.start == 18
    assert r.end == 20
def test_GenCoor_relocate5():
    g1 = GenCoor(chrom="chr1", start=10, end=20, name="test", strand="-")
    r = g1.relocate(mode='center as center', width=2, inplace=False)
    assert r.start == 14
    assert r.end == 16
    r = g1.relocate(mode='5end as center', width=2, inplace=False)
    assert r.start == 19
    assert r.end == 21
    r = g1.relocate(mode='3end as center', width=2, inplace=False)
    assert r.start == 9
    assert r.end == 11
    r = g1.relocate(mode='5end as 5end', width=2, inplace=False)
    assert r.start == 18
    assert r.end == 20
    r = g1.relocate(mode='3end as 3end', width=2, inplace=False)
    assert r.start == 10
    assert r.end == 12

#################################################################
#### Tests on GenCoorSet
#################################################################

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
    bedfile = os.path.join(os.getenv("HOME"), "gencoor_data/hg38/genes_hg38.bed")
    bedfile2 = os.path.join(os.getenv("HOME"), "gencoor_data/hg38/genes_hg38_test.bed")
    genset.load(filename=bedfile, filetype="BED")
    genset.save(filename=bedfile2, filetype="BED")
    os.remove(bedfile2)
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
    bedfile = os.path.join(os.getenv("HOME"), "gencoor_data/hg38/genes_hg38.bed")
    genset.load(filename=bedfile, filetype="BED")
    res = genset.split_by_strands()
    assert set([g.strand for g in res["+"]]) == set(["+"])
    assert set([g.strand for g in res["-"]]) == set(["-"])
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
    assert len(res) == 3
    assert res[0].start == 10
    assert res[0].end == 20
    assert res[0].strand == "+"
def test_intersect_1():
    """
    Two empty sets
    A : none
    B : none
    R : none
    """
    genset1 = GenCoorSet(name="Test_set")
    genset2 = GenCoorSet(name="Test_set")
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0
def test_intersect_2():
    """
    One empty set
    A :   -----
    B : none
    R : none
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=10, end=20, name="test", strand="+"))
    genset2 = GenCoorSet(name="Test_set")
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0
def test_intersect_3():
    """
    A : none
    B :   -----
    R : none
    """
    genset1 = GenCoorSet(name="Test_set")
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=10, end=20, name="test", strand="+"))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0
def test_intersect_4():
    """
    No overlapping
    A : ------      ---------               -------
    B :        ----          ------  ------
    R : none
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=11, end=20, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=33, end=38, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=7, end=9, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=20, end=25, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=26, end=31, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0

def test_intersect_5():
    """
    End-to-end attach
    A : ------      ------
    B :       ------
    R : none
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=11, end=20, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=5, end=11, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0

def test_intersect_6():
    """
    No length attach
    A : .      .
    B :    .   .
    R : none
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=2, end=2, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=20, end=20, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=5, end=5, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=20, end=20, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 1
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 1
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 1

def test_intersect_7():
    """
    Perfect overlapping
    A : ------
    B : ------
    R : ------
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=10, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=500, end=550, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=600, end=650, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=700, end=750, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=725, end=800, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=1, end=10, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=500, end=550, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=600, end=650, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=700, end=750, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=725, end=800, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 6
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 5
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 5

def test_intersect_8():
    """
    One overlapping region
    A : ------
    B :     --------
    R1:     --       (overlap)
    R2: ------       (original)
    R3:              (comp_incl)
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=10, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=7, end=20, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 1
    assert res[0].start == 7
    assert res[0].end == 10
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 1
    assert res[0].start == 1
    assert res[0].end == 10
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0

def test_intersect_9():
    """
    Two simple overlapping regions
    A : -------      --------
    B :     -------------
    R1:     ---      ----     (overlap)
    R2: -------      -------- (original)
    R3:                       (comp_incl)
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=10, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=26, end=35, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=7, end=30, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 2
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 2
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0

def test_intersect_10():
    """
    Two separately overlapping regions
    A : -------      --------
    B :     -----        --------
    R1:     ---          ----     (overlap)
    R2: -------      --------     (original)
    R3:                           (comp_incl)
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=10, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=26, end=35, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=7, end=15, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=30, end=40, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 2
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 2
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0

def test_intersect_11():
    """
    Many various overlapping (mixed)
    A :   ------------------            --------   ---------
    B : ----   -------    ------            ----------
    R1:   --   -------    --                ----   ---       (overlap)
    R2:   ------------------            --------   --------- (original)
    R3:                                                      (comp_incl)
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=3, end=30, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=50, end=60, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=70, end=85, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=10, end=19, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=27, end=35, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=55, end=75, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 5
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 3
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0

def test_intersect_12():
    """
    Different chromosomes
    A : chr1  -------
    B : chr2  -------
    R : none
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=10, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr2", start=1, end=10, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 0
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0

def test_intersect_13():
    """
    Completely included overlapping
    A : ---------------------------
    B : ----    ------       -----------
    R1: ----    ------       ------      (overlap)
    R2: ---------------------------      (original)
    R3:                                  (comp_incl)
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=50, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=10, end=19, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=45, end=60, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 3
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 1
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 0

def test_intersect_14():
    """
    A : ----    ------       -----------
    B : ---------------------------
    R1: ----    ------       ------      (overlap)
    R2: ----    ------       ----------- (original)
    R3: ----    ------                   (comp_incl)
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=10, end=19, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=45, end=60, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=1, end=50, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 3
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 3
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 2

def test_intersect_15():
    """
    A : --------------         -------
            ------
    B :       -----          ----------------
    R1:       -----            -------      (overlap)
              ----
    R2: --------------         -------      (original)
            ------
    R3:                        -------      (comp_incl)
    """
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=50, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=20, end=40, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=70, end=80, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=25, end=45, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=65, end=95, name="test", strand="."))
    res = genset1.intersect(genset2, mode="overlap")
    assert len(res) == 3
    res = genset1.intersect(genset2, mode="original")
    assert len(res) == 3
    res = genset1.intersect(genset2, mode="complete_included")
    assert len(res) == 1

def test_standard_chromosome1():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=50, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr2", start=20, end=40, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr3_random", start=70, end=80, name="test", strand="."))
    genset1.standard_chromosome()
    assert len(genset1) == 2
def test_standard_chromosome2():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=50, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr2", start=20, end=40, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr3_random", start=70, end=80, name="test", strand="."))
    res = genset1.standard_chromosome(inplace=False)
    assert len(genset1) == 3
    assert len(res) == 2
def test_total_coverage1():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr3_random", start=1, end=80, name="test", strand="."))
    cov = genset1.total_coverage()
    assert cov == 85
def test_total_coverage2():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    cov = genset1.total_coverage()
    assert cov == 4
def test_rm_duplicates1():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr3_random", start=1, end=80, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test2", strand="."))
    res = genset1.rm_duplicates(inplace=False)
    assert len(res) == 3
def test_rm_duplicates2():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr3_random", start=1, end=80, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr2", start=2, end=4, name="test2", strand="-"))
    res = genset1.rm_duplicates(inplace=False)
    assert len(res) == 4

def test_distance1():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=20, end=24, name="test", strand="."))
    g = GenCoor(chrom="chr1", start=7, end=10, name="test", strand=".")
    res = genset1.distance(g, sign=False)
    assert res == 2

def test_distance2():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=20, end=24, name="test", strand="."))
    g = GenCoor(chrom="chr1", start=7, end=19, name="test", strand=".")
    res = genset1.distance(g, sign=False)
    assert res == 1

def test_distances1():
    genset1 = GenCoorSet(name="Test_set")
    genset1.add(GenCoor(chrom="chr1", start=1, end=5, name="test", strand="."))
    genset1.add(GenCoor(chrom="chr1", start=80, end=84, name="test", strand="."))
    genset2 = GenCoorSet(name="Test_set")
    genset2.add(GenCoor(chrom="chr1", start=10, end=15, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=60, end=64, name="test", strand="."))
    genset2.add(GenCoor(chrom="chr1", start=160, end=164, name="test", strand="."))
    res = genset1.distances(genset2, sign=False)
    print(res)
    assert len(res) == 2
    assert res == [5, 16]

os.system("bash <(curl -s https://codecov.io/bash) -t f34d6638-3e35-4782-98ac-8f355151c2fe")