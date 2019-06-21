import sys, os
myPath = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, myPath + '/../')

import pytest
from gencoor import GenCoor


def test_GenCoor_init():
    # Normal case
    g1 = GenCoor(chrom="chr1", start=1, end=100, name="g1", strand=".", score=None, data=None)
    with pytest.raises(Exception) as e_info:
        g1 = GenCoor(chrom=2, start=1, end=100, name="g1", strand=".", score=None, data=None)
    assert e_info.value.message == 'chromosome is not a string: 2'
    # integer to chrom
    # with pytest.raises(Exception) as e_info:

