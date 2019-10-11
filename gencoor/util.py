class OverlapType:
    """Class of overlap type constants.

    *Constants:*

        - OVERLAP -- Return new GenomicRegionSet including only the overlapping regions.
        - ORIGINAL -- Return the regions of original GenomicRegionSet which have any intersections.
        - COMP_INCL -- Return region(s) of the GenomicRegionSet which are 'completely' included.
    """

    OVERLAP = 0
    ORIGINAL = 1
    COMP_INCL = 2