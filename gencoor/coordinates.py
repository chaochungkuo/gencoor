from .exceptions import ChromosomeNotStrError, PositionNotIntegerError, NameNotStrError, \
                        StrandNotStrError, CoordinateFlipError


class GenCoor:
    """A Python class for handling a genomic coordinate.

    Parameters
    ----------
    chrom : string
    Chromosome name of the genomic coordinate
    start : int
    Start position of the genomic coordinate
    end : int
    End position of the genomic coordinate
    name : string
    Name of the genomic coordinate
    strand : string
    Strand of the genomic coordinate, '+', '-' or '.'
    score : float
    Score for any numerical value, such as fold change or p-value
    data : list
    Data stores a list including any extra information of the genomic coordinate

    Examples
    --------
    >>> g = GenCoor(chrom='chr1', start=100, end=200, name='region_1', strand='+')
    >>> len(g)


    """

    def __init__(self, chrom, start, end, name="", strand=".", score=None, data=None):

        if not isinstance(chrom, str):
            print(f"chromosome is not a string: {chrom}")
            raise ChromosomeNotStrError
        if not isinstance(start, int):
            print(f"start is not a integer: {start}")
            raise PositionNotIntegerError
        if not isinstance(end, int):
            print(f"end is not a integer: {end}")
            raise PositionNotIntegerError
        if not isinstance(name, str):
            print(f"name is not a string: {name}")
            raise NameNotStrError
        if not isinstance(strand, str):
            print(f"strand is not a string: {name}")
            raise StrandNotStrError
        if start > end:
            print(f"start position ({start}) is larger than end position ({end})")
            raise CoordinateFlipError

        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.strand = strand
        self.score = score
        self.data = data

    def __len__(self):
        """Return the length of the coordinate."""
        return self.end - self.start

    def __str__(self):
        """Return a readable string."""
        return f"{str(self.start)}-{str(self.end)}"

    def ___repr__(self):
        """Return all information as a string."""
        return f"{self.chrom}:{str(self.start)}-{str(self.end)} {self.strand}"

    def capital_name(self):
        """Capitalize the name of the coordinate."""
        self.name = self.name.upper()

    def overlap(self, a_gencoor, strandness=False):
        """Return True if it overlaps with the given genomic coordinate, otherwise False."""
        res = False
        # region against region
        if self.chrom == a_gencoor.chrom and \
                self.start < a_gencoor.end and \
                self.end > a_gencoor.start:
            if strandness:
                if self.strand == a_gencoor.strand:
                    res = True
            else:
                res = True
        # single point to single point
        elif self.chrom == a_gencoor.chrom and \
                self.start == a_gencoor.end and \
                self.end == a_gencoor.start:
            res = True
        return res

    def distance(self, a_gencoor, sign=False):
        """Return the distance between two genomic coordinates. '0' means overlapping and None means they are \
        on different chromosomes and their distance cannot be calculated. If sign is activated, \
        positive values mean that the given query coordinate is on the downstream of the reference coordinate. \
        Negative values mean opposite, which is upstream."""

        if self.overlap(a_gencoor):
            return 0
        else:
            if self.chrom == a_gencoor.chrom:
                if self.start < a_gencoor.start:
                    dis = a_gencoor.start - self.end
                else:
                    dis = a_gencoor.end - self.start
            else:
                return None

        if sign:
            return dis
        else:
            return abs(dis)


class GenCoorSet:
    """A Python class for handling a set of genomic coordinates.

        Parameters
        ----------
        name : string
        Name of this set

        Attributes
        ----------
        list : list
        A list containing all the genomic coordinates
        sort : boolean
        A label showing whether this list is sorted or not

        Examples
        --------
        >>> regions = GenCoorSet(name="A_set")
        >>> len(regions)


        """

    def __init__(self, name):
        self.list = []
        self.sort = False
        self.name = name

    def add(self, gencoor):
        """Add a genomic coordinate into the list"""
        self.list.append(gencoor)
        self.sort = False

    def len(self):
        """Return the length of the list of genomic coordinates"""
        return len(self.list)

