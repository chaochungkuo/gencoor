

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
        assert isinstance(chrom, str), f"chromosome is not a string: {chrom}"
        assert isinstance(start, int), f"start is not a integer: {start}"
        assert isinstance(end, int), f"start is not a integer: {end}"
        assert isinstance(name, str), f"chromosome is not a string: {name}"
        assert isinstance(strand, str), f"chromosome is not a string: {strand}"
        assert start <= end, f"start position ({start}) is larger than end position ({end})"

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

    # """A Python class for storing and calculating genomic coordinates"""
    # def __init__(self, name):
    #     self.name = name
    #     self.coordinates = None
    #     self.sorted = False
    #     # self.coordinates = np.array([[1, 12]], dtype=object)
    #
    # def load_bed(self, filepath):
    #     self.sorted = False
    #
    #     coordinates = []
    #     with open(filepath) as f:
    #         for line in f:
    #             if not line.startswith("#"):
    #                 l = line.strip().split()
    #                 if len(l) == 6:
    #                     coordinates.append([l[0], int(l[1]), int(l[2]), l[3],
    #                                         float(l[4]), l[5]])
    #                 elif len(l) < 4:
    #                     coordinates.append([l[0], int(l[1]), int(l[2]), l[3],
    #                                         float(l[4]), l[5]])
    #
    #     self.coordinates = np.array(coordinates, dtype=object)
    #
    # def sort(self):
    #     self.sorted = True
    #     self.coordinates = self.coordinates[self.coordinates[:, 1].argsort()]
    #     self.coordinates = self.coordinates[self.coordinates[:, 0].argsort()]
    #
    # def add(self, coordinate):
    #     coordinate[1] = int(coordinate[1])
    #     coordinate[2] = int(coordinate[2])
    #     coordinate[4] = float(coordinate[4])
    #     self.coordinates = np.vstack((self.coordinates, np.array(coordinate, dtype=object)))
    #     self.sorted = False
    #
    # # def write(self, filepath):
    # #
    # def __len__(self):
    #     return self.coordinates.shape[0]



