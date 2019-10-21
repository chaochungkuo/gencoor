from .exceptions import OverlapTypeError

class GenCoorFileIO:
    class Bed:
        """
        Each row maps to a GenCoor.

        Note: Chrom (1), start (2), end (2), name (4) and orientation (6) is used for GenCoor.
              All other columns (5, 7, 8, ...) are put to the data attribute of the GenCoor.
              The numbers in parentheses are the columns of the BED format.
        """

        @staticmethod
        def read_to_gcs(gcs, str filename):
            cdef str line
            cdef list l
            with open(filename, "r") as f:
                for line in f.readlines():
                    line = line.strip()
                    if not line.startswith("#"):
                        l = line.split()
                        gcs.add(GenCoor(chrom=l[0], start=int(l[1]), end=int(l[2]),
                                        name=l[3], score=l[4], strand=l[5], data="/t".join(l[6:])))

        @staticmethod
        def write_from_gcs(gcs, str filename, str mode="w"):
            with open(filename, mode) as f:
                for gc in gcs:
                    print(gc, file=f)

    class Bed12:
        """
        Bed file with "block information", eg exons.
        """

        @staticmethod
        def read_to_gcs(grs, str filename):
            cdef str line
            cdef list l
            with open(filename, "r") as f:
                for line in f.readlines():
                    line = line.strip()
                    if not line.startswith("#"):
                        l = line.split()
                        gencoor = GenCoor(chrom=l[0], start=int(l[1]), end=int(l[2]),
                                          name=l[3], score=l[4], strand=l[5], data="/t".join(l[6:]))
                        gencoor.extract_bed12()
                        grs.add(gencoor)

        @staticmethod
        def write_from_gcs(gcs, filename, mode="w"):
            with open(filename, mode) as f:
                for gc in gcs:
                    print(gc, file=f)
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

    def __init__(self, str chrom, int start, int end, str name="",
                 str strand=".", score=0, str data=""):

        # if not isinstance(chrom, str):
        #     print(f"chromosome is not a string: {chrom}")
        #     raise ChromosomeNotStrError
        # if not isinstance(start, int):
        #     print(f"start is not a integer: {start}")
        #     raise PositionNotIntegerError
        # if not isinstance(end, int):
        #     print(f"end is not a integer: {end}")
        #     raise PositionNotIntegerError
        # if not isinstance(name, str):
        #     print(f"name is not a string: {name}")
        #     raise NameNotStrError
        # if not isinstance(strand, str):
        #     print(f"strand is not a string: {name}")
        #     raise StrandNotStrError
        # if start > end:
        #     print(f"start position ({start}) is larger than end position ({end})")
        #     raise CoordinateFlipError

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

    def __hash__(self):
        return hash((self.chrom, self.start, self.end, self.strand))

    def __eq__(self, other):
        return (self.chrom, self.start, self.end, self.strand) == \
               (other.chrom, other.start, other.end, other.strand)

    def __str__(self):
        """Return a readable string."""
        # return f"{str(self.start)}-{str(self.end)}"
        if not self.data:
            return("\t".join([self.chrom, str(self.start), str(self.end),
                              self.name, str(self.score), self.strand]))
        else:
            return ("\t".join([self.chrom, str(self.start), str(self.end),
                               self.name, str(self.score), self.strand,
                               self.data]))

    def __repr__(self):
        """Return all information as a string."""
        return f"{self.chrom}:{str(self.start)}-{str(self.end)} {self.strand}"

    def __cmp__(self, region):
       """Return negative value if x < y, zero if x == y and strictly positive if x > y."""
       if self.chrom < region.chrom:
           return -1
       elif self.chrom > region.chrom:
           return 1
       else:
           if self.start < region.start:
               return -1
           elif self.start > region.start:
               return 1
           else:
               if self.end < region.end:
                   return -1
               elif self.end > region.end:
                   return 1
               else:
                   return 0

    def __lt__(self, other):
        return self.__cmp__(other) < 0

    def __gt__(self, other):
        return self.__cmp__(other) > 0

    def __eq__(self, other):
        return self.__cmp__(other) == 0

    def __le__(self, other):
        return self.__cmp__(other) <= 0

    def __ge__(self, other):
        return self.__cmp__(other) >= 0

    def __ne__(self, other):
        return self.__cmp__(other) != 0

    def str_split(self, filetype="BED"):
        if filetype=="BED":
            return("/t".join([self.chrom, str(self.start), str(self.end),
                              self.name, str(self.score), self.strand]))
        else:
            return ("/t".join([self.chrom, str(self.start), str(self.end),
                               self.name, str(self.score), self.strand,
                               self.data]))

    def capital_name(self):
        """Capitalize the name of the coordinate."""
        self.name = self.name.upper()

    def overlap(self, region, strand_specific=False):
        """Return True if it overlaps with the given genomic coordinate, otherwise False."""
        # return self.overlap_pairing(self.chrom, self.start, self.end, self.strand,
        #                             region.chrom, region.start, region.end, region.strand,
        #                             strandness)
        res = False
        # region against region
        if self.chrom == region.chrom and \
                self.start < region.end and \
                self.end > region.start:
            if strand_specific:
                if self.strand == region.strand:
                    res = True
            else:
                res = True
        # single point to single point
        elif self.chrom == region.chrom and \
                self.start == region.end and \
                self.end == region.start:
            res = True
        return res

    @staticmethod
    def overlap_pairing(str chrom_1, int start_1, int end_1, str strand_1,
                        str chrom_2, int start_2, int end_2, str strand_2,
                        strand_specific=False):
        """Return True if it overlaps with the given genomic coordinate, otherwise False."""
        res = False
        # region against region
        if chrom_1 == chrom_2 and start_1 < end_2 and end_1 > start_2:
            if strand_specific:
                if strand_1 == strand_2:
                    res = True
            else:
                res = True
        # single point to single point
        elif chrom_1 == chrom_2 and start_1 == end_2 and end_1 == start_2:
            res = True
        return res

    def distance(self, region, sign=False):
        """Return the distance between two genomic coordinates. '0' means overlapping and None means they are \
        on different chromosomes and their distance cannot be calculated. If sign is activated, \
        positive values mean that the given query coordinate is on the downstream of the reference coordinate. \
        Negative values mean opposite, which is upstream."""
        cdef int dis
        if self.overlap(region):
            return 0
        else:
            if self.chrom == region.chrom:
                if self.start < region.start:
                    dis = region.start - self.end
                else:
                    dis = region.end - self.start
            else:
                return None

        if sign:
            return dis
        else:
            return abs(dis)

    def extract_bed12(self):
        """Return a list which contains all the genomic coordinates extracted from the given \
        region in BED12 format."""
        l = self.data.split("/t")
        blocksizes = [int(x.strip()) for x in l[4].split(",")]
        blockstarts = [int(x.strip()) for x in l[5].split(",")]
        res = []
        cdef int i
        for i in range(int(l[3])):
            res.append(GenCoor(chrom=self.chrom,
                               start=self.start + blockstarts[i],
                               end=self.start + blockstarts[i] + blocksizes[i],
                               name=self.name,
                               strand=self.strand))
        return res


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

    def __init__(self, str name):
        self.list = []
        self.sorted = False
        self.name = name

    def __len__(self):
        """Return the length of the list of genomic coordinates"""
        return len(self.list)

    def __iter__(self):
        return iter(self.list)

    def __getitem__(self, key):
        return(self.list[key])

    def add(self, gencoor):
        """Add a genomic coordinate into the list"""
        self.list.append(gencoor)
        self.sorted = False

    def sort(self, inplace=True):
        """Sort the genomic coordinates by their genomic position in the order of chromosome, \
        start position and then end position. If inplace is set True, it sorts the genomic coordinates in place, \
        otherwise, it returns a sorted genomic coordinate set."""
        res = sorted(self.list, key=lambda x: x.start)
        res = sorted(res, key=lambda x: x.chrom)
        if inplace:
            self.list = res
        else:
            res_gencoorset = GenCoorSet(name=self.name)
            res_gencoorset.sorted = True
            res_gencoorset.list = res
            return res_gencoorset

    def load(self, str filename, filetype="BED"):
        """Load the genomic coordinates from the given file. The options for file format are: BED and BED12."""
        self.sorted = False
        if filetype == "BED":
            self = GenCoorFileIO.Bed.read_to_gcs(self, filename)
        elif filetype == "BED12":
            self = GenCoorFileIO.Bed12.read_to_gcs(self, filename)

    def save(self, str filename, filetype="BED"):
        """Save the genomic coordinates to a file. The options for file format are: BED and BED12."""

        if filetype == "BED":
            GenCoorFileIO.Bed.write_from_gcs(self, filename)
        elif filetype == "BED12":
            GenCoorFileIO.Bed12.write_from_gcs(self, filename)

    def extend(self, str mode, int length, inplace=False):
        """Extend the genomic coordinates by certain length according to the given mode.

           Mode
           ----------
           left
           right
           5end
           3end
           both

        """
        new_gcs = GenCoorSet(self.name)

        def ext_start(gc):
            ngc = GenCoor(chrom=gc.chrom, name=gc.name, strand=gc.strand, score=gc.score,
                          start=gc.start - length, end=gc.end, data=gc.data)
            return ngc

        def ext_end(gc):
            ngc = GenCoor(chrom=gc.chrom, name=gc.name, strand=gc.strand, score=gc.score,
                          start=gc.start, end=gc.end + length, data=gc.data)
            return ngc

        def ext_both(gc):
            ngc = GenCoor(chrom=gc.chrom, name=gc.name, strand=gc.strand, score=gc.score,
                          start=gc.start - length, end=gc.end + length, data=gc.data)
            return ngc

        for gc in self:
            if mode == "left":
                ngc = ext_start(gc)
            elif mode == "right":
                ngc = ext_end(gc)
            elif mode == "5end":
                if gc.strand == "+":
                    ngc = ext_start(gc)
                elif gc.strand == "-":
                    ngc = ext_end(gc)
            elif mode == "3end":
                if gc.strand == "+":
                    ngc = ext_end(gc)
                elif gc.strand == "-":
                    ngc = ext_start(gc)
            else:
                ngc = ext_both(gc)

            new_gcs.add(ngc)

        if inplace:
            self.list = new_gcs.list
        else:
            return new_gcs

    def to_lists(self):
        cdef list chroms=[]
        cdef list starts=[]
        cdef list ends=[]
        cdef list names=[]
        cdef list scores=[]
        cdef list strands=[]
        cdef list datas=[]
        for r in self:
            chroms.append(r.chrom)
            starts.append(r.start)
            ends.append(r.end)
            names.append(r.name)
            scores.append(r.score)
            strands.append(r.strand)
            datas.append(r.data)
        return chroms, starts, ends, names, scores, strands, datas

    def split_by_strands(self):
        """Return a dictionary with unique strand as the keys and GenCoorSet as the values."""
        res = {}
        for r in self:
            if r.strand not in res.keys():
                res[r.strand] = GenCoorSet(name=r.strand)
                res[r.strand].add(r)
            else:
                res[r.strand].add(r)
        return res

    def split_by_chromosome(self):
        """Return a dictionary with unique chromosome as the keys and GenCoorSet as the values."""
        res = {}
        for r in self:
            if r.chrom not in res.keys():
                res[r.chrom] = GenCoorSet(name=r.chrom)
                res[r.chrom].add(r)
            else:
                res[r.chrom].add(r)
        return res

    def merge(self, w_return=False, strand_specific=False):
        """Merge the regions within the GenCoorSet

            Parameters
            ----------
            w_return : If TRUE, it returns a GenCoorSet; if FALSE, it merges the regions in place.
            strand_specific : If TRUE, only the regions on the same strand will be merged.
        """
        def add_pre_r(z, pre_r):
            z.add(GenCoor(chrom=pre_r[0],
                          start=pre_r[1],
                          end=pre_r[2],
                          name=pre_r[4],
                          strand=pre_r[3]))

        cdef int i
        cdef int len_self = len(self)

        if not self.sorted:
            self.sort()
        if len(self) in [0, 1]:
            if w_return:
                return self
            else:
                return
        else:
            z = GenCoorSet(name=self.name)
            chroms, starts, ends, names, scores, strands, datas = self.to_lists()
            pre_r = None
            if not strand_specific:
                ########################################
                ## Merging without strand specific
                for i in range(0, len_self):
                    if not pre_r:
                        pre_r = [chroms[i], starts[i], ends[i],
                                       strands[i], names[i]]
                        continue
                    if GenCoor.overlap_pairing(pre_r[0],
                                               pre_r[1],
                                               pre_r[2],
                                               pre_r[3],
                                               chroms[i],
                                               starts[i],
                                               ends[i],
                                               strands[i]):
                        pre_r[1] = min(pre_r[1], starts[i])
                        pre_r[2] = max(pre_r[2], ends[i])
                        pre_r[4] = ",".join([pre_r[4], names[i]])
                        if not pre_r[3] == strands[i]:
                            pre_r[3] = "."
                        if i == len_self-1:
                            add_pre_r(z, pre_r)
                    else:
                        add_pre_r(z, pre_r)
                        if i == len_self-1:
                            z.add(GenCoor(chrom=chroms[i],
                                          start=starts[i],
                                          end=ends[i],
                                          name=names[i],
                                          strand=strands[i]))
            else:
                ########################################
                ## Merging with strand specific
                strand_type = list(set([g.strand for g in self]))
                pre_r = {}
                for s in strand_type:
                    pre_r[s] = None
                for i in range(0, len_self):
                    if not pre_r[strands[i]]:
                        pre_r[strands[i]] = [chroms[i], starts[i], ends[i], strands[i], names[i]]
                    if GenCoor.overlap_pairing(pre_r[strands[i]][0],
                                               pre_r[strands[i]][1],
                                               pre_r[strands[i]][2],
                                               pre_r[strands[i]][3],
                                               chroms[i], starts[i], ends[i], strands[i],
                                               strand_specific=strand_specific):
                        pre_r[strands[i]][1] = min(pre_r[strands[i]][1], starts[i])
                        pre_r[strands[i]][2] = max(pre_r[strands[i]][2], ends[i])
                        pre_r[strands[i]][4] = ",".join([pre_r[strands[i]][4], names[i]])
                    else:
                        add_pre_r(z, pre_r[strands[i]])
                        pre_r[strands[i]] = [chroms[i], starts[i], ends[i], strands[i], names[i]]
                    if i == len_self-1:
                        for s in strand_type:
                            add_pre_r(z, pre_r[s])
            z.sort()
            if w_return:
                return z
            else:
                self.list = z

    def intersect(self, gcs, mode="overlap", strand_specific=False):
        """Get intersection regions between two GenCoorSets

            Parameters
            ----------
            mode : ['overlap', 'original', 'complete_included']
            strand_specific : If TRUE, only the regions on the same strand will be calcuated.
        """
        cdef int i = 0
        cdef int j = 0
        cdef int prej = 0
        cdef int last_i = len(self) - 1
        cdef int last_j = len(gcs) - 1
        cdef int cont_loop = 1
        cdef int cont_overlap = 0

        def check_i_end(i, last_i, cont_loop):
            if i < last_i:
                i += 1
            else:
                cont_loop = 0
            return i, cont_loop

        z = GenCoorSet(self.name)
        if len(self) == 0 or len(gcs) == 0:
            return z
        else:
            if not self.sorted:
                self.sort()
            if not gcs.sorted:
                gcs.sort()
            chroms1, starts1, ends1, names1, scores1, strands1, datas1 = self.to_lists()
            chroms2, starts2, ends2, names2, scores2, strands2, datas2 = gcs.to_lists()
            ####################### OverlapType.OVERLAP ###############################
            if mode == "overlap":
                while cont_loop == 1:
                    # When the regions overlap
                    if GenCoor.overlap_pairing(chroms1[i], starts1[i], ends1[i], strands1[i],
                                               chroms2[j], starts2[j], ends2[j], strands2[j],
                                               strand_specific=strand_specific):
                        z.add(GenCoor(chrom=chroms1[i],
                                      start=max(starts1[i], starts2[j]),
                                      end=min(ends1[i], ends2[j]),
                                      name=names1[i],
                                      strand=strands1[i],
                                      data=datas1[i]))
                        if cont_overlap == 0:
                            prej = j
                        if j == last_j:
                            i, cont_loop = check_i_end(i, last_i, cont_loop)
                        else:
                            j += 1
                        cont_overlap = 1

                    elif self[i] < gcs[j]:
                        if i < last_i:
                            i += 1
                            if chroms1[i] == chroms2[j]:
                                j = prej
                            cont_overlap = 0
                        else:
                            cont_loop = 0
                    elif self[i] > gcs[j]:
                        if j == last_j:
                            cont_loop = 0
                        else:
                            j += 1
                            cont_overlap = 0
                    else:
                        i, cont_loop = check_i_end(i, last_i, cont_loop)

            ####################### OverlapType.ORIGINAL ###############################
            elif mode == "original":
                while cont_loop:
                    # When the regions overlap
                    if GenCoor.overlap_pairing(chroms1[i], starts1[i], ends1[i], strands1[i],
                                               chroms2[j], starts2[j], ends2[j], strands2[j],
                                               strand_specific=strand_specific):
                        z.add(self[i])
                        i, cont_loop = check_i_end(i, last_i, cont_loop)
                    elif self[i] < gcs[j]:
                        i, cont_loop = check_i_end(i, last_i, cont_loop)
                    elif self[i] > gcs[j]:
                        if j == last_j:
                            cont_loop = 0
                        else:
                            j += 1
                    else:
                        i, cont_loop = check_i_end(i, last_i, cont_loop)
            ####################### OverlapType.COMP_INCL ###############################
            elif mode == "complete_included":
                while cont_loop:
                    # When the regions overlap
                    if GenCoor.overlap_pairing(chroms1[i], starts1[i], ends1[i], strands1[i],
                                               chroms2[j], starts2[j], ends2[j], strands2[j],
                                               strand_specific=strand_specific):
                        if starts1[i] >= starts2[j] and ends1[i] <= ends2[j]:
                            z.add(self[i])
                        if cont_overlap==0:
                            prej = j
                        if j == last_j:
                            i, cont_loop = check_i_end(i, last_i, cont_loop)
                        else:
                            j += 1
                        cont_overlap = 1

                    elif self[i] < gcs[j]:
                        if i < last_i:
                            if chroms1[i] == chroms2[j] and prej > 0:
                                j = prej
                            i += 1
                            cont_overlap = 0
                        else:
                            cont_loop = 0
                    elif self[i] > gcs[j]:
                        if j == last_j:
                            cont_loop = 0
                        else:
                            j += 1
                            cont_overlap = 0
                    else:
                        i, cont_loop = check_i_end(i, last_i, cont_loop)
            else:
                print("Please define the mode as one of the three options: overlap, original, or complete_included.")
                raise OverlapTypeError
            return z

    # def standard_chromosome(self, organism):
    #
    # def total_coverage(self):
    #
    # def rm_duplicates(self):
    #
    # def
