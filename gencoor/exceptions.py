# define Python user-defined exceptions
class Error(Exception):
    """Base class for other exceptions"""
    pass


class CoordinateFlipError(Error):
    """Raised when the coordinate is flipped, e.g. start position > end position"""
    pass


class NegativePositionError(Error):
    """Raised when the position on the genome is negative"""
    pass


class ChromosomeNotStrError(Error):
    """Raised when the input chromosome is not a string"""
    pass


class PositionNotIntegerError(Error):
    """Raised when the input position is not an integer"""
    pass


class NameNotStrError(Error):
    """Raised when the input name is not a string"""
    pass


class StrandNotStrError(Error):
    """Raised when the input strand is not a string"""
    pass
