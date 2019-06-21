# define Python user-defined exceptions
class Error(Exception):
   """Base class for other exceptions"""
   pass

class CoordinateFlipError(Error):
   """Raised when the coordinate is flipped, e.g. start position > end position"""
   pass

class ValueTooLargeError(Error):
   """Raised when the input value is too large"""
   pass