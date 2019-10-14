import setuptools
import io
import os
import re
from sys import platform
from Cython.Build import cythonize

"""
Installs gencoor with standard setuptools options.

Authors: Joseph Chao-Chung Kuo
"""

def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

current_version = find_version("gencoor", "__version__.py")

###################################################################################################
# Unsupported Platforms
###################################################################################################

supported_platforms = ["linux", "linux2", "darwin"]
if platform not in supported_platforms:
    print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
    exit(0)


with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gencoor",
    version="0.0.1",
    author="Joseph Chao-Chung Kuo",
    author_email="jovesus@gmail.com",
    description="A package for storing and analyzing genomic coordinates.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jovesus/gencoor",
    packages=setuptools.find_packages(),
    classifiers=["Programming Language :: Python :: 3",
                 "License :: OSI Approved :: MIT License",
                 "Operating System :: OS Independent"],
    ext_modules=cythonize("gencoor/coordinates.pyx",
                          compiler_directives={'language_level': "3"},
                          annotate=True)
)
