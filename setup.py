import setuptools
import io
import os
import re
import shutil
from sys import platform
from Cython.Build import cythonize
from setuptools.command.install import install

"""
Installs gencoor with standard setuptools options.

Authors: Joseph Chao-Chung Kuo
"""

# if the environment variable is set, use it; otherwise use the home directory as a default
data_path = os.path.expanduser(os.getenv("GENCOORDATA", os.path.join(os.getenv("HOME"), "gencoor_data")))
# Creating Data Path
if not os.path.exists(data_path):
    os.makedirs(data_path)

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

class CustomInstallCommand(install):
    """Customized setuptools install command."""

    def write_genome_data(self, config, genome, latest_GTFs, genome_organism):
        config.write("[" + genome + "]\n")
        config.write("genome: " + os.path.join(data_path, genome, "genome_" + genome + ".fa\n"))
        config.write("chromosome_sizes: " + os.path.join(data_path, genome, "chrom.sizes." + genome + "\n"))
        config.write("genes: " + os.path.join(data_path, genome, "genes_" + genome + ".bed\n"))
        config.write("annotation: " + os.path.join(data_path, genome, latest_GTFs[genome] + "\n"))
        config.write("gene_alias: " + os.path.join(data_path, genome, "alias_" + genome_organism[genome] + ".txt\n\n"))

    def run(self):
        ###################################################################################################
        # Creating Data Path
        ###################################################################################################
        latest_GTFs = {"mm9": "gencode.vM1.annotation.gtf",
                       "mm10": "gencode.vM23.annotation.gtf",
                       "hg19": "gencode.v19.annotation.gtf",
                       "hg38": "gencode.v32.annotation.gtf"}
        genome_organism = {"mm9": "mouse",
                           "mm10": "mouse",
                           "hg19": "human",
                           "hg38": "human"}

        # if the environment variable is set, use it; otherwise use the home directory as a default
        self.data_path = os.path.expanduser(
            os.getenv("GENCOORDATA", os.path.join(os.getenv("HOME"), "gencoor_data")))

        # Creating data.config
        data_config_file_name = os.path.join(self.data_path, "data.config")
        data_config_file = open(data_config_file_name, "w")
        data_config_file.write(
            "# Configuration file loaded at rgt startup. CAREFUL: any changes shall be overwritten\n"
            "# whenever rgt is (re)installed. Use data.config.user for permanent changes.\n\n")
        for genome in latest_GTFs.keys():
            self.write_genome_data(data_config_file, genome, latest_GTFs, genome_organism)
        data_config_file.close()

        # Creating data.config.user, but only if not already present
        user_config_file_name = os.path.join(self.data_path, "data.config.user")
        if not os.path.isfile(user_config_file_name):
            user_config_file = open(user_config_file_name, "w")

            user_config_file.write(
                "# Here you can overwrite any property set in the data.config file. It shall not be\n"
                "# be overwritten in any case, so if you are experiencing problems rename or remove this\n"
                "# file. See data.config for how the file should be formatted.\n\n")
            genome = "self_defined"
            user_config_file.write("# Template to add a genomic section.\n")
            user_config_file.write("#[" + genome + "]\n")
            user_config_file.write("#genome: undefined\n")
            user_config_file.write("#chromosome_sizes: undefined\n")
            user_config_file.write("#gene_regions: undefined\n")
            user_config_file.write("#annotation: undefined\n")
            user_config_file.write("#gene_alias: undefined\n\n")
        install.run(self)

current_version = find_version("gencoor", "__version__.py")

###################################################################################################
# Unsupported Platforms
###################################################################################################

supported_platforms = ["linux", "linux2", "darwin"]
if platform not in supported_platforms:
    print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
    exit(0)

###################################################################################################
# Copying the data files
###################################################################################################
data_dir = os.path.join(os.path.dirname(__file__), "data")

shutil.copy(os.path.join(data_dir, "setup_genome_data.py"),
            os.path.join(data_path, "setup_genome_data.py"))
for entry in os.listdir(data_dir):
    genome_dir = os.path.join(data_dir, entry)
    dest_genome_dir = os.path.join(data_path, entry)

    if os.path.isdir(genome_dir):
        if not os.path.exists(dest_genome_dir):
            os.makedirs(dest_genome_dir)
        for file in os.listdir(genome_dir):
            source_f = os.path.join(data_dir, entry, file)
            destin_f = os.path.join(data_path, entry, file)
            shutil.copy(source_f, destin_f)
#############################

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
                          annotate=True),
    cmdclass={'install': CustomInstallCommand}
)
