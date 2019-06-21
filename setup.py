import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Gen-Coor",
    version="0.0.1",
    author="Joseph Chao-Chung Kuo",
    author_email="jovesus@gmail.com",
    description="A package for storing and calculating genomic coordinates.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jovesus/gen-coor",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
