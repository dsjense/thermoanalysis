#!/usr/bin/env python3

from setuptools import setup, find_packages
import sys

if sys.version_info.major < 3:
    raise SystemExit("Python 3 is required!")

setup(
    name="thermoanalysis",
    version="0.1",
    description="Thermochemistry from QC program logs with cclib.",
    url="https://github.com/eljost/thermoanalysis",
    maintainer="Johannes Steinmetzer",
    maintainer_email="johannes.steinmetzer@uni-jena.de",
    license="GPL 3",
    platforms=["unix"],
    packages=find_packages(),
    install_requires=[
        "cclib",
        "numpy",
        "h5py",
        "pandas",
        "pytest",
        "tabulate",
    ],
    entry_points={
        "console_scripts": [
            "thermo = thermoanalysis.main:run",
        ]
    },
)
