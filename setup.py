#!/usr/bin/env python

import setuptools

VER = "0.0.2"

reqs = ["numpy",
        "scipy",
        "sympy",
        "datetime",
        "matplotlib",
        "particle",
        "h5py"
        ]

setuptools.setup(
    name="GAMPixTools",
    version=VER,
    author="Tom Shutt, Daniel Douglas, Henry Purcell and others",
    author_email="dougl215@slac.stanford.edu",
    description="A package for simulating the GAMPix readout for LArTPCs",
    packages=setuptools.find_packages(),
    install_requires=reqs,
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    python_requires='>=3.2',
)
