Ferro
===========

[![Documentation Status](https://readthedocs.org/projects/ferro/badge/?version=latest)](https://ferro.readthedocs.io/en/latest/?badge=latest)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FJAnderson419%2FFerro.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FJAnderson419%2FFerro?ref=badge_shield)

Master
[![Build Status](https://travis-ci.org/JAnderson419/Ferro.svg?branch=master)](https://travis-ci.org/JAnderson419/Ferro)
[![codecov](https://codecov.io/gh/JAnderson419/Ferro/branch/master/graph/badge.svg)](https://codecov.io/gh/JAnderson419/Ferro)

Dev
[![Build Status](https://travis-ci.org/JAnderson419/Ferro.svg?branch=dev)](https://travis-ci.org/JAnderson419/Ferro)
[![codecov](https://codecov.io/gh/JAnderson419/Ferro/branch/dev/graph/badge.svg)](https://codecov.io/gh/JAnderson419/Ferro)


Ferro aims to be a python package to ease manipulation of ferroelectric (and perhaps ferromagnetic) test data. It includes a HysteresisData class for read-in, storage, and display of PV/PUND/IV measurements as well as the beginning of a Landau modeling class (work in progress) that currently has a simple Presiach model implemented. If you would like to use some portion of this for your research (or add on features), I would be happy to help you get started.

Docstrings are in numpy style. For an example of a complete device analysis, see Ferro/bin/multidomainAnalysis.py

For an overview of the role of each function beyond the docstrings, see chapter 3 of "Measurement of Ferroelectric Films in MFM and MFIS Structures" by J. Anderson, available online: http://scholarworks.rit.edu/theses/9547/. Note that this work uses old syntax - the code has since been refactored to use more standard python syntax for classes and function/method names. Functionality, however, remains the same.

Developed on python 3.5 and 3.6 using Anaconda. If downloading from Github, be sure to use pip to install the package to ensure you have the proper prerequisites. If you have any questions, please contact me. 

**Jackson**


## License
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2FJAnderson419%2FFerro.svg?type=large)](https://app.fossa.com/projects/git%2Bgithub.com%2FJAnderson419%2FFerro?ref=badge_large)