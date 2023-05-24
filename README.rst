====================
RNAseq normalization
====================

|build| |black| |docs| |pypi_version| |pypi_pyversions| |pypi_downloads|

.. |build| image:: https://github.com/genialis/RNAnorm/workflows/build/badge.svg?branch=master
    :target: https://github.com/genialis/RNAnorm/actions?query=branch%3Amaster
    :alt: Build Status

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: Code Style Black

.. |docs| image:: https://readthedocs.org/projects/rnanorm/badge/?version=latest
    :target: http://rnanorm.readthedocs.io/
    :alt: Documentation Status

.. |pypi_version| image:: https://img.shields.io/pypi/v/rnanorm.svg
    :target: https://pypi.org/project/rnanorm
    :alt: Version on PyPI

.. |pypi_pyversions| image:: https://img.shields.io/pypi/pyversions/rnanorm.svg
    :target: https://pypi.org/project/rnanorm
    :alt: Supported Python versions

.. |pypi_downloads| image:: https://pepy.tech/badge/rnanorm
    :target: https://pepy.tech/project/rnanorm
    :alt: Number of downloads from PyPI


Python implementation of common RNAseq normalization methods:

- CPM (Counts per million)
- FPKM_ (Fragments per kilobase million)
- TPM_ (Transcripts per million)
- UQ_ (Upper quartile)
- CUF_ (Counts adjusted with UQ factors)
- TMM_ (Trimmed mean of M-values)
- CTF_ (Counts adjusted with TMM factors)


.. _FPKM: https://www.nature.com/articles/nmeth.1226
.. _TPM: https://link.springer.com/article/10.1007/s12064-012-0162-3
.. _UQ: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94
.. _CUF: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9/
.. _TMM: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25
.. _CTF: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9/

Features
========

- Pure Python implementation (no need for R, etc.)
- Scikit-learn_ compatible
- Command line interface
- Verbose documentation_ (at least we hope so...)
- Tested methods


.. _Scikit-learn: https://scikit-learn.org/
.. _documentation: https://rnanorm.readthedocs.io/


Install
=======

We recommend installing RNAnorm with pip::

   pip install rnanorm


Quick start
===========

Implemented methods can be used from Python or from the command line.

Normalize from Python
---------------------

Most commonly normalization methods are run from Python. E.g.::

   >>> from rnanorm.datasets import load_rnaseq_toy
   >>> from rnanorm import FPKM
   >>> dataset = load_rnaseq_toy()
   >>> dataset.exp
          G1     G2      G3      G4       G5
   S1  200.0  300.0   500.0  2000.0   7000.0
   S2  400.0  600.0  1000.0  4000.0  14000.0
   S3  200.0  300.0   500.0  2000.0  17000.0
   S4  200.0  300.0   500.0  2000.0   2000.0
   >>> fpkm = FPKM(dataset.gtf_path).set_output(transform="pandas")
   >>> fpkm.fit_transform(dataset.exp)
             G1        G2        G3        G4        G5
   S1  100000.0  100000.0  100000.0  200000.0  700000.0
   S2  100000.0  100000.0  100000.0  200000.0  700000.0
   S3   50000.0   50000.0   50000.0  100000.0  850000.0
   S4  200000.0  200000.0  200000.0  400000.0  400000.0


Normalize from command line
---------------------------

Often it is handy to do normalization from the command line::

   rnanorm fpkm exp.csv --gtf annotation.gtf --out exp_fpkm.csv


Contribute
==========

To learn about contributing to the code base, read the  Contributing_ section.

.. _Contributing: docs/contributing.rst
