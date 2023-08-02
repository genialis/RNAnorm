====================================
RNA-seq data normalization in Python
====================================

|build| |black| |docs| |pypi_version| |pypi_pyversions| |pypi_downloads|

.. |build| image:: https://github.com/genialis/RNAnorm/actions/workflows/build.yml/badge.svg?branch=main
    :target: https://github.com/genialis/RNAnorm/actions?query=branch%3Amain
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


Python implementation of common RNA-seq normalization methods:

- CPM (Counts per million)
- FPKM_ (Fragments per kilobase million)
- TPM_ (Transcripts per million)
- UQ_ (Upper quartile)
- CUF_ (Counts adjusted with UQ factors)
- TMM_ (Trimmed mean of M-values)
- CTF_ (Counts adjusted with TMM factors)

For in-depth description of methods see documentation_.

.. _FPKM: https://www.nature.com/articles/nmeth.1226
.. _TPM: https://link.springer.com/article/10.1007/s12064-012-0162-3
.. _UQ: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-94
.. _CUF: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9/
.. _TMM: https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25
.. _CTF: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9/
.. _documentation: https://rnanorm.readthedocs.io/


Features
========

- Pure Python implementation (no need for R, etc.)
- Compatible with Scikit-learn_
- Command line interface
- Verbose documentation_
- Validated method implementation


.. _Scikit-learn: https://scikit-learn.org/
.. _documentation: https://rnanorm.readthedocs.io/


Install
=======

We recommend installing RNAnorm with pip::

   pip install rnanorm


Quick start
===========

The implemented methods can be executed from Python or from the command line.

Normalize from Python
---------------------

The most common use case is to run normalization from Python::

   >>> from rnanorm.datasets import load_toy_data
   >>> from rnanorm import FPKM
   >>> dataset = load_toy_data()
   >>> # Expressions need to have genes in columns and samples in rows
   >>> dataset.exp
             Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
   Sample_1     200     300     500    2000    7000
   Sample_2     400     600    1000    4000   14000
   Sample_3     200     300     500    2000   17000
   Sample_4     200     300     500    2000    2000
   >>> fpkm = FPKM(dataset.gtf_path).set_output(transform="pandas")
   >>> fpkm.fit_transform(dataset.exp)
                Gene_1    Gene_2    Gene_3    Gene_4    Gene_5
   Sample_1   100000.0  100000.0  100000.0  200000.0  700000.0
   Sample_2   100000.0  100000.0  100000.0  200000.0  700000.0
   Sample_3    50000.0   50000.0   50000.0  100000.0  850000.0
   Sample_4   200000.0  200000.0  200000.0  400000.0  400000.0


Normalize from command line
---------------------------

Normalization from the command line is also supported. To list available
methods and general help::

    rnanorm --help

Get info about a particular method, e.g., CPM::

    rnanorm cpm --help

To normalize with CPM::

   rnanorm cpm exp.csv --out exp_cpm.csv

File ``exp.csv`` needs to be comma separated file with genes in columns and
samples in rows. Values should be raw counts. The output is saved to
``exp_cpm.csv``. Example of input file::

    cat exp.csv
    ,Gene_1,Gene_2,Gene_3,Gene_4,Gene_5
    Sample_1,200,300,500,2000,7000
    Sample_2,400,600,1000,4000,14000
    Sample_3,200,300,500,2000,17000
    Sample_4,200,300,500,2000,2000

One can also provide input through standard input::

   cat exp.csv | rnanorm cpm --out exp_cpm.csv

If file specified with ``--out`` already exists the command will fail. If you
are sure that you wish to overwrite, use ``--force`` flag::

   cat exp.csv | rnanorm cpm --force --out exp_cpm.csv

If no file is specified with ``--out`` parameter, output is printed to standard
output::

   cat exp.csv | rnanorm cpm > exp_cpm.csv

Methods TPM and FPKM require gene lengths. These can be provided either with GTF_
file or with "gene lengths" file. The later is a two columns file. The first
column should include the genes in the header of ``exp.csv`` and the second
column should contain gene lengths computed by union exon model::

    # Use GTF file
    rnanorm tpm exp.csv --gtf annotations.gtf > exp_out.csv
    # Use gene lengths file
    rnanorm tpm exp.csv --gene-lengths lenghts.csv > exp_out.csv
    # Example of gene lengths file
    cat lenghts.csv
    gene_id,gene_length
    Gene_1,200
    Gene_2,300
    Gene_3,500
    Gene_4,1000
    Gene_5,1000


.. _GTF: https://www.ensembl.org/info/website/upload/gff.html


Contribute
==========

To learn about contributing to the code base, read the  Contributing_ section.

.. _Contributing: docs/contributing.rst

Citing
======

If you are using RNAnorm in your research, please cite as suggested by "Cite this repository" section in the side panel of this page.
