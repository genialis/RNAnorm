=====================
RNA-Seq Normalization
=====================

|build| |black| |pypi_version| |pypi_pyversions| |pypi_downloads|

.. |build| image:: https://github.com/genialis/rnaseq-normalization/workflows/build/badge.svg?branch=master
    :target: https://github.com/genialis/rnaseq-normalization/actions?query=branch%3Amaster
    :alt: Build Status

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: Code Style Black

.. |pypi_version| image:: https://img.shields.io/pypi/v/rnanorm.svg
    :target: https://pypi.org/project/rnanorm
    :alt: Version on PyPI

.. |pypi_pyversions| image:: https://img.shields.io/pypi/pyversions/rnanorm.svg
    :target: https://pypi.org/project/rnanorm
    :alt: Supported Python versions

.. |pypi_downloads| image:: https://pepy.tech/badge/rnanorm
    :target: https://pepy.tech/project/rnanorm
    :alt: Number of downloads from PyPI

Normalization of RNA-seq gene expression data. Supported methods:

* Counts per million (CPM)
* Transcript per kilobase million (TPM)
* Fragments per kilobase million (FPKM)
* Quantile normalization to average distribution

The TPM / FPKM normalization can either accept pre-computed gene lengths
on the input or compute gene lengths from gene annotation in GTF format,
using the union exon-based approach. The computed gene lengths are
identical to the lengths reported by featureCounts (validated for *Homo
sapiens*, *Mus musculus*, *Rattus norvegicus*, and *Macaca mulatta* of
ENSEMBL and UCSC annotations).

Quantile normalization is implemented as described on Wikipedia_. First, we
compute an average distribution by sorting each sample (column) and taking the
mean over rows to determine the rank values. Second, we compute ranks over
columns (samples) and substitute the rank with the rank value (average
expression for each rank).

.. _Wikipedia: https://en.wikipedia.org/wiki/Quantile_normalization

Usage
=====

Install ``rnanorm`` Python package:

.. code::

    pip install rnanorm

See ``rnanorm`` command help:

.. code::

    rnanorm --help

Run ``rnanorm`` with pre-computed gene lengths:

.. code::

    rnanorm expr.tsv --cpm-output=expr.cpm.tsv --tpm-output=expr.tpm.tsv --gene-lengths=lengths.tsv

Run ``rnanorm`` with genome annotation - gene lengths will be computed on the fly:

.. code::

    rnanorm expr.tsv --cpm-output=expr.cpm.tsv --tpm-output=expr.tpm.tsv --annotation=annot.gtf

For quantile normalization we suggest using TPM expressions on the input:

.. code::

    rnanorm expr.tpm.tsv --quantile-output=expr.quantile.tsv

Contributing
============

Install ``rnanorm`` Python package for development:

.. code::

    flit install --deps=all --symlink

Run all tests and linters:

.. code::

    tox
