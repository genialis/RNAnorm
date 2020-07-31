====================
RNASeq normalization
====================

Normalization of RNA-seq gene expression. Supported methods: Transcript per kilobase million (TPM) and Counts per million (CPM).

Usage
=====

Install ``rnanorm`` Python package:

.. code::

    pip install rnanorm

Run ``rnanorm`` command line tool:

.. code::

    rnanorm expr.tsv --cpm-output=expr.cpm.tsv --tpm-output=expr.tpm.tsv --gene-lengths=lengths.tsv

Contributing
============

Install ``rnanorm`` Python package for development:

.. code::

    flit install --deps=all --symlink

Run all tests and linters:

.. code::

    tox

