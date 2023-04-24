====================
RNAseq normalization
====================

Python implementation of common RNAseq normalization methods:

- CPM (Counts per million)
- FPKM (Fragments per kilobase million)
- TPM (Transcripts per million)
- UQ (Upper quartile)
- CUF (Counts adjusted with UQ factors)
- TMM (Trimmed mean of M-values)
- CTF (Counts adjusted with TMM factors)

**WARNING**: The content of this repository is under substantial restructuring. Some parts do not work yet and some features are not documented.


Installation
============

This repository is getting as new structure and API. However, it is not
yet mature enough to be published on Pypi. This means that the usual
``pip install rnanorm`` will install the old version of ``rnanorm``.

To use the newest code please install from source::

   git clone git@github.com:genialis/RNAnorm.git
   cd RNAnorm
   pip install .


Normalize from Python
=====================

Example::

   import pandas as pd
   from rnanorm import CPM

   exp_raw = pd.DataFrame({"G1": [100000, 300000], "G2": [400000, 700000]}, index=["S1", "S2"])
   CPM().fit_transform(exp_raw)


Normalize from command line
===========================

Coming soon.

.. TODO

Contributing
============

To learn about contributing to the code base, read the  Contributing_ section
of our documentation.

.. _Contributing: docs/contributing.rst
