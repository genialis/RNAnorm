##########
Change Log
##########

All notable changes to this project are documented in this file.


==================
2.2.0 - 2025-01-29
==================

Added
-----
- Add support for Python 3.13

Changed
-------
- Remove support for Python 3.8
- Update deprecated ``sklearn`` validation methods


==================
2.1.0 - 2023-11-03
==================

Added
-----
- Add support for Python 3.12


==================
2.0.1 - 2023-09-27
==================

Added
-----
- Provide an example of gene lengths file in the ``README.rst``

Fixed
-----
- Fix ``fit`` methods to match sklearn signature
- Fix ``LibrarySize`` to always return ``np.array`` in private functions
- Fix ``get_norm_factors`` to follow the ``set_output`` config
- Remove leftovers of ``as_frame`` parameter


==================
2.0.0 - 2023-06-21
==================

Added
-----
- Implementation of the following methods:

    - CPM
    - FPKM
    - TPM
    - UQ
    - CUF
    - TMM
    - CTF

- Add a "toy" and GTEx dataset
- Add command line interface for all of the above methods
- Add tests
- Support calculation of gene lengths from GTf or gene lengths file in TPM /
  FPKM
