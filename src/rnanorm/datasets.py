"""Datasets."""

from pathlib import Path

import pandas as pd
from sklearn.utils import Bunch

FILES_DIR = Path(__file__).parent / "files"


def load_toy_data() -> Bunch:
    """
    Load a minimal, made-up raw count RNA-seq dataset.

    This dataset is not representative of any real RNA-seq dataset. However,
    it is small and has just enough complexity to showcase the effects of
    methods. It is also simple enough to make calculation by hand possible.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_toy_data
    >>> dataset = load_toy_data()
    >>> dataset.exp
              Gene_1  Gene_2  Gene_3  Gene_4  Gene_5
    Sample_1     200     300     500    2000    7000
    Sample_2     400     600    1000    4000   14000
    Sample_3     200     300     500    2000   17000
    Sample_4     200     300     500    2000    2000
    >>> # TPM and FPKM normalization also require GTF file
    >>> dataset.gtf_path
    PosixPath('/Users/me/.../toy.gtf')

    """
    ds = Bunch()
    ds.exp = pd.read_csv(FILES_DIR / "toy_exp.csv", index_col=0)
    ds.gtf_path = FILES_DIR / "toy.gtf"

    return ds


def load_gtex() -> Bunch:
    """
    Load a real RNA-seq dataset from GTFx project.

    Dataset is sampled to contain just chr21 and first 30 samples from GTEx
    lung V8.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_gtex
    >>> dataset = load_gtex()
    >>> dataset.exp
                              ENSG00000141956.13  ENSG00000141959.16   ...
    GTEX-111CU-0326-SM-5GZXO                 871                8129   ...
    GTEX-111FC-1126-SM-5GZWU                 852                7076   ...
    GTEX-111VG-0726-SM-5GIDC                 912               11016   ...
    ...                                      ...                 ...   ...
    >>> # TPM and FPKM normalization also require GTF file
    >>> dataset.gtf_path
    PosixPath('/Users/.../gencode.v26.annotation.chr21.gtf.gz')

    """
    ds = Bunch()
    ds.exp = pd.read_csv(FILES_DIR / "gtex_lung.first30.chr21.csv.gz", index_col=0)
    ds.gtf_path = FILES_DIR / "gencode.v26.annotation.chr21.gtf.gz"

    return ds
