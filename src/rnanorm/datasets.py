"""Datasets."""
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.utils import Bunch


def load_rnaseq_toy(as_frame: bool = False) -> Bunch:
    """
    Load an artificial toy dataset representative of RNAseq data.

    :param as_frame: Return expression as pandas.DataFrame instead of Numpy
                     array.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_rnaseq_toy
    >>> dataset = load_rnaseq_toy()
    >>> dataset.exp
    array([[  200.,   300.,   500.,  2000.,  7000.],
           [  400.,   600.,  1000.,  4000., 14000.],
           [  200.,   300.,   500.,  2000., 17000.],
           [  200.,   300.,   500.,  2000.,  2000.]])
    >>> # This can be used for TPM and FPKM normalization
    >>> # since they require GTF file
    >>> dataset.gtf_path
    PosixPath('/var/folders/tmp8odklw36/example.gtf')

    """

    ds = Bunch()

    # Define expression matrix
    exp = pd.DataFrame(
        [
            [200, 300, 500, 2000, 7000],  # The ref sample
            [400, 600, 1000, 4000, 14000],  # Doubled all counts of ref
            [200, 300, 500, 2000, 17000],  # Doubled library size of ref
            [200, 300, 500, 2000, 2000],  # Halved library size of ref
        ],
        index=[f"S{i}" for i in range(1, 5)],
        columns=[f"G{i}" for i in range(1, 6)],
        dtype=np.float64,
    )
    if not as_frame:
        exp = exp.to_numpy()
    ds.exp = exp

    # GTF file
    gtf_data = [
        ["#!genome-build GRCh38.p12"],
        ["#!genome-version GRCh38"],
        ["1", ".", "gene", "1001", "1200", ".", "+", ".", 'gene_id "G1";'],
        ["1", ".", "exon", "1001", "1200", ".", "+", ".", 'gene_id "G1";'],
        ["1", ".", "gene", "2001", "3000", ".", "+", ".", 'gene_id "G2";'],
        ["1", ".", "exon", "2001", "2200", ".", "+", ".", 'gene_id "G2";'],
        ["1", ".", "exon", "2901", "3000", ".", "+", ".", 'gene_id "G2";'],
        ["1", ".", "gene", "2001", "3000", ".", "-", ".", 'gene_id "G3";'],
        ["1", ".", "exon", "2001", "2400", ".", "-", ".", 'gene_id "G3";'],
        ["1", ".", "exon", "2901", "3000", ".", "-", ".", 'gene_id "G3";'],
        ["1", ".", "gene", "2001", "3000", ".", "+", ".", 'gene_id "G4";'],
        ["1", ".", "exon", "2001", "2700", ".", "+", ".", 'gene_id "G4";'],
        ["1", ".", "exon", "2501", "3000", ".", "+", ".", 'gene_id "G4";'],
        ["2", ".", "gene", "1001", "2000", ".", "+", ".", 'gene_id "G5";'],
        ["2", ".", "exon", "1001", "2000", ".", "+", ".", 'gene_id "G5";'],
    ]
    tmp_dir = tempfile.mkdtemp(suffix=None, prefix=None, dir=None)
    gtf_path = Path(tmp_dir) / "example.gtf"
    with open(gtf_path, "wt") as handle:
        for row in gtf_data:
            handle.write("\t".join(row) + "\n")

    ds.gtf_path = gtf_path

    return ds


def load_gtex(as_frame: bool = False) -> Bunch:
    """
    Load a real RNAseq dataset from GTFx project.

    Dataset is reduced to chr21 and first 30 samples from GTEx lung V8.

    :param as_frame: Return expression as pandas.DataFrame instead of Numpy
                     array.

    .. rubric:: Examples

    >>> from rnanorm.datasets import load_gtex
    >>> dataset = load_gtex()
    >>> dataset.exp
    array([[  871,  8129,   101, ...,     0,     0,     0],
          [  852,  7076,    72, ...,    12,     0,     0],
          [  912, 11016,   174, ...,     5,     0,     0],
          ...,
          [ 1082,  6941,    50, ...,     5,     0,     0],
          [ 1248,  7052,   130, ...,     1,     0,     0],
          [ 2207, 10937,    24, ...,     0,     0,     0]])
    >>> # TPM and FPKM normalization also require GTF file
    >>> dataset.gtf_path
    PosixPath('/Users/.../gencode.v26.annotation.chr21.gtf.gz')

    """
    ds = Bunch()

    files_dir = Path(__file__).parent / "files"

    exp = pd.read_csv(files_dir / "gtex_lung.first30.chr21.csv.gz", index_col=0)
    if not as_frame:
        exp = exp.to_numpy()
    ds.exp = exp

    ds.gtf_path = files_dir / "gencode.v26.annotation.chr21.gtf.gz"

    return ds
