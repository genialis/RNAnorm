import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


class Dataset:
    exp = pd.DataFrame()
    gtf_path = Path()


def load_example() -> Dataset:
    """Load minimal example dataset with expressions and a GTF file."""
    ds = Dataset()

    # Define expressions matrix
    ds.exp = pd.DataFrame(
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
