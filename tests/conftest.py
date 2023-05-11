import numpy as np
import pandas as pd
import pytest

N_SAMPLES = 4
N_GENES = 5


@pytest.fixture
def exp():
    return pd.DataFrame(
        [
            [200, 300, 500, 2000, 7000],  # The ref sample
            [400, 600, 1000, 4000, 14000],  # Doubled all counts of ref
            [200, 300, 500, 2000, 17000],  # Doubled library size of ref
            [200, 300, 500, 2000, 2000],  # Halved library size of ref
        ],
        index=[f"S{i}" for i in range(1, N_SAMPLES + 1)],
        columns=[f"G{i}" for i in range(1, N_GENES + 1)],
        dtype=np.float64,
    )


@pytest.fixture
def gtf_file(tmp_path):
    data = [
        ["#!genome-build GRCh38.p12"],
        ["#!genome-version GRCh38"],
        # Gene 1 - just one exon
        ["1", ".", "gene", "1001", "1200", ".", "+", ".", 'gene_id "G1";'],
        ["1", ".", "exon", "1001", "1200", ".", "+", ".", 'gene_id "G1";'],
        # Gene 2 - two exons
        ["1", ".", "gene", "2001", "3000", ".", "+", ".", 'gene_id "G2";'],
        ["1", ".", "exon", "2001", "2200", ".", "+", ".", 'gene_id "G2";'],
        ["1", ".", "exon", "2901", "3000", ".", "+", ".", 'gene_id "G2";'],
        # Gene 3 - two exons on the opposite strand of Gene 2
        ["1", ".", "gene", "2001", "3000", ".", "-", ".", 'gene_id "G3";'],
        ["1", ".", "exon", "2001", "2400", ".", "-", ".", 'gene_id "G3";'],
        ["1", ".", "exon", "2901", "3000", ".", "-", ".", 'gene_id "G3";'],
        # Gene 4 - two overlapping exons
        ["1", ".", "gene", "2001", "3000", ".", "+", ".", 'gene_id "G4";'],
        ["1", ".", "exon", "2001", "2700", ".", "+", ".", 'gene_id "G4";'],
        ["1", ".", "exon", "2501", "3000", ".", "+", ".", 'gene_id "G4";'],
        # Gene 5 is deliberately missing
    ]

    gtf_path = tmp_path / "annotation.gtf"
    with open(gtf_path, "wt") as handle:
        for row in data:
            handle.write("\t".join(row) + "\n")

    return gtf_path


@pytest.fixture
def expected_factors():
    return [1.0, 1.0, 0.5, 2.0]
