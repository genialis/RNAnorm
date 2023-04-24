import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def between_sample_data():
    """Toy, example data for between sample normalization."""
    n_samples = 3
    n_genes = 6
    genes = [f"G{i + 1}" for i in range(n_genes)]
    samples = [f"S{i + 1}" for i in range(n_samples)]

    data = pd.DataFrame(
        [
            # S1, lib_size = 11.111
            [0, 1, 10, 100, 1000, 10000],
            # Same as S1, but double all the counts, lib_size = 22.222
            # This needs normalization by library size (CPM)
            [0, 2, 20, 200, 2000, 20000],
            # Same as S1, only one gene is super higher so that lib_size = S2
            # This one needs correction of effective library size
            [0, 1, 10, 100, 1000, 21111],
        ],
        index=samples,
        columns=genes,
        dtype=np.float64,
    )

    expected_norm_factors = pd.Series([1.2599, 1.2599, 0.62996], index=samples)

    expected_data1 = pd.DataFrame(
        [
            [0, 71.4337, 714.337, 7143.37, 71433.7, 714337],
            [0, 71.4337, 714.337, 7143.37, 71433.7, 714337],
            [0, 71.4337, 714.337, 7143.37, 71433.7, 1508038],
        ],
        index=samples,
        columns=genes,
        dtype=np.float64,
    )

    expected_data2 = pd.DataFrame(
        [
            [0.0, 0.793701, 7.937005, 79.370053, 793.700, 7937],
            [0.0, 1.5874, 15.874, 158.74, 1587.4, 15874],
            [0.0, 1.5874, 15.874, 158.74, 1587.4, 33511],
        ],
        index=samples,
        columns=genes,
        dtype=np.float64,
    )
    return (
        data,
        expected_norm_factors,
        expected_data1,
        expected_data2,
    )


@pytest.fixture
def gtf_file(tmp_path):
    # This is first gene in ENSEMBL 92
    data = [
        ["#!genome-build GRCh38.p12"],
        ["#!genome-version GRCh38"],
        ["#!genome-date 2013-12"],
        ["#!genome-build-accession NCBI:GCA_000001405.27"],
        ["#!genebuild-last-updated 2018-01"],
        ["1", "havana", "gene", "11869", "14409", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        [
            "1",
            "havana",
            "transcript",
            "11869",
            "14409",
            ".",
            "+",
            ".",
            'gene_id "ENSG00000223972";',
        ],
        ["1", "havana", "exon", "11869", "12227", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12613", "12721", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "13221", "14409", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        [
            "1",
            "havana",
            "transcript",
            "12010",
            "13670",
            ".",
            "+",
            ".",
            'gene_id "ENSG00000223972";',
        ],
        ["1", "havana", "exon", "12010", "12057", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12179", "12227", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12613", "12697", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12975", "13052", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "13221", "13374", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "13453", "13670", ".", "+", ".", 'gene_id "ENSG00000223972";'],
    ]

    gtf_path = tmp_path / "annotation.gtf"
    with open(gtf_path, "wt") as handle:
        for row in data:
            handle.write("\t".join(row) + "\n")

    return gtf_path
