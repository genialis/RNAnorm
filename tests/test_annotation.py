import pandas as pd

from rnanorm.annotation import GTF


def test_union_exon_lengths(gtf_file, exp):
    gene_lengths = GTF(gtf_file).length
    assert len(gene_lengths) == 4
    pd.testing.assert_series_equal(
        gene_lengths,
        pd.Series(
            [200, 300, 500, 1000],
            index=pd.Series(["G1", "G2", "G3", "G4"], name="gene_id"),
            name="gene_length",
        ),
    )
