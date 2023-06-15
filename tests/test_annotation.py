import pandas as pd

from rnanorm.annotation import GTF


def test_union_exon_lengths(gtf_file, exp):
    gene_lengths = GTF(gtf_file).length
    assert len(gene_lengths) == 5
    pd.testing.assert_series_equal(
        gene_lengths,
        pd.Series(
            [200, 300, 500, 1000, 1000],
            index=pd.Series(["Gene_1", "Gene_2", "Gene_3", "Gene_4", "Gene_5"], name="gene_id"),
            name="gene_length",
        ),
    )
