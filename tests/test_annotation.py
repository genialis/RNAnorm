from rnanorm.annotation import GTF


def test_union_exon_lengths(gtf_file):
    gene_lengths = GTF(gtf_file).length
    assert len(gene_lengths) == 1
    assert gene_lengths.loc["ENSG00000223972"] == 1735
