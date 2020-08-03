import pytest

from rnanorm.annotation import union_exon_lengths


@pytest.fixture(scope="session")
def create_gtf_file(tmpdir_factory):
    # This is first gene in ENSEMBL 92
    data = [
        ["#!genome-build GRCh38.p12"],
        ["#!genome-version GRCh38"],
        ["#!genome-date 2013-12"],
        ["#!genome-build-accession NCBI:GCA_000001405.27"],
        ["#!genebuild-last-updated 2018-01"],
        ["1", "havana", "gene", "11869", "14409", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "transcript", "11869", "14409", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "11869", "12227", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12613", "12721", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "13221", "14409", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "transcript", "12010", "13670", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12010", "12057", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12179", "12227", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12613", "12697", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "12975", "13052", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "13221", "13374", ".", "+", ".", 'gene_id "ENSG00000223972";'],
        ["1", "havana", "exon", "13453", "13670", ".", "+", ".", 'gene_id "ENSG00000223972";'],
    ]

    output_dir = tmpdir_factory.mktemp("output")
    gtf_path = output_dir.join("annotation.gtf")
    with open(gtf_path, "wt") as handle:
        for row in data:
            handle.write("\t".join(row) + "\n")

    return gtf_path


def test_union_exon_lengths(create_gtf_file):
    gene_lengths = union_exon_lengths(str(create_gtf_file))
    assert len(gene_lengths) == 1
    assert gene_lengths.loc["ENSG00000223972", "GENE_LENGTHS"] == 1735
