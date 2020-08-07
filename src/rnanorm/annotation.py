"""Annotation parsing methods."""

from collections import defaultdict
from operator import itemgetter

import pandas as pd
import pybedtools as pbt

GTF_EXTENSIONS = (
    ".gtf",
    ".gtf.gz",
)


def _compute_union_exon_length(starts_ends):
    """Compute gene length based on the union of its exons.

    Before summing up exon lengths, one needs to merge all overlapping
    exons. To do this, ``pybedtools.merge()`` could be used, but simple
    hand-built merge works much faster.

    This function works only if start_ends is a list of lists. It will
    not work with list of tuples.

    """
    starts_ends.sort(key=itemgetter(0))
    merged = [starts_ends[0]]
    for current in starts_ends:
        previous = merged[-1]
        if current[0] <= previous[1]:
            previous[1] = max(previous[1], current[1])
        else:
            merged.append(current)

    return sum(end - start for (start, end) in merged)


def union_exon_lengths(annotation, gene_id_attr="gene_id"):
    """Compute gene lengths based on union exon model of genome annotation.

    Group exon start & end coordinates by gene ID (level 1) and chrom &
    strand (level 2). Then perfrom merge and length calculation for each
    chrom & strand separately. The latter is needed since
    ``gene_id_attr`` is not unique in some annotations (e.g. RefSeq).

    This function is confirmed to produce identical output as
    featureCounts output (column "Length") for the following species /
    annotations:

      - Homo sapiens:
        - UCSC hg19
        - UCSC hg38
        - ENSEMBL 92
        - ENSEMBL 100
      - Mus musculus:
        - UCSC mm10
        - ENSEMBL 92
        - ENSEMBL 100
      - Rattus norvegicus:
        - UCSC rn6
        - ENSEMBL 92
        - ENSEMBL 100
      - Macaca mulatta:
        - ENSEMBL 97
        - ENSEMBL 100

    """
    if not annotation.endswith(GTF_EXTENSIONS):
        raise ValueError(f"Input file ({annotation}) should be in GTF format")

    data = defaultdict(lambda: defaultdict(list))
    for segment in pbt.BedTool(annotation):
        if segment[2] != "exon":
            continue
        if gene_id_attr not in segment.attrs:
            raise ValueError(
                "Gene ID attribute is missing in the segment {segment[:]}. Please "
                "supply correct gene ID attribute with --gene-id-attr parameter."
            )

        data[segment.attrs[gene_id_attr]][(segment.chrom, segment.strand)].append([segment.start, segment.end])

    gene_lengths = defaultdict(int)
    for gene_id, by_chrom_strand in data.items():
        for starts_ends in by_chrom_strand.values():
            gene_lengths[gene_id] += _compute_union_exon_length(starts_ends)

    df = pd.DataFrame.from_dict(gene_lengths, orient="index", columns=["GENE_LENGTHS"])
    df.index.names = ["FEATURE_ID"]
    return df
