"""GTF file parsing."""

import io
from functools import lru_cache
from pathlib import Path
from typing import Union

import pandas as pd


class GTF:
    """GTF file parser."""

    def __init__(self, gtf: Union[str, Path]) -> None:
        """Initialize class.

        :param gtf: GTF file
        """
        self.gtf = gtf

    @property
    def length(self) -> pd.Series:
        """Gene length by union exon length model."""
        return self._get_df()["gene_length"]

    @lru_cache()
    def _get_df(self) -> pd.DataFrame:
        """Load data from disc or fetch it from server and cache it on disc."""
        gtf_df = self._parse_gtf(self.gtf)

        # Compute gene length
        gene_length_df = self._gene_length(gtf_df)

        # Use only rows that represent genes
        gtf_df = gtf_df[gtf_df["feature_type"] == "gene"]
        gtf_df["gene_length"] = gene_length_df

        return gtf_df

    @staticmethod
    def _parse_gtf(file: Union[io.BytesIO, str, Path]) -> pd.DataFrame:
        """Parse GTF file."""

        # Read columns we need: 0=chromosome, 2=feature type, 3=genomic start
        # location, 4=genomic stop location, 6=genomic strand, 8=attributes
        df = pd.read_csv(
            file,
            sep="\t",
            comment="#",
            header=None,
            usecols=[0, 2, 3, 4, 6, 8],
            low_memory=False,
        )

        return pd.DataFrame(
            {
                "feature_type": df[2],
                "gene_id": df[8].str.extract(r'gene_id "(.+?)"', expand=False).values,
                "gene_name": df[8].str.extract(r'gene_name "(.+?)"', expand=False).values,
                "chromosome": df[0],
                "strand": df[6],
                "start": df[3],
                "end": df[4] + 1,  # Note: GTF is 1-based and includes stop position
            },
        ).set_index("gene_id")

    def _gene_length(self, gtf_df: pd.DataFrame, gene_id_attr: str = "gene_id") -> pd.Series:
        """Compute gene lengths based on union exon model.

        Group exon start & end coordinates by gene ID & chromosome &
        strand. Then perfrom merge and length calculation for each
        group separately. The latter is needed since ``gene_id_attr``
        is not unique in some annotations (e.g., RefSeq).
        """
        gtf_df = gtf_df[gtf_df["feature_type"] == "exon"]

        group_exons = gtf_df.groupby([gene_id_attr, "chromosome", "strand"])[["start", "end"]]
        group_length = group_exons.apply(self._compute_union_exon_length)

        gene_length = group_length.groupby(gene_id_attr).sum()

        return gene_length

    def _compute_union_exon_length(self, starts_ends: pd.DataFrame) -> pd.Series:
        """Compute gene length based on the union of its exons.

        Before summing up exon lengths, one needs to merge all
        overlapping exons.
        """
        exons = starts_ends.values
        exons = exons[exons[:, 0].argsort()]
        length = 0  # length of a rolling union of exons
        end = -1  # end of the rolling union of exons
        for exon in exons:
            if exon[1] > end:
                if exon[0] > end:
                    length += exon[1] - exon[0]
                else:
                    length += exon[1] - end
                end = exon[1]
        return length
