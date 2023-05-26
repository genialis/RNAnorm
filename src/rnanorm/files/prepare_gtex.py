# type: ignore
"""
Prepare data for GTEx dataset.

Data needs to be sub-sampled, otherwise it:
- Inflates GitHub repo size
- Inflates the size of binary and source distributions
- Increase the load times of the load_gtex function
- Increase times of tests that rely on this dataset

Reductions:
- Dataset is reduced to only genes from chr21 and first 30 samples
- GTF file is reduced to just chr21 entries. Additionally, redundant attributes
  are removed to save space

"""
import csv
from pathlib import Path

import pandas as pd

# Downloaded from here
# https://www.gencodegenes.org/human/release_26.html
# https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
GTF_INPUT_PATH = "/Users/jure/Downloads/gencode.v26.annotation.gtf.gz"

# Downloaded from here
# https://gtexportal.org/home/datasets
# https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_lung.gct.gz
EXP_INPUT_PATH = "/Users/jure/Downloads/gene_reads_2017-06-05_v8_lung.gct.gz"

GTF_PATH = Path(__file__).parent / "gencode.v26.annotation.chr21.gtf.gz"
EXP_PATH = Path(__file__).parent / "gtex_lung.first30.chr21.csv.gz"


def prepare_gtf():
    gtf = pd.read_csv(GTF_INPUT_PATH, sep="\t", comment="#", header=None)

    # Take only chr21
    gtf = gtf[gtf[0] == "chr21"]
    # Take only gene and exon lines, since only this is needed for union exon
    # length computation
    gtf = gtf[gtf[2].isin(["gene", "exon"])]
    # Remove redundant attributes

    def keep_gene_id_name(x):
        elements = x.split(";")
        elements = [e.strip() for e in elements]
        elements = [e for e in elements if e.startswith("gene_name") or e.startswith("gene_id")]
        out = "; ".join(elements) + ";"
        return out

    gtf[8] = gtf[8].apply(keep_gene_id_name)

    gtf.to_csv(
        GTF_PATH,
        sep="\t",
        header=False,
        index=False,
        quoting=csv.QUOTE_NONE,
    )

    gene_ids = gtf[8].str.extract(r'gene_id "([^"]+)";', expand=False).unique()
    return gene_ids


def prepare_exp(gene_ids):
    exp = pd.read_csv(EXP_INPUT_PATH, sep="\t", skiprows=2, index_col=1)
    exp.drop(columns=["Description", "id"], inplace=True)
    exp.index.name = "ENSEMBL_ID"
    # Sort to make it reproducible and cast to list
    common_gene_ids = sorted(set(gene_ids).intersection(exp.index))
    exp = exp.loc[common_gene_ids]
    exp = exp.iloc[:, :30]

    # Save the transposed version
    exp.T.to_csv(EXP_PATH)


if __name__ == "__main__":
    # First prepare GTF, since we need gene ID's from chr21 to extract data
    # from expression matrix
    gene_ids = prepare_gtf()

    prepare_exp(gene_ids)
