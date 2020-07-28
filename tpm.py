"""TPM normalization"""
import argparse


__version__ = "1.0.0"

parser = argparse.ArgumentParser(description='TPM normalization')
parser.add_argument('expressions_file', type=str, help='tab file with gene expression data (genes in rows, samples in cols)')
parser.add_argument('gene_lengths_file', type=str, help='tab file with gene lengths')


def main():
    args = parser.parse_args()
    print(args)

if __name__ == "__name__":
    main()
