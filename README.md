Normalization of RNA-seq gen expression. Supported methods: Transcript per kilobase million (TPM) and Counts per million (CPM).

# Usage

Install `rnanorm` Python package for development:

```
flit install
```

Run `rnanorm` command line tool:

```
rnanorm gene_expressions.tsv gene_lengths.tsv
```

# Contributing

Install `rnanorm` Python package for development:

```
flit install --symlink
```

We suggest using linters:

```
black --check rnanorm
flake8 rnanorm
isort --recursive --check-only --diff rnanorm
```
