Normalization of RNA-seq gene expression. Supported methods: Transcript per kilobase million (TPM) and Counts per million (CPM).

# Usage

Install `rnanorm` Python package:

```
flit install
```

Run `rnanorm` command line tool:

```
rnanorm expr.tsv --cpm-output=expr.cpm.tsv --tpm-output=expr.tpm.tsv --gene-lengths=lengths.tsv
```

# Contributing

Install `rnanorm` Python package for development:

```
flit install --symlink
```

Run all tests:

```
tox
```

Run linters:

```
tox -e linters
```
