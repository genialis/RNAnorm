# Usage

Install `tpm` Python package for development:

```
flit install
```

Run `tpm` command line tool:

```
tpm gene_expressions.tsv gene_lengths.tsv
```

# Contributing

Install `tpm` Python package for development:

```
flit install --symlink
```

We suggest using linters:

```
black --check oncologie
flake8 oncologie
isort --recursive --check-only --diff oncologie
```
