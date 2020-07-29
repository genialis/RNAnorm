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

Run all tests:

```
tox -e
```

Run linters:

```
tox -e linters
```
