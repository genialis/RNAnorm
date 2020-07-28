# Install

Install `tpm` Python package for development:

```
flit install --symlink
```

# Contributing

We suggest using linters:

```
black --check oncologie
flake8 oncologie
isort --recursive --check-only --diff oncologie
```
