# S2D2: Sentinel-2 data deepening

[![github repo badge](https://img.shields.io/badge/github-repo-000.svg?logo=github&labelColor=gray&color=blue)](https://github.com/space-accountants/s2d2)
[![github license badge](https://img.shields.io/github/license/space-accountants/s2d2)](https://github.com/space-accountants/s2d2)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10654893.svg)](https://doi.org/10.5281/zenodo.10654893)
[![OpenSSF badge](https://bestpractices.coreinfrastructure.org/projects/8399/badge)](https://bestpractices.coreinfrastructure.org/projects/8399)
[![fair-software badge](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)
[![workflow scc badge](https://sonarcloud.io/api/project_badges/measure?project=space-accountants_s2d2&metric=coverage)](https://sonarcloud.io/dashboard?id=space-accountants_s2d2)
[![Documentation Status](https://readthedocs.org/projects/s2d2/badge/?version=latest)](https://s2d2.readthedocs.io/en/latest/?badge=latest)

## How to use s2d2

sentinel-2 data deepening


## Installation

Download and access the package folder using `git`:

```console
git clone https://github.com/space-accountants/s2d2.git
cd s2d2
```

The dependencies are most easily installed with `conda` from the `conda-forge` channel (see
[Miniforge installers](https://github.com/conda-forge/miniforge/releases) for a minimal Conda
installation). Create and activate a virtual environment with all the required dependencies:

```console
conda env create -n s2d2 -f environment.yml
conda activate s2d2
```

Install `s2d2` using `pip` (add the `-e` option to install in development mode):

```console
pip install .
```

## Documentation

Read the full project documentation [here](https://s2d2.readthedocs.io).

## Contributing

If you want to contribute to the development of s2d2,
have a look at the [contribution guidelines](CONTRIBUTING.md).

## Credits

This package was created with [Cookiecutter](https://github.com/audreyr/cookiecutter) and the [NLeSC/python-template](https://github.com/NLeSC/python-template).
