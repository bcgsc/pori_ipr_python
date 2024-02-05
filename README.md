
# IPR

![build](https://github.com/bcgsc/pori_ipr_python/workflows/build/badge.svg) [![PyPi](https://img.shields.io/pypi/v/ipr.svg)](https://pypi.org/project/ipr) [![codecov](https://codecov.io/gh/bcgsc/pori_ipr_python/branch/master/graph/badge.svg)](https://codecov.io/gh/bcgsc/pori_ipr_python) [![PyPI - Downloads](https://img.shields.io/pypi/dm/ipr)](https://pypistats.org/packages/ipr) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5730671.svg)](https://doi.org/10.5281/zenodo.5730671)

This repository is part of the [Platform for Oncogenomic Reporting and Interpretation (PORI)](https://bcgsc.github.io/pori/).

This python tool takes in variant inputs as tab-delimited files and annotates them using GraphKB.
The resulting output is uploaded to IPR as a report. Additional report content such as images and
metadata can be passed to be included in the report upload.

For documentation on how to create reports using this adaptor, see the [main documentation site](https://bcgsc.github.io/pori/) for the platform.

## Getting Started

### Install (For developers)

clone this repository

```bash
git clone https://github.com/bcgsc/pori_ipr_python.git
cd pori_ipr_python
```

create a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate
```

install the package and its development dependencies

```bash
pip install -U pip setuptools
pip install -e .[dev]
```

Run the tests

```bash
pytest tests
```

## Documentation

The user documentation for this tool is hosted with the [main documentation site](https://bcgsc.github.io/pori/).

Developers: Any updates to this tool should be edited and reflected in the main site documentation as well.


# GraphKB (Python)

![build](https://github.com/bcgsc/pori_graphkb_python/workflows/build/badge.svg) [![PyPi](https://img.shields.io/pypi/v/graphkb.svg)](https://pypi.org/project/graphkb) [![codecov](https://codecov.io/gh/bcgsc/pori_graphkb_python/branch/master/graph/badge.svg)](https://codecov.io/gh/bcgsc/pori_graphkb_python) [![PyPI - Downloads](https://img.shields.io/pypi/dm/graphkb)](https://pypistats.org/packages/graphkb) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5730523.svg)](https://doi.org/10.5281/zenodo.5730523)

This repository is part of the [platform for oncogenomic reporting and interpretation](https://github.com/bcgsc/pori).

Python adapter package for querying the GraphKB API. See the [user manual](https://bcgsc.github.io/pori/graphkb/scripting/)

- [Getting Started](#getting-started)
  - [Install (For developers)](#install-for-developers)
  - [Run Tests](#run-tests)
- [Generating the Documentation](#generating-the-documentation)
- [Deployment (Publishing)](#deployment-publishing)

## Getting Started

### Install (For developers)

clone this repository

```bash
git clone https://github.com/bcgsc/pori_graphkb_python
cd pori_graphkb_python
```

create a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate
```

install the package and its development dependencies

```bash
pip install -U pip setuptools
pip install -e .[dev]
```

### Run Tests

```bash
pytest tests
```

## Generating the Documentation

User documentation for this repository is hosted in the [central PORI repository](https://github.com/bcgsc/pori/)

## Deployment (Publishing)

Install the deployment dependencies

```bash
pip install .[deploy]
```

Build the distribution files

```bash
python setup.py install sdist bdist_wheel
```

Upload the distibutions to the package server (`-r` is defined in your pypirc)

```bash
twine upload -r bcgsc dist/*
```
