
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
