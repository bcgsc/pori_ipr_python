
# IPR

![build](https://github.com/bcgsc/pori_ipr_python/workflows/build/badge.svg) [![PyPi](https://img.shields.io/pypi/v/ipr.svg)](https://pypi.org/project/ipr) [![codecov](https://codecov.io/gh/bcgsc/pori_ipr_python/branch/master/graph/badge.svg)](https://codecov.io/gh/bcgsc/pori_ipr_python) [![PyPI - Downloads](https://img.shields.io/pypi/dm/ipr)](https://pypistats.org/packages/ipr)

This python tool takes in variant inputs as tab-delimited files and annotates them using GraphKB.
The resulting output is uploaded to IPR as a report. Additional report content such as images and
metadata can be passed to be included in the report upload.


## Getting Started

### Install (For developers)

clone this repository

```
git clone https://github.com/bcgsc/pori_ipr_python.git
cd pori_ipr_python
```

create a virtual environment

```
python3 -m venv venv
source venv/bin/activate
```

install the package and its development dependencies

```
pip install -e .[dev]
```

Run the tests

```
pytest tests
```

## Generating the Documentation

This documentation is generated using [mkdocs](https://www.mkdocs.org), [mkdocs-material](https://pypi.org/project/mkdocs-material), and [markdown_refdocs](https://pypi.org/project/markdown-refdocs).

First install the documentation dependencies

```bash
pip install .[doc]
```

Then generate the user manual files

```bash
markdown_refdocs ipr -o docs/reference
mkdocs build
```

There should now be static html files under `build-docs`. To view the files, serve the folder using
the built-in python http server

```bash
python3 -m http.server -d build-docs
```
