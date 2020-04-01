
# Genomic Report

This python tool takes in variant inputs as tab-delimited files and annotates them using GraphKB.
The resulting output is uploaded to IPR as a report. Additional report content such as images and
metadata can be passed to be included in the report upload.


## Getting Started

### Install (For developers)

clone this repository

```
git clone ssh://git@svn.bcgsc.ca:7999/sdev/genomic_report.git
cd genomic_report
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

## Deployment (Publishing)

Install the deployment dependencies

```bash
pip install .[deploy]
```

Build the distribution files

```bash
python setup.py install sdist bdist_wheel
```

Upload the distibutions to the package server (-r defined in your pypirc)

```bash
twine upload -r bcgsc dist/*
```

## Deployment (Scripts)

A buildout config is included by default which will build all console scripts defined
in the package.

create a virtual environment and install buildout

```bash
python3 -m venv venv
source venv/bin/activate
pip install -U pip setuptools zc.buildout
```

run buildout

```bash
buildout
```

This will create a directory `bin` with the executable scripts
