# Genomic Report User Manual

The Genomic Report python package is a tool for generating an IPR report with GraphKB annotations.
Lists of variants are parsed, annotated, and then uploaded as part of a report.

## Getting Started

Install the package with pip

```bash
pip install genomic_report
```

This will require python 3.6 or greater.

this can now be used as a command line tool

```bash
genomic_report -h
```

or as part of a script

```python
from argparse import Namespace

from genomic_report.main import create_report

create_report(...)
```


