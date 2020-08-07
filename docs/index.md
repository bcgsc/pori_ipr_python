# IPR Adapater User Manual

The IPR python package is a tool for generating an IPR report with GraphKB annotations.
Lists of variants are parsed, annotated, and then uploaded as part of a report. See the
user manual [here](https://bcgsc.github.io/pori_ipr_python)

## Getting Started

Install the package with pip

```bash
pip install ipr
```

This will require python 3.6 or greater.

this can now be used as a command line tool

```bash
ipr -h
```

or as part of a script

```python
from argparse import Namespace

from ipr.main import create_report

create_report(...)
```
