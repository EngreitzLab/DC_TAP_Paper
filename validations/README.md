# DC-TAP-seq Validations

> [!WARNING]
> This directory is a work in progress and may not be complete or fully accurate.

This folder contains analysis code and notebook for DC-TAP-seq validation experiments.

## Overview

The analyses are organized in subdirectories named by date (`YYYY-MM-DD`) under the
`analyses/` directory. `README.md` file inside this folder contains a description
of the analysis performed on each date, as well as the code and results.

Results are stored in the `results/` directory, organized by date corresponding to
the analysis that generated them.

## Requirements

### Software Dependencies

| Software | Version  |
| -------- | -------- |
| python   | ≥ 3.12   |
| bedtools | ≥ 2.31.1 |
| blat     | ≥ 35     |

### Environment Setup

It is recommended to use a Python virtual environment.  
To create and activate an environment and install the required libraries:

```bash
cd validations
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Running notebooks

```bash
marimo edit analyses/YYYY-MM-DD/my_notebook.py
```
