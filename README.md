# OmicsExplorer

OmicsExplorer is an interactive Shiny-based application for the visualization and exploration of genomic and transcriptomic analysis results.

The platform has been developed in the context of rare disease diagnostics and bioinformatics workflows, integrating results from:

* RNA sequencing (RNA-seq)
* Whole exome sequencing (WES)
* Whole genome sequencing (WGS)
* Aberrant expression and splicing analyses
* Variant prioritization workflows
* Coverage and expression visualization

The application is designed to facilitate the interpretation of omics data by non-computational users through an intuitive graphical interface.

---

# Features

Current modules include:

* Interactive gene expression visualization
* RNA-seq expression comparison plots
* Dynamic variant tables
* Coverage visualization using bigWig files
* Integration of DROP results
* Interactive filtering and search tools
* Transcriptomic and genomic result exploration
* QC and analysis result visualization

The application is modular and can be extended with additional analysis and visualization components.

---

# Technologies

OmicsExplorer is mainly built with:

* R
* Shiny
* Bioconductor
* DuckDB
* Plotly
* DT
* dplyr

The project is intended to integrate smoothly with bioinformatics workflows such as:

* nf-core/rnaseq
* nf-core/sarek
* DROP
* Nextflow-based pipelines

---

# Repository Structure

```text
OmicsExplorer/
│
├── app/
│   ├── app.R
│   ├── modules/
│   ├── www/
│   └── data/
│
├── scripts/
│   ├── run_app.R
│   └── install_packages.R
│
├── launcher/
│   └── Abrir_OmicsExplorer.bat
│
├── README.md
├── LICENSE
└── renv.lock
```

---

# Installation

## 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/OmicsExplorer.git
cd OmicsExplorer
```

## 2. Install R dependencies

Install required packages from R:

```r
install.packages("renv")
renv::restore()
```

Alternatively, dependencies can be installed manually.

---

# Running the Application

## Option 1 — Run from R

```r
library(shiny)
runApp("app")
```

## Option 2 — Windows launcher

Double click:

```text
Abrir_OmicsExplorer.bat
```

This automatically starts the Shiny application and opens it in the default web browser.

---

# Portable Version

A portable version of OmicsExplorer can be distributed internally including:

* Portable R installation
* Required R packages
* Application launcher
* Shiny application

This allows non-technical users to run the application without installing R or additional dependencies.

---

# Data

This repository does not include:

* BAM/CRAM files
* VCF files
* bigWig files
* Large databases
* Patient-sensitive data

These files should be stored separately due to size and privacy constraints.

---

# Reproducibility

The project uses:

* Git version control
* Modular application design
* `renv` dependency management
* Standardized bioinformatics workflows

to improve reproducibility and maintainability.

---

# License

This project is licensed under the MIT License.

You are free to:

* Use
* Modify
* Distribute
* Reuse

the software, provided that the original copyright and license notice are included.

See the LICENSE file for details.

---

# Recommended LICENSE file

Create a file named:

```text
LICENSE
```

with the following content:

```text
MIT License

Copyright (c) 2026 Sandra Sánchez Pinilla

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

# Citation

If you use this project in academic work, please cite the repository or the associated thesis/project.

---

# Author

Sandra Sánchez Pinilla
Biomedical Bioinformatics and Omics Analysis

---

# Acknowledgements

This project has been developed using concepts and workflows inspired by:

* DROP
* nf-core
* Bioconductor
* Shiny
* Nextflow

and related open-source bioinformatics tools.
