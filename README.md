# OmicsExplorer

OmicsExplorer is an interactive Shiny-based application for the visualization and exploration of genomic and transcriptomic analysis results.

The platform was developed in the context of rare disease diagnostics and omics data interpretation, providing an intuitive graphical interface for non-computational users.

Supported data sources include:

- RNA sequencing (RNA-seq)
- Whole exome sequencing (WES)
- Whole genome sequencing (WGS)
- DROP analyses:
  - Aberrant Expression (AE)
  - Aberrant Splicing (AS)
  - Mono-Allelic Expression (MAE)
- Coverage tracks (bigWig)

---

# Key Features

## Genomic exploration

- Dynamic variant tables
- Interactive filtering and search
- WES/WGS result exploration
- Coverage visualization from bigWig files

## Transcriptomic exploration

- Gene expression visualization
- RNA-seq sample comparison plots
- Transcript-level exploration
- Expression profile inspection

## Functional evidence integration

- Integration of DROP analyses
- Aberrant expression (AE)
- Aberrant splicing (AS)
- Mono-allelic expression (MAE)

## Quality control and reporting

- Interactive quality control summaries
- QC plots and metrics visualization
- Multi-omics result overview

---

# Quick Start

## Clone repository

```bash
git clone https://github.com/sandrasnchz/OmicsExplorer.git
cd OmicsExplorer
```

## Restore dependencies

```r
install.packages("renv")
renv::restore()
```

## Run application

```r
library(shiny)

runApp("app")
```

---

# Portable Release

Preconfigured portable versions are available in the GitHub Releases section.

Portable releases include:

✔ Portable R environment  
✔ Preinstalled package dependencies  
✔ Application launcher  
✔ No additional installation required  

### Quick start

1. Download the latest portable release from **Releases**
2. Extract the `.zip` file
3. Open:

```text
launcher/Abrir_OmicsExplorer.bat
```

4. The application will automatically start and open in your default web browser.

---

# Technologies

OmicsExplorer is mainly built with R.

The application is designed to integrate with bioinformatics workflows such as:

- nf-core/sarek
- nf-core/rnaseq
- DROP

Dependency management and reproducibility are supported through:

- `renv`
- Git version control
- Modular application architecture

---

# Repository Structure

```text
OmicsExplorer/
│
├── app/
│   ├── app.R
│   ├── modules/       # Shiny modules
│   ├── queries/       # Database queries
│   └── www/           # Static resources
│
├── data/              # Processed datasets
├── db/                # Database file
│
├── renv/              # Dependency management
├── renv.lock
│
├── README.md
├── LICENSE
└── .gitignore
```

---

# Data Availability

This repository intentionally excludes:

- BAM / CRAM files
- VCF files
- bigWig files
- Patient-derived datasets
- Large genomic databases

These files are excluded due to privacy requirements and repository size limitations.

---

# Citation

If you use OmicsExplorer in academic work, please cite the repository and/or the associated master thesis project.

---

# License

This project is distributed under the MIT License.

See the `LICENSE` file for details.

---

# Author

**Sandra Sánchez Pinilla**

Bioinformatics Unit (BU-ISCIII)  
Instituto de Salud Carlos III (ISCIII)
