# OmicsExplorer

OmicsExplorer is an interactive Shiny-based application for the visualization and exploration of genomic and transcriptomic analysis results.

The platform was developed in the context of rare disease diagnostics and omics data interpretation, providing an intuitive graphical interface for non-computational users.

It integrates and visualizes information from:

- RNA sequencing (RNA-seq)
- Whole exome sequencing (WES)
- Whole genome sequencing (WGS)
- Aberrant expression analyses
- Aberrant splicing analyses
- Variant prioritization workflows
- Coverage and expression data

---

# Features

Current functionalities include:

- Interactive gene expression visualization
- RNA-seq expression comparison plots
- Dynamic variant tables
- Coverage visualization from bigWig files
- Integration of DROP analysis results
- Interactive filtering and search tools
- Transcriptomic and genomic result exploration
- Quality control and summary visualizations

The application follows a modular architecture and can be extended with additional analysis and visualization modules.

---

# Technologies

OmicsExplorer is mainly built with:

- R
- Shiny
- Bioconductor
- DuckDB
- Plotly
- DT
- dplyr

The application is designed to integrate with bioinformatics workflows such as:

- nf-core/rnaseq
- nf-core/sarek
- DROP
- Nextflow-based pipelines

---

# Repository Structure

```text
OmicsExplorer/
│
├── app/
│   ├── app.R
│   ├── modules/
│   ├── queries/
│   └── www/
│
├── data/
├── db/
│
├── renv/
├── renv.lock
│
├── README.md
├── LICENSE
└── .gitignore
```

---

# Installation

## 1. Clone the repository

```bash
git clone https://github.com/sandrasnchz/OmicsExplorer.git
cd OmicsExplorer
```

## 2. Restore dependencies

```r
install.packages("renv")
renv::restore()
```

---

# Running the Application

Launch the application from R:

```r
library(shiny)
runApp("app")
```

---

# Releases

Preconfigured portable versions are available in the GitHub Releases section.

The portable release includes:

- Portable R environment
- Preinstalled package dependencies
- Application launcher
- No additional installation required

This allows non-technical users to run the application with a simple double click.

---

# Data

This repository does not include:

- BAM/CRAM files
- VCF files
- bigWig files
- Patient-sensitive data
- Large genomic databases

These files are excluded due to size and privacy constraints.

---

# Reproducibility

OmicsExplorer uses:

- Git version control
- Modular application design
- `renv` dependency management
- Standardized bioinformatics workflows

to improve reproducibility and maintainability.

---

# License

This project is distributed under the MIT License.

See the `LICENSE` file for details.

---

# Citation

If you use OmicsExplorer in academic work, please cite the repository and/or the associated thesis project.

---

# Author

**Sandra Sánchez Pinilla**  

---

# Acknowledgements

This project has been developed using concepts and workflows inspired by:

- DROP
- nf-core
- Bioconductor
- Shiny
- Nextflow

and related open-source bioinformatics tools.
