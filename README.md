# NGS Trim - a generic module for read trimming

[![Snakemake](https://img.shields.io/badge/snakemake-≥9.13.5-brightgreen.svg)](https://snakemake.readthedocs.io/en/stable/)
![CI](https://github.com/GeoGenetics/ngs-trim/actions/workflows/ci.yml/badge.svg)

This module performs NGS trimming steps:
- Read trimming:
  - [AdapterRemoval2](https://adapterremoval.readthedocs.io/en/latest/) / [BBduk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/) / [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) / [Fastp](https://github.com/OpenGene/fastp) / [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

And QC steps:
- [MultiQC](https://multiqc.info/) (aggregates QC from several of the tools above, plus):
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (on both raw and trimmed reads)

This module can be used directly, but is designed to be used together with other downstream modules.

## Authors

* Filipe G. Vieira

## Usage

For an example on how to use this module, check repo [aeDNA](https://github.com/GeoGenetics/aeDNA).
