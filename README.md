# AstrocyteHeterogeneity
Contains snRNAseq preprocessing and analysis scripts, CellProfiler pipelines, and Fiji and MATLAB scripts for analysis of ExR and multiExR images.

Accompanies the manuscript by Schroeder et al., 2025 ["A transcriptomic atlas of astrocyte heterogeneity across space and time in mouse and marmoset"](https://www.biorxiv.org/content/10.1101/2024.10.11.617802v3).

This respository is organized topically based on our manuscript. There is a different subfolder for each type of analysis. Environment .yaml files with necessary package requirements are provided in each subfolder.

System requirements to run these scripts with data of our size:
* \>128G RAM (most code was run on a Linux machine with 256G RAM)
* at least 1 GPU

Before running these scripts, please create a fresh conda environment and install the packages in "general_scanpy_scvi_environment.yml" via pip. This environment will work for most analyses, except those in folders where a sepeparate environmnet .yml file is listed (e.g. SATURN, LIANA, pertpy, etc.).
