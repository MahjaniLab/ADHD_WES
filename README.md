## Overview
This repository contains codes used for the analyses in our manuscript:
"Rare Variant Analyses in Ancestrally Diverse Cohorts Reveal Novel ADHD Risk Genes."
doi: https://doi.org/10.1101/2025.01.14.25320294

In this study, we investigated the role of ultra-rare deleterious genetic variants in ADHD risk, using whole-exome sequencing (WES) and whole-genome sequencing (WGS) data from multiple datasets, including All of Us, SPARK, and other published cohorts. We identified 15 ADHD risk genes, performed pathway enrichment, protein-protein interaction (PPI), and ancestry-specific analyses, and examined shared genetic architecture with autism spectrum disorder (ASD).

This repository provides the scripts and workflows used in the analysis.


## Contents
1. Rare Variant Association and Gene Discovery
"Gamma_calculation_for_TADA.R", and "Run_TADA.R": Implements TADA (Transmission and De Novo Association Test) for rare variant association testing.

2. Differential Expression Analysis
"Run_differential_gene_expression_analysis_using_GEO.R": Analyzes ADHD risk gene expression using RNA-seq data from GEO.

3. Ancestry-Specific Analyses
Ancestry_analysis_15genes.R: Compares variant burden across ancestries in ADHD cases and population reference datasets.

4. ADHD-ASD Shared Genetic Architecture
Run_Heterogeneity_test_between_ADHD_and_ASD.R: Tests for genetic heterogeneity between ADHD and ASD.


## Data Availability
All of Us controlled tier data is available for registered institutions (https://allofus.nih.gov/).
SPARK dataset is available through Simons Foundation (https://sparkforautism.org/).
Summary statistics from Olfson et al. (2024) and Satterstrom et al. (2019) were used.


## Contact
For questions about the code, please contact:

Behrang Mahjani, PhD
Email: behrang.mahjani@mssm.edu

Seulgi Jung, PhD
Email: seulgi.jung@mssm.edu
