# Pediatric Hyperinflammation Proteomics

This repository contains analysis code associated with the manuscript:

**Comparative plasma proteomics reveals shared and divergent inflammatory signatures in pediatric sepsis, MIS-C, and CAR-T cytokine release syndrome**
Steven D. Ham, Rawan Shraim, Caroline Diorio, Apoorva Babu, Edward M. Behrens, David T. Teachey, Sarah E. Henrickson, Robert B. Lindell

## Overview

Sepsis, multisystem inflammatory syndrome in children (MIS-C), and cytokine release syndrome (CRS) after CAR-T therapy are clinically distinct pediatric syndromes that may converge on shock, organ dysfunction, and systemic hyperinflammation. This project uses high-dimensional Olink proteomic data to evaluate whether these conditions represent a shared inflammatory state or distinct host-response programs.

The analysis integrates proximity extension assay-based proteomic data across pediatric cohorts of sepsis, MIS-C, CRS, COVID-19, and healthy controls. Plasma and serum measurements were harmonized using validated serum-to-plasma transformation factors. The final analytic dataset included 701 proteins from 150 pediatric participants.

## Repository contents

This repository includes code used to generate the analyses and figures for the manuscript, including:

* Harmonized proteomic data preprocessing
* Principal component analysis and UMAP projection
* Protein correlations with UMAP dimensions
* Hallmark pathway single-sample gene set variation analysis
* Heatmap and radar plot visualization of pathway enrichment
* XGBoost-based feature prioritization comparing MIS-C and severe CRS with sepsis

## Data availability

The clinical and proteomic datasets analyzed in this study are not included in this repository. Data may be available from the corresponding author upon reasonable request and in accordance with institutional review board approvals, data use agreements, and participant privacy protections.

This repository is intended to provide transparent and reproducible code for the analyses reported in the manuscript.

## Contact

For questions about the manuscript or analysis code, please contact:

Steven D. Ham, MD, MPH
Children’s Hospital of Philadelphia
[HamS1@chop.edu](mailto:HamS1@chop.edu)

Robert B. Lindell, MD
Children’s Hospital of Philadelphia
University of Pennsylvania Perelman School of Medicine
[LindellR@chop.edu](mailto:LindellR@chop.edu)
