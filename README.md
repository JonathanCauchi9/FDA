# Investigating Different Weighting Approaches for Functional PCA

**Author:** Jonathan Cauchi  
**Degree:** B.Sc. (Hons.), Department of Statistics & OR, University of Malta  
**Date:** May 2025  

---

# Overview

This repository houses all code and data pipelines for:

> **Investigating Different Weighting Approaches for Functional Principal Components Analysis**  
> Jonathan Cauchi, May 2025

We apply **Functional Data Analysis (FDA)**—with a focus on **Functional PCA (FPCA)**—to two real-world time series:

1. **Radiation levels** at Dutch KNMI weather stations  
2. **Daily closing prices** of assets in a hedge-fund portfolio

The folders in this repository contain the following material:

1. **Raw Data**  
   Contains all the original datasets used in this dissertation.

2. **Figures & Results**  
   Includes the relevant figures and results generated during the analysis.

3. **R Scripts**  
   Houses the R scripts used to produce the figures and results.
   
Four weighting schemes are compared:

- Unweighted  
- Auto-Regressive (AR)  
- Rolling-variance  
- GARCH-based  

We demonstrate how each affects smoothing, eigencomponent estimation, and FPCA results.

---
