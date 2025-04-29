# Investigating Different Weighting Approaches for Functional PCA

**Author:** Jonathan Cauchi  
**Degree:** B.Sc. (Hons.), Department of Statistics & OR, University of Malta  
**Date:** May 2025  

---

## 1. Overview

This repository houses all code and data pipelines for:

> **Investigating Different Weighting Approaches for Functional Principal Components Analysis**  
> Jonathan Cauchi, May 2025

We apply **Functional Data Analysis (FDA)**—with a focus on **Functional PCA (FPCA)**—to two real-world time series:

1. **Radiation levels** at Dutch KNMI weather stations  
2. **Daily closing prices** of assets in a hedge-fund portfolio  

Four weighting schemes are compared:

- Unweighted  
- Auto-Regressive (AR)  
- Rolling-variance  
- GARCH-based  

We demonstrate how each affects smoothing, eigencomponent estimation, and FPCA results.

---

## 2. Repository Layout

```text
/
├── data/
│   ├── radiation_levels.csv        # Raw KNMI radiation data
│   └── hedge_fund_returns.RDA      # Hedge-fund returns (R data file)
│
├── scripts/
│   ├── 01_preprocessing.R          # Imputation & weight construction
│   ├── 02_smoothing.R              # Basis setup & GCV-driven smoothing
│   ├── 03_fpca_analysis.R          # FPCA under each weighting scheme
│   └── 04_derivative_estimation.R  # First/second derivative estimation
│
├── figures/                        # Generated plots (.png, .pdf)
├── tables/                         # LaTeX-ready summary tables
└── README.md                       # This file
