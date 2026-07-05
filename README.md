# Optimizing Training Sets in Genomic Prediction based on Bayesian Optimization

[![Language: R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains the R implementation of the core algorithms developed in the master's thesis: *Optimizing Training Sets in Genomic Prediction based on Bayesian Optimization*.

To protect unpublished academic data and original genomic datasets, this repository provides the complete algorithmic framework, simulation pipelines, and parallel computing architecture, accompanied by randomly generated **dummy data** for logical verification. 

The core simulation pipeline adopts a data-driven approach. It estimates variance components from empirical phenotypic data and evaluates the performance of various Bayesian Optimization acquisition functions (e.g., EI, UCB) across different training set sizes.

## Authors
* **Kun-Hong Liao** M.S., Biometry Division, Department of Agronomy, National Taiwan University
* **Advisor**: Dr. Chen-Tuo Liao

## Tech Stack & Dependencies
This project is built primarily with R for statistical computing and predictive modeling, leveraging the `future` framework for parallel processing. 
Please ensure the following key packages are installed:
* `BGLR` (Bayesian Generalized Linear Regression)
* `future`, `future.apply` (Parallel and distributed processing)
* `dplyr`, `purrr`, `tibble` (Data manipulation)
* `progressr`, `cli` (CLI progress visualization)
* `MASS`

## Repository Structure
The repository is modularized, separating core functions from execution scripts:

```text
├── src/                                  # Core algorithms and custom functions
│   ├── AcquisitionFunction.R             # Defines acquisition functions (e.g., EI, UCB) for Bayesian Optimization
│   ├── DesignOptimalTrainingFunction.R   # Core logic for constructing the optimized training set
│   ├── GvFunction.R                      # Defines A-optimality-like criteria
│   ├── SimulationFunction.R              # Main function for parallel simulation and metric calculation
│   └── TwoStage.R                        # Implements the two-stage prediction and evaluation model (discussion)
│
└── README.md
```

## Installation

Install required R packages:

```r
install.packages(c(
  "BGLR",
  "future",
  "future.apply",
  "progressr",
  "cli",
  "MASS",
  "dplyr",
  "purrr",
  "tibble"
))
```

Load core scripts: 

```r
source("src/AcquisitionFunction.R")
source("src/GvFunction.R")
source("src/DesignOptimalTrainingFunction.R")
source("src/SimulationFunction.R")
```

## Input Data Structure
### Additive kinship matrix
```r
kinshipA
```
A square symmetric matrix representing additive genetic relationships among individuals.

### Dominance kinship matrix (Optional)
```r
kinshipD
```
A square symmetric matrix representing dominance genetic relationships among individuals.

### Requirements:

* Must be square matrices
* Must have identical dimensions
* Must share the same row/column ordering
* Row/column names correspond to individuals

### Example:
```r
dim(kinshipA)
# 500 x 500
```

## Main Workflow
1. Simulate phenotypic values using kinship matrices
2. Split individuals via cross-validation
3. Fit RKHS models using `BGLR`
4. Compute acquisition scores (EI, UCB) and GV
5. Rank candidate individuals

## Example Usage
### Additive Kernel Only
```r
result <- OPTtrain(
  kinshipA = KA,
  nsim = 100,
  folds = 5,
  h = 0.5,
  mu = 100,
  varA = 20
)
```

### Additive and Dominance Kernel Only
```r
result <- OPTtrain(
  kinshipA = KA,
  kinshipD = KD,
  nsim = 100,
  folds = 5,
  h = 0.5,
  mu = 100,
  varA = 20,
  rho = 0.5
)
```

### Output
The function returns a list containing:

| Object | Description |
|--------|-------------|
| EI     | Augmented Expected Improvement ranking |
| UCB    | Augmented Upper Confidence Bound ranking |
| GV     | Genetic variance-based baseline ranking |

Example:
```r
head(result$EI)
```

