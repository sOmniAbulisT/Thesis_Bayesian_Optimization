# Optimizing Training Sets in Genomic Prediction based on Bayesian Optimization

[![Language: R](https://img.shields.io/badge/Language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains the R implementation of the core algorithms developed in the master's thesis: *Optimizing Training Sets in Genomic Prediction based on Bayesian Optimization*.

To protect unpublished academic data and original genomic datasets, this repository provides the complete algorithmic framework, simulation pipelines, and parallel computing architecture, accompanied by randomly generated **dummy data** for logical verification. 

The core simulation pipeline adopts a data-driven approach. It estimates variance components from empirical phenotypic data and evaluates the performance of various Bayesian Optimization acquisition functions (e.g., EI, UCB) across different training set sizes.

## Authors
* ** Kun-Hong Liao ** M.S., Biometry Division, Department of Agronomy, National Taiwan University
* **Advisor**: Dr. Chen-Tuo Liao

## Tech Stack & Dependencies
This project is built primarily with R for statistical computing and predictive modeling, leveraging the `future` framework for parallel processing. 
Please ensure the following key packages are installed:
* `BGLR` (Bayesian Generalized Linear Regression)
* `future`, `future.apply` (Parallel and distributed processing)
* `dplyr`, `purrr`, `tibble` (Data manipulation)
* `progressr`, `cli` (CLI progress visualization)
* `MASS`