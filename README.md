# Predicting Camera Trap Detection Zones  

This repository provides an R script to **predict camera trap detection zones** (effective detection distance and angle) for wildlife species based on species- and camera-specific traits. The approach is based on the paper:  

> Wiegers et al. (2025). *Inferring camera trap detection zones for rare species using species- and camera-specific traits: a meta-level analysis*. Remote Sensing in Ecology and Evolution.  

Accurate estimation of a camera trap’s **effective detection radius (EDD)** and **effective detection angle (EDA)** is essential for wildlife monitoring, especially when calculating densities of unmarked species. However, direct estimation requires large sample sizes of detection events, which is often not feasible for rare species.  
This script provides a modeling framework that combines **species traits, technical covariates, and limited empirical data** to predict detection zones in three practical scenarios.

---

## Features  

The script allows users to:  
1. **Estimate EDD and EDA for species without sufficient detection data** using mixed-effects models (manuscript scenario 1).    
2. **Combine limited empirical measurements with model-based predictions** to improve estimates for rare species (manuscript scenario 2).  

---

##  Repository Structure  

- **`predict_detection_zones.R`** → Main script for running detection zone predictions.  
- **`LOOCV_data 2025.csv`** → Training dataset for EDD model.  
- **`LOOCV_data_theta.csv`** → Training dataset for EDA model.  
- **`predict.csv`** → Example input file for predicting EDD.  
- **`predict_theta.csv`** → Example input file for predicting EDA.  
- **`example_measurements.csv`** → Example dataset of limited field measurements (used in hybrid scenario 2).  

---

## Dependencies  

The script requires the following R packages:  

```r
library(readr)
library(nlme)
library(R.utils)
library(Distance)   # for detection function fitting
library(plyr)       # for rbind.fill
library(dplyr)      # for joins and data manipulation
```

---

## Input Data  

Users must provide a CSV file (`predict.csv` or `predict_theta.csv`) with the following columns:  

- **body_mass** (kg)  
- **ct_height** (m)  
- **ct_brand** (factor: `Browning`, `Reconyx`, `Bushnell`, `Little Acorn`, `Spromise`)  
- **T** (snapshot interval in seconds, or `"contact"`)  
- **data.order** (factor: `Aves`, `Carnivora`, `Primates`, `Rodentia`, `Ungulata`)  
- **max_r_site** (right-truncated distance, m; for EDD)  
- **max_theta_site** (right-truncated angle, degrees; for EDA)  

---

## Usage  

1. Place (or edit) the prediction input file (`predict.csv` or `predict_theta.csv`) in the working directory.  
2. Open and run the R script `predict_detection_zones.R`.  
3. The script will:  
   - Train mixed-effects models on the provided LOOCV datasets.  
   - Generate predictions for your species/site combination.  
   - Save results with estimated detection radius/angle.  

---

### **Scenario 2 (Hybrid)**  
- Combines **limited empirical measurements** (e.g., <50 detections) with model predictions.  
- Uses weighted averaging between the two sources to improve precision.  
- Only implemented for **EDD only** since we show in the manuscript that this does not improve the EDA estimates.  

---

## 📧 Contact  

For questions or contributions, please contact:  
**J.N. Wiegers** – j.n.wiegers@uu.nl  

We would love to add more training data to the models so feel free to send me (peer-reviewed) data or measurements of detection zones.



