# DGEA for Alzheimer's with ML for Biomarker Identification and Validation

This project explores the use of differential gene expression analysis (DGEA) and machine learning (ML) techniques for the identification and validation of biomarkers related to Alzheimer's disease (AD). It utilizes RNA-seq data and bioinformatics tools to detect differentially expressed genes (DEGs) and applies machine learning algorithms to assess the potential of these DEGs in distinguishing Alzheimer’s patients from control groups.

## Table of Contents

- [Introduction](#introduction)
- [Project Structure](#project-structure)
- [Data Sources](#data-sources)
- [Analysis Workflow](#analysis-workflow)
- [Requirements](#requirements)
- [Installation](#installation)
  
## Introduction

Alzheimer’s disease (AD) is a neurodegenerative disorder that leads to cognitive decline. Understanding the molecular mechanisms underlying AD is critical for developing diagnostic biomarkers and therapeutic strategies. This project utilizes differential gene expression analysis (DGEA) and machine learning (ML) to identify and validate potential biomarkers related to AD, aiming to improve early diagnosis and understanding of the disease’s pathophysiology.

By using RNA-seq data from Alzheimer's patients and control groups, the project compares different tools for gene expression analysis and investigates the predictive power of biomarkers using machine learning algorithms like Support Vector Machine (SVM) and Random Forest.

## Project Structure

- **Data/**: Raw and processed data, including the GEO153873 Alzheimer's dataset and output files from differential expression analysis.
- **Scripts/**: R and Python script files for data processing, analysis, and machine learning model building.
- **Plots/**: Figures and visualizations, including heatmaps, EDA, and machine learning performance curves.
- **README.md**: Project overview and instructions.

## Data Sources

- **GEO153873**: A publicly available RNA-seq dataset for Alzheimer's disease from the Gene Expression Omnibus (GEO).
  - **Dataset Link**: [GEO153873](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GEO153873)
  - **Data Format**: CSV (gene expression counts)
  - **Number of Samples**: 30
 
## Analysis Workflow

**Data Import and Cleaning** 
- Load and preprocess the gene expression data.
- Normalize the data and filter out lowly expressed genes.
- Conducted in R

**Differential Gene Expression Analysis**
- Use DESeq2 and edgeR to identify differentially expressed genes (DEGs) between Alzheimer's and control samples.
- Conducted in R

**Pathway Enrichment Analysis**
- Conduct pathway enrichment analysis using DEGs to explore biological pathways related to neuroinflammation, synaptic dysfunction, and Alzheimer’s disease.
- Conducted in R

**Machine Learning Model Training and Testing**
- Train and evaluate machine learning models, such as SVM and Random Forest, to assess the predictive power of the identified biomarkers.
- Generate performance metrics (AUC, accuracy, F1 score) to compare model effectiveness.
- Conducted in Python

**Visualizations**
- Generate plots, such as box plots, heatmaps, and performance curves, to visualize results from the DGEA and machine learning models.
- Conducted in R and Python
  
## Requirements

- **R** (version 4.0 or higher)
  - R libraries: DESeq2, edgeR, limma, ggplot2, dplyr, pheatmap, reshape2, gridExtra, clusterProfiler, org.Hs.eg.db
- **Python** (version 3.6 or higher)
  - Python libraries: pandas, numpy, scikit-learn, matplotlib, seaborn, biopython, scipy

## Installation
Clone the repository and install the necessary dependencies:

```bash
git clone https://github.com/your-username/DGEA-for-Alzheimer-s-with-ML-for-Biomarker-Identification-and-Validation.git
cd DGEA-for-Alzheimer-s-with-ML-for-Biomarker-Identification-and-Validation
