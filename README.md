# DGEA for Alzheimer's with ML for Biomarker Identification and Validation

This project explores the use of differential gene expression analysis (DGEA) and machine learning (ML) techniques for the identification and validation of biomarkers related to Alzheimer's disease (AD). It utilizes RNA-seq data and bioinformatics tools to detect differentially expressed genes (DEGs) and applies machine learning algorithms to assess the potential of these DEGs in distinguishing Alzheimer’s patients from control groups.

## Project Overview
The primary goal of this project is to evaluate different bioinformatics tools (such as DESeq2 and edgeR) for differential gene expression analysis and identify key biomarkers linked to Alzheimer's disease. After identifying DEGs, machine learning models such as Support Vector Machine (SVM) and Random Forest are employed to assess the predictive power of these biomarkers.

Key Objectives:
- Perform differential gene expression analysis using DESeq2 and edgeR.
- Conduct pathway enrichment analysis to identify biological processes linked to Alzheimer's.
- Apply machine learning algorithms to validate the identified biomarkers.

## Dataset
The dataset used in this project is the **GEO153873** dataset, which contains RNA-seq data for Alzheimer's disease patients and control subjects. The dataset is publicly available on the Gene Expression Omnibus (GEO) platform.

- **Dataset Source**: [GEO153873 on NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GEO153873)
- **Data Format**: CSV / Gene Expression Counts
- **Number of Samples**: [30]
- **Key Variables**: Gene Expression Levels, Sample Type (AD vs. Control)

## Installation
To run this project locally, clone the repository and install the necessary dependencies:

```bash
git clone https://github.com/your-username/DGEA-for-Alzheimer-s-with-ML-for-Biomarker-Identification-and-Validation.git
cd DGEA-for-Alzheimer-s-with-ML-for-Biomarker-Identification-and-Validation
