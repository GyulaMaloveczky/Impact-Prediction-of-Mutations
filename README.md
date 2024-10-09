# Impact Prediction of Mutations
Impact prediction of mutations using a Gradient Boosting model, with training data created using NCBI's and UniProt's REST APIs.

# Predictive Model for Protein Mutations

This repository contains a machine learning model for predicting the pathogenicity of missense mutations in human proteins. The model is based on features extracted from amino acid sequences and is compared against two popular tools: **SIFT** and **PolyPhen2**. The training data was obtained from Uniprot and NCBI, using **REST API**.

## Table of Contents
- [Project Overview](#project-overview)
- [Key Features](#key-features)
- [Model Description](#model-description)
- [Data Acquisition](#data-acquisition)
- [Discussion](#discussion)

## Project Overview
The project applies a machine learning pipeline to classify missense mutations as **benign** or **pathogenic** using a dataset obtained from ClinVar. The model leverages amino acid sequence features, Bayesian hyperparameter tuning, and XGBoost for classification. It also compares performance against the **SIFT** and **PolyPhen2** tools.

## Key Features
- **Bayesian Hyperparameter Tuning:** Optimal model parameters were selected using Bayesian optimization. The tuning process reached saturation, meaning that further optimization would yield diminishing returns. Below is a plot showing how the model’s performance plateaued over successive iterations:

![tune_plot](https://github.com/user-attachments/assets/0728b69c-53c2-47a5-9ae1-f87866dcd543)




- **Comparison with SIFT and PolyPhen2:** The model's ROC curve is compared with these popular tools for protein mutation prediction for comparison the test dataset was used:
![comparison](https://github.com/user-attachments/assets/c3bbdf23-2058-4800-a2be-93e003b64bc4)


## Model Description
This project uses the **XGBoost** algorithm, a highly efficient gradient boosting model, to predict whether a given mutation is **pathogenic** or **benign**. The model's input features include:
- Neighboring amino acids around the mutation site
- Ratios of various amino acids in the protein sequence
- Hydrophobicity scores and other structural properties of the amino acids

The model undergoes hyperparameter tuning using **Bayesian optimization**, ensuring that it explores the parameter space efficiently to find the best settings.

### Selected Features
From the dataset, the most important features for prediction were selected using **Recursive Feature Elimination**:
- The original and the mutated amino acid and their hydrophobicity properties
- A score assigned based on the surrounding aminoacids likelyhood of appearing in beta sheets (negative score) vs. alpha helice (positive score)


## Data Acquisition
The training data for this project was obtained using **REST APIs** to query the **ClinVar** and **UniProt** databases.

### ClinVar Data Retrieval
We used the **NCBI Entrez API** to fetch missense mutation data from ClinVar, specifically filtering for pathogenic and benign variants:
- Pathogenic mutations: Retrieved using the query `clinsig_pathogenic[Properties] AND missense_variant[molecular consequence]`
- Benign mutations: Retrieved using the query `clinsig_benign[Properties] AND missense_variant[molecular consequence]`

To achieve this, the script uses the following steps:
1. **esearch API call**: This retrieves a list of ClinVar variant IDs matching the query.
2. **esummary API call**: For each variant ID, a summary of the variant (including its protein change) is retrieved.
3. **efetch API call**: Fetches more detailed data for each variant, including its specific amino acid changes.

### Protein Sequence Data Retrieval
For each mutation, the affected protein's sequence is retrieved using the **UniProt REST API**:
- The UniProt API is queried using the relevant mRNA or protein accession number from ClinVar.
- The sequence is parsed and used to calculate features such as amino acid ratios and neighboring residues around the mutation site.

The complete process is automated using Python's `requests` library to communicate with the APIs, and `xmltodict` and `pandas` to parse and manipulate the data.



## Discussion

The model can not outperform popular impact prediction models like SIFT or PolyPhen2, however it's predictive power is markedly higher then that of the Baseline model, which is based on the BLOSUM62 matrix.

