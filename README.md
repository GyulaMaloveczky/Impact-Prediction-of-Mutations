# Impact Prediction of Mutations
Impact prediction of mutations using a Gradient Boosting model, with training data created using NCBI's and UniProt's REST APIs.

# Predictive Model for Protein Mutations

This repository contains a machine learning model for predicting the pathogenicity of missense mutations in human proteins. The model is based on features extracted from amino acid sequences and is compared against two popular tools: **SIFT** and **PolyPhen2**. The training data was obtained from UniProt and NCBI using **REST APIs**.

## Table of Contents
- [Abstract](#abstract)
- [Results](#results)
- [Conclusions](#conclusions)
- [Model Description](#model-description)
- [Data Acquisition](#data-acquisition)


## Abstract
The project applies a machine learning pipeline to classify missense mutations as **benign** or **pathogenic** using a dataset obtained from ClinVar. The model leverages amino acid sequence features, Bayesian hyperparameter tuning, and XGBoost for classification. It also compares performance against **SIFT** and **PolyPhen2**.

## Results
- **Bayesian Hyperparameter Tuning:** Optimal model parameters were selected using Bayesian optimization. The tuning process reached saturation, meaning that further optimization would yield diminishing returns. Below is a plot showing how the model’s performance plateaued over successive iterations, evaluated using 5-fold Cross Validation:

![tune_plot](https://github.com/user-attachments/assets/0728b69c-53c2-47a5-9ae1-f87866dcd543)

- **Comparison with SIFT and PolyPhen2:** The model's ROC curve is compared with these popular tools for protein mutation prediction:

![comparison](https://github.com/user-attachments/assets/c3bbdf23-2058-4800-a2be-93e003b64bc4)

## Conclusions
The model does not outperform popular impact prediction models like SIFT or PolyPhen2, however, its predictive power is markedly higher than that of the Baseline model, which is based on the BLOSUM62 matrix.

## Model Description
This project uses the **XGBoost** algorithm, a highly efficient gradient boosting model, to predict whether a given mutation is **pathogenic** or **benign**. The model's input features include:
- Neighboring amino acids around the mutation site
- Ratios of various amino acids in the protein sequence
- Hydrophobicity scores (each amino acid was label encoded using their Eisenberg Hydrophobicity score, as this uniquely labels each amino acid compared to other metrics)
- Engineered features (see below)

The model undergoes hyperparameter tuning using **Bayesian optimization**, ensuring efficient exploration of the parameter space to find the best settings.

### Selected Features
Engineered features were derived from the mutation and surrounding sequence. For each mutation, we calculated the relative frequency of the reference and mutated amino acid in the sequence (normalized by background frequencies). We also computed an alpha score for the amino acids proximal to the mutation site. This score was derived by subtracting the count of amino acids frequently present in beta sheets from those in alpha helices. This was performed in the immediate proximity (+/-2 amino acids) of the mutation position.

While this is a oversimplified metric of secondary structure, it was employed for the sake of simplicity.

From the dataset, the most important features for prediction were selected using **Recursive Feature Elimination**:
- The original and mutated amino acid, including their hydrophobicity properties
- The relative frequency of the mutated and reference amino acids in the sequence
- A score based on the likelihood of surrounding amino acids appearing in beta sheets (negative score) versus alpha helices (positive score)

## Data Acquisition
The training data for this project was obtained using **REST APIs** to query the **ClinVar** and **NCBI's protein** or **UniProt** databases.

### ClinVar Data Retrieval
We used the **NCBI Entrez API** to fetch missense mutation data from ClinVar, specifically filtering for pathogenic and benign variants and applying the filters mentioned in the project manual:
- Pathogenic mutations query: `clinvar_snp[Filter] AND clinsig_pathogenic[Properties] AND missense_variant[molecular consequence] NOT moi_mitochondrial[prop] AND single_gene[prop]`
- Benign mutations query: `clinvar_snp[Filter] AND clinsig_benign[Properties] AND missense_variant[molecular consequence] NOT moi_mitochondrial[prop] AND single_gene[prop]`

To achieve this, the script follows these steps:
1. **esearch API call**: Retrieves a list of ClinVar variant IDs matching the query.
2. **esummary API call**: Retrieves a summary for each variant ID, including its protein change.
3. **efetch API call**: Fetches detailed data for each variant, including the corresponding protein's ID.

### Protein Sequence Data Retrieval
For each mutation, the affected protein's sequence is retrieved using the **NCBI protein database** or, when not found, the **UniProt REST API**:
- The NCBI/UniProt API is queried using the relevant mRNA or protein accession number from ClinVar.
- Validation is performed to ensure that the reference amino acid (from ClinVar) matches the amino acid at the corresponding position in the sequence.
- The sequence is parsed and used in further feature engineering.

The complete process is automated using Python's `requests` library to communicate with the APIs, and `xmltodict` and `pandas` for data parsing and manipulation.

### Cleaning
Further data cleaning was done to remove: Entries with the same genomic location or ones that were close to each other (15 nucleic acid range). Training examples with the same genomic locations as test examples were also removed to avoid data leakage.
### Final Train data
The training data included 29548 examples in total. The train data is balanced containing 50% Bening and 50% Pathogenic labelled entries.
