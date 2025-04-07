# Predictive Model for Protein Mutations

## Table of Contents
- [Abstract](#abstract)
- [Results](#results)
- [Conclusions](#conclusions)
- [Methods](#methods)
- [Model Description](#model-description)
- [Data Acquisition](#data-acquisition)

## Abstract
The project applies a machine learning pipeline to classify missense mutations as **benign** or **pathogenic** using a dataset obtained from ClinVar [<a href="#ref1">1</a>]. The model leverages amino acid sequence features, Bayesian hyperparameter tuning, and XGBoost for classification. It also compares performance against **SIFT** [<a href="#ref2">2</a>] and **PolyPhen2** [<a href="#ref3">3</a>].

## Results
- **Bayesian Hyperparameter Tuning:** Optimal model parameters were selected using Bayesian optimization [<a href="#ref4">4</a>]. The tuning process reached saturation, meaning that further optimization would yield diminishing returns. Below is a plot showing how the model’s performance plateaued over successive iterations, evaluated using 5-fold Cross Validation:

![tune_plot](https://github.com/user-attachments/assets/0728b69c-53c2-47a5-9ae1-f87866dcd543)

- **Comparison with SIFT and PolyPhen2:** The model's ROC curve is compared with these popular tools for protein mutation prediction:

![comparison](https://github.com/user-attachments/assets/c3bbdf23-2058-4800-a2be-93e003b64bc4)

## Conclusions
The model does not outperform popular impact prediction models like SIFT [<a href="#ref2">2</a>] or PolyPhen2 [<a href="#ref3">3</a>], however, its predictive power is markedly higher than that of the Baseline model, which is based on the BLOSUM62 matrix [<a href="#ref5">5</a>].


## Methods
## Model Description
This project uses the **XGBoost** algorithm, a highly efficient gradient boosting model, to predict whether a given mutation is **pathogenic** or **benign** [<a href="#ref6">6</a>]. The model's input features include:
- Neighboring amino acids around the mutation site
- Ratios of various amino acids in the protein sequence
- Hydrophobicity scores (each amino acid was label encoded using their Eisenberg Hydrophobicity score, as this uniquely labels each amino acid) [<a href="#ref7">7</a>]
- Engineered features (see below)

The model undergoes hyperparameter tuning using **Bayesian optimization** [<a href="#ref4">4</a>], ensuring efficient exploration of the parameter space to find the best parameters.

### Selected Features
Engineered features were derived from the mutation and surrounding sequence. For each mutation, we calculated the relative frequency of the reference and mutated amino acid in the sequence (normalized by background frequencies). We also computed an alpha score for the amino acids proximal to the mutation site. This score was derived by subtracting the count of amino acids frequently present in beta sheets from those in alpha helices (frequently present was defined by having higher than 1 propensit according to: [<a href="#ref8">8</a>]). This was performed in the immediate proximity (+/-2 amino acids) of the mutation position.

While this is an oversimplified metric of secondary structure, it was employed for the sake of simplicity.

From the dataset, the most important features for prediction were selected using **Recursive Feature Elimination**:
- The original and mutated amino acid, including their hydrophobicity properties
- The relative frequency of the mutated and reference amino acids in the sequence
- A score based on the likelihood of surrounding amino acids appearing in beta sheets (negative score) versus alpha helices (positive score)

## Data Acquisition
The training data for this project was obtained using **REST APIs** to query the **ClinVar** [<a href="#ref1">1</a>] and **NCBI's protein** or **UniProt** [<a href="#ref9">9</a>] databases.

### ClinVar Data Retrieval
We used the **NCBI Entrez API** to fetch missense mutation data from ClinVar [<a href="#ref1">1</a>], specifically filtering for pathogenic and benign variants and applying the filters mentioned in the project manual:
- Pathogenic mutations query: `clinvar_snp[Filter] AND clinsig_pathogenic[Properties] AND missense_variant[molecular consequence] NOT moi_mitochondrial[prop] AND single_gene[prop]`
- Benign mutations query: `clinvar_snp[Filter] AND clinsig_benign[Properties] AND missense_variant[molecular consequence] NOT moi_mitochondrial[prop] AND single_gene[prop]`

To achieve this, the script follows these steps:
1. **esearch API call**: Retrieves a list of ClinVar variant IDs matching the query.
2. **esummary API call**: Retrieves a summary for each variant ID, including its protein change.
3. **efetch API call**: Fetches detailed data for each variant, including the corresponding protein's ID.

### Protein Sequence Data Retrieval
For each mutation, the affected protein's sequence is retrieved using the **NCBI protein database** or, when not found, the **UniProt REST API** [<a href="#ref9">9</a>]:
- The NCBI/UniProt API is queried using the relevant mRNA or protein accession number from ClinVar.
- Validation is performed to ensure that the reference amino acid (from ClinVar) matches the amino acid at the corresponding position in the sequence.
- The sequence is parsed and used in further feature engineering.

The complete process is automated using Python's `requests` library to communicate with the APIs, and `xmltodict` and `pandas` for data parsing and manipulation.

### Cleaning
Further data cleaning was done to remove entries with the same genomic location or ones that were close to each other (within 15 nucleic acids). Training examples with the same genomic locations as test examples were also removed to avoid data leakage.

### Final Train Data
The training data included 29,548 examples in total. The train data is balanced, containing 50% benign and 50% pathogenic labeled entries.

## References
<a id="ref1">1. ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/</a>  
<a id="ref2">2. SIFT: Ng, P. C., & Henikoff, S. (2003). SIFT: Predicting amino acid changes that affect protein function. *Nucleic Acids Research*, 31(13), 3812-3814. https://doi.org/10.1093/nar/gkg509</a>  
<a id="ref3">3. PolyPhen-2: Adzhubei, I. A., Schmidt, S., Peshkin, L., Ramensky, V. E., Gerasimova, A., Bork, P., Kondrashov, A. S., & Sunyaev, S. R. (2010). A method and server for predicting damaging missense mutations. *Nature Methods*, 7(4), 248-249. https://doi.org/10.1038/nmeth0410-248</a>  
<a id="ref4">4. Bayesian Optimization: Snoek, J., Larochelle, H., & Adams, R. P. (2012). Practical Bayesian optimization of machine learning algorithms. *Advances in Neural Information Processing Systems*, 25.</a>  
<a id="ref5">5. BLOSUM62: Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences*, 89(22), 10915-10919.</a>  
<a id="ref6">6. XGBoost: Chen, T., & Guestrin, C. (2016). XGBoost: A scalable tree boosting system. *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*.</a>  
<a id="ref7">7. Eisenberg Hydrophobicity Scale: Eisenberg, D., Schwarz, E., Komaromy, M., & Wall, R. (1984). Analysis of membrane and surface protein sequences with the hydrophobic moment plot. *Journal of Molecular Biology*, 179(1), 125-142.</a>  
<a id="ref8">8. Chou, P. Y., & Fasman, G. D. (1974). Prediction of protein conformation. *Biochemistry*, **13**(2), 222–245.</a>
<a id="ref9">9. UniProt: https://www.uniprot.org/</a>


