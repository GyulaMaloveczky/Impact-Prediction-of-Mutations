# -*- coding: utf-8 -*-


import pandas as pd
import numpy as np
import sklearn
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier

data_train = pd.read_csv('output/cleaned_train.csv')

data_train = data_train.drop_duplicates(subset= [ 'ref_aa', 'neighbor_minus_1','neighbor_minus_2', 'neighbor_minus_3', 'neighbor_minus_4', 'neighbor_plus_1','neighbor_plus_2', 'neighbor_plus_3', 'neighbor_plus_4' ])
data_train = data_train.drop_duplicates(subset= ['W_ratio', 'length', 'I_ratio', 'C_ratio', 'ref_aa','mut_aa'])

test_data = pd.read_csv('output/test.tsv', delimiter = '\t')
test_data = test_data.drop(columns = ['self'])




X_train = data_train.iloc[:,2:-1]
y_train = data_train.iloc[:,-1]


X_test = test_data.iloc[:,1:]



mask = ~(data_train.isin([ 'U', 'o', '+'])).any(axis=1)

# Apply the mask to both X_train and y_train
X_train = X_train[mask].reset_index(drop=True)
y_train = y_train[mask].reset_index(drop=True)


# Background frequencies of amino acids
amino_acid_frequencies = {
'A': 0.074,  # Alanine
'R': 0.042,  # Arginine
'N': 0.044,  # Asparagine
'D': 0.059,  # Aspartic Acid
'C': 0.033,  # Cysteine
'E': 0.058,  # Glutamic Acid
'Q': 0.037,  # Glutamine
'G': 0.074,  # Glycine
'H': 0.029,  # Histidine
'I': 0.038,  # Isoleucine
'L': 0.076,  # Leucine
'K': 0.072,  # Lysine
'M': 0.018,  # Methionine
'F': 0.04,  # Phenylalanine
'P': 0.050,  # Proline
'S': 0.0810,  # Serine
'T': 0.062,  # Threonine
'W': 0.013,  # Tryptophan
'Y': 0.033,  # Tyrosine
'V': 0.068,  # Valine
}



# Amino acids common in beta sheets



beta_sheet_aa = [
    "I",  # Isoleucine
    "V",  # Valine
    "F",  # Phenylalanine
    "Y",  # Tyrosine
    "W",  # Tryptophan
    "N",  # Asparagine
    "Q",  # Glutamine
    "T"   # Threonine
]

# Amino acids common in alpha helices
alpha_helix_aa = [
    "A",  # Alanine
    "E",  # Glutamic Acid
    "L",  # Leucine
    "M",  # Methionine
    "R",  # Arginine
    "K",  # Lysine
    "H",  # Histidine
    "S"   # Serine
]




for index, row in X_train.iterrows():
    # Ensure that 'neighbor_plus_X' and 'neighbor_minus_X' are iterable (like a list of strings)
# Initialize the sums for beta and alpha
    beta_count = 0
    alpha_count = 0
    
    # Iterate through the relevant neighbors
    for neighbor in ['ref_aa','neighbor_plus_1', 
                     'neighbor_minus_1','neighbor_minus_2','neighbor_plus_2']:
        
        
        # Check if the current neighbor is in the beta and alpha sets
    
        beta_count += row[neighbor] in beta_sheet_aa
        alpha_count += row[neighbor] in alpha_helix_aa
    
    # Assign the computed values to the DataFrame
    X_train.at[index, 'beta'] = beta_count - alpha_count
    X_train.at[index, 'alpha'] = alpha_count - beta_count

    x= row['mut_aa']
    z= row['ref_aa']
    
    X_train.at[index, 'aa_ratio_mut'] = row[f'{x}_ratio']/amino_acid_frequencies[x]
   
         
    X_train.at[index, 'aa_ratio_ref'] = row[f'{z}_ratio']/amino_acid_frequencies[z]
   
   
    

# For X_test
for index, row in X_test.iterrows():
    
    # Iterate through the relevant neighbors
    for neighbor in ['ref_aa','neighbor_plus_1', 'neighbor_plus_2',
                     'neighbor_minus_1', 'neighbor_minus_2']:
        
        # Check if the current neighbor is in the beta and alpha sets
        beta_count += row[neighbor] in beta_sheet_aa
        alpha_count += row[neighbor] in alpha_helix_aa
    
    # Assign the computed values to the DataFrame
    X_test.at[index, 'beta'] = beta_count//3
    X_test.at[index, 'alpha'] = alpha_count//3
    x= row['mut_aa']
    z= row['ref_aa']
    if True:
        X_test.at[index, 'aa_ratio_mut'] = row[f'{x}_ratio']/amino_acid_frequencies[x]
   # elif row[f'{x}_ratio']/amino_acid_frequencies[x] > 1:
      #  X_test.at[index, 'aa_ratio_mut'] = 1
  #  else:
       # X_test.at[index, 'aa_ratio_mut'] = 0
    if True:
         
         X_test.at[index, 'aa_ratio_ref'] = row[f'{z}_ratio']/amino_acid_frequencies[z]
    #else:
     #   X_test.at[index, 'aa_ratio_ref'] = 0
   ## if row[f'{z}_ratio']/amino_acid_frequencies[z] > 3:
     #   X_test.at[index, 'aa_ratio_ref'] = 3
   
    
  

# Label encoding each amino acid as their Eisenber hydrophobicity
hydrophobicity_dict_eisenberg = {
    'A': 0.62,   # Alanine
    'C': 0.29,   # Cysteine
    'D': -0.90,  # Aspartic acid
    'E': -0.74,  # Glutamic acid
    'F': 1.19,   # Phenylalanine
    'G': 0.48,   # Glycine
    'H': -0.40,  # Histidine
    'I': 1.38,   # Isoleucine
    'K': -1.50,  # Lysine
    'L': 1.06,   # Leucine
    'M': 0.64,   # Methionine
    'N': -0.78,  # Asparagine
    'P': 0.12,   # Proline
    'Q': -0.85,  # Glutamine
    'R': -2.53,  # Arginine
    'S': -0.18,  # Serine
    'T': -0.05,  # Threonine
    'V': 1.08,   # Valine
    'W': 0.81,   # Tryptophan
    'Y': 0.26,   # Tyrosine
    '*': -10     # Stop codon (arbitrary value)
}

X_train = X_train.replace(hydrophobicity_dict_eisenberg)
X_test = X_test.replace(hydrophobicity_dict_eisenberg)


# Amino acids common in beta sheets
# Amino acids common in beta sheets


params = {
    'criterion': 'entropy',           # or 'gini'
    'max_depth': 15,                  # Limit the depth to avoid overfitting
    'min_samples_split': 20,          # Prevent overfitting by setting a higher value
    'min_samples_leaf': 10,           # Each leaf should have at least 10 samples
    'max_features': 'sqrt',           # Use a subset of features at each split
    'max_leaf_nodes': 200,            # Limit the total number of leaf nodes
    'random_state': 42                # Ensure reproducibility
}

# Optimal parameters obtained from hyperparameter tuning
params ={'colsample_bytree': 0.6842052413176569, 'gamma': 4.546976958278643, 'learning_rate': 0.002341619497778759, 'max_depth': 14, 'min_child_weight': 1.0, 'n_estimators': 10000, 'reg_alpha': 0.2596870025918193, 'reg_lambda': 4.161182806030448, 'scale_pos_weight': 3.626399287792276, 'subsample': 0.9018541597688445}


# Creating the classifier
tree = XGBClassifier(**params)

# Using only a subset of features
selected_features = ['ref_aa', 'mut_aa', 'alpha']


X_train = X_train.loc[:,selected_features]
X_test = X_test.loc[:,selected_features]


y_train = y_train.map({'Benign':0, 'Pathogenic': 1})


X_train =X_train.apply(pd.to_numeric, errors='coerce')
X_test =X_test.apply(pd.to_numeric, errors='coerce')

# Training the model and creating prediction for the test set
tree.fit(X_train, y_train)
prediction = tree.predict_proba(X_test)


# Deriving and plotting importances
importances = tree.feature_importances_
indices = np.argsort(importances)[::-1]
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.title("Feature Importance")
plt.bar(range(X_train.shape[1]), importances[indices])
plt.xticks(range(X_train.shape[1]), X_train.columns[indices], rotation=90)
plt.savefig('output/importance.png')

final_df = pd.concat([test_data.iloc[:,0], pd.Series(prediction[:,1])], axis = 1)


# Saving the output as tsv
out_path = 'xgb_restults.tsv'
final_df.to_csv(out_path, sep='\t', index = False, header = False)
