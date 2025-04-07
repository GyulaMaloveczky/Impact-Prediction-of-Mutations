
import pandas as pd
from skopt import BayesSearchCV
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import make_scorer, roc_auc_score
import matplotlib.pyplot as plt
import numpy as np
import sys
from xgboost import XGBClassifier
from skopt.space import Real, Integer

def prepare_data():   
    
    # Creating the corrected dictionary with background frequencies in decimal format
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



    
    # Importing the train data and dropping anything that is included in the test data
    data_train = pd.read_csv('output/cleaned_train.csv')
    #data_test = pd.read_csv(r'\\wsl.localhost\Ubuntu\home\gyuszi\data/HGVS_2014_benchmark.tsv', delimiter= '\t')
 
    
    
    # Dropping U and * containing records (They are not in the test data, U is Selenocysteine and * is a stop codon)
    data_train = data_train.where(data_train != '*')
    data_train = data_train.where(data_train != 'U')
    data_train = data_train.where(data_train != 'o')
    data_train = data_train.where(data_train != '+')
    
    data_train = data_train.dropna()
    
    
    # List of amino acids common in beta helices
    beta_sheet_aa = [
     "I",  # Ile
     "V",  # Val
     "F",  # Phe
     "Y",  # Tyr
     "W",  # Trp
     "L"   # Leu
     "C"   # Cys
     "Q"   #Gln
     "T"   # Thr
     "M"   # Met
]

# Amino acids common in alpha helices
alpha_helix_aa = [
    "A",  # Ala
    "Q",  # Gln
    "E",  # Glu
    "L",  # Leu
    "M",  # Met
    "R",  # Arginine
    "K",  # Lys
    "H",  # His
    "S"   # Serine
    "F",  # Phe
    "W"   # Trp
    "V"   #Val
]

    # Label encoding the amino acids by their hidrophobicity (Eisenberg hydrophobicity was used as it gives a unique value for each aa, need for proper label encoding) and the stop codon is labeled as 1000
    
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
    
    
    # Splitting labels from features
    

    data_train = data_train.drop_duplicates(subset= [ 'ref_aa', 'neighbor_minus_1','neighbor_minus_2', 'neighbor_minus_3', 'neighbor_minus_4', 'neighbor_plus_1','neighbor_plus_2', 'neighbor_plus_3', 'neighbor_plus_4' ])
    data_train = data_train.drop_duplicates(subset= ['W_ratio', 'length', 'I_ratio', 'C_ratio'])
   
    X_train = data_train.iloc[:,2:-1]
    y_train = data_train.iloc[:, -1].replace({'Benign' : 1, 'Pathogenic' : 0})
   
    for index, row in X_train.iterrows():
        beta_count = 0
        alpha_count = 0
        for neighbor in ['ref_aa','neighbor_plus_1',
                         'neighbor_minus_1', 'neighbor_plus_2', 'neighbor_minus_2']:
            
            # Check if the current neighbor is in the beta and alpha sets
            beta_count += row[neighbor] in beta_sheet_aa
            alpha_count += row[neighbor] in alpha_helix_aa
        
        # Assign the computed values to the DataFrame
        X_train.at[index, 'beta'] = beta_count - alpha_count
        X_train.at[index, 'alpha'] = alpha_count - beta_count
        x= row['mut_aa']
        z= row['ref_aa']
        if True: #row[f'{x}_ratio']/amino_acid_frequencies[x] > 1:
            X_train.at[index, 'aa_ratio_mut'] = row[f'{x}_ratio']/amino_acid_frequencies[x]
        else:
            X_train.at[index, 'aa_ratio_mut'] = 0
       # if row[f'{x}_ratio']/amino_acid_frequencies[x] > 2:
           # X_train.at[index, 'aa_ratio_mut'] = 2
            
        if True: #row[f'{z}_ratio']/amino_acid_frequencies[z] > 1:
            
            X_train.at[index, 'aa_ratio_ref'] = row[f'{z}_ratio']/amino_acid_frequencies[z]
        else:
            X_train.at[index, 'aa_ratio_ref'] = 0
        #if row[f'{z}_ratio']/amino_acid_frequencies[z] > 2:
            #X_train.at[index, 'aa_ratio_ref'] = 2
    X_train = X_train.replace(hydrophobicity_dict_eisenberg)
    X_train = X_train.loc[:,['ref_aa', 'mut_aa', 'alpha']]

    return X_train, y_train
    
    
    
    
# A class to save data ROC AUCs during the hyperparameter tuning. This will enable us to plot it later
class SaveScoreCallback:
    def __init__(self):
        self.scores = []
        
    def __call__(self, result):
        # Store the current best score at each iteration
        self.scores.append(-result['fun'])  # 'fun' holds the best score from the current iteration


# Function for hyperparameter tuning/ searching the hyperparameter space
def search(X,y):
    X = X.apply(pd.to_numeric, errors='coerce')
    # Creating a decision tree classifier
    tree = XGBClassifier()

    # Hyperparameter space for hyperparameter tuning 
    

    param_space = {
    'n_estimators': Integer(50, 10000),               # Number of boosting rounds (trees)
    'max_depth': Integer(2, 15),                     # Maximum depth of a tree
    'learning_rate': Real(0.001, 0.5, prior='log-uniform'),  # Learning rate (shrinkage)
    'subsample': Real(0.5, 1.0),                     # Fraction of samples used for fitting individual trees
    'colsample_bytree': Real(0.3, 1.0),              # Fraction of features used for fitting individual trees
    'gamma': Real(0, 10),                            # Minimum loss reduction required to make a split
    'reg_lambda': Real(1e-3, 10, prior='log-uniform'),  # L2 regularization term on weights
    'reg_alpha': Real(1e-3, 10, prior='log-uniform'),  # L1 regularization term on weights
    'min_child_weight': Real(1, 100),                 # Minimum sum of instance weight (hessian) needed in a child
    'scale_pos_weight': Real(0.1, 10)                # Control the balance of positive and negative weights
    }

    # Stratified k-fold makes sure the train test splits are right
    stratified_kfold = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    
    # Model for hyperparameter tuning, bayes search is more effective in searching a big hyperparameter space then other models. High iteration was set because of the size of the space.
    
    bayes_search = BayesSearchCV(tree, param_space, n_iter=200, cv=stratified_kfold, n_jobs=-1, verbose=0, random_state=42, scoring= make_scorer(roc_auc_score, needs_proba=True) )
    save_score_cb = SaveScoreCallback()
    bayes_search.fit(X, y, callback = save_score_cb)
    tree = XGBClassifier(**bayes_search.best_params_)
    tree.fit(X, y)
  
    # Plotting the feature importances
    importances = tree.feature_importances_
    indices = np.argsort(importances)[::-1]
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.title("Feature Importance")
    plt.bar(range(X.shape[1]), importances[indices])
    plt.xticks(range(X.shape[1]), X.columns[indices], rotation=90)
    plt.savefig('output/importance.png')
    print("Best hyperparameters:", bayes_search.best_params_)
    return save_score_cb.scores
def plot(scores): 
    # Plotting the hyperparameter tuning process. X: iterations Y: Score (ROC AUC)
    plt.close()
    plt.clf()
    plt.plot(np.arange(len(scores)), scores, label='ROC AUC')
    plt.xlabel('Iteration')
    plt.ylabel('ROC AUC')
    plt.title('Bayesian Search ROC AUC Progress')
    plt.legend()
    plt.savefig('output/tune_plot.png')
    plt.close()



def main():
    X, y = prepare_data() # Preprocesses the data
    
    scores = search(X, y)  # prints best parameters, returns scores for iterations
    plot(scores) #plots the scores against iterations

if __name__ == '__main__':
    main()



