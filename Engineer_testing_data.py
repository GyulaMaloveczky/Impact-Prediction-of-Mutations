# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 17:18:49 2024

@author: gyusz
"""

import pandas as pd 
import requests
import xmltodict
import argparse
import xml.etree.ElementTree as ET
import sys
import time


hgvs_df = pd.read_csv('data/vep/HGVS_2014_VEP_baseline.tsv', delimiter= '\t', header= 0)

def get_ids():
    ids = []
    change = []
    poses = []
    for i in range(0,40):
        print(i)
        
        
        
        # Define the Ensembl VEP API endpoint for batch querying
        VEP_API = "https://rest.ensembl.org/vep/human/hgvs"
    
        # List of HGVS IDs
        hgvs_ids = list(hgvs_df.iloc[i*10:i*10+10,0])
        len(hgvs_ids)
        
        
        # Set headers for the API request
        headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
    
        # Create the payload with multiple HGVS IDs
        payload = {
            "hgvs_notations": hgvs_ids
        }
    
        # Make a POST request to the VEP API for batch processing
        response = requests.post(VEP_API, headers=headers, json=payload)
        
        # Check if the request was successful
        if response.status_code == 200:
            vep_data = response.json()
            print(len(vep_data))
            # Iterate through the response for each HGVS ID
            for k,result in enumerate(vep_data):
                for h in result['transcript_consequences']:
                    if 'protein_end' in h and 'amino_acids' in h and 'transcript_id' in h:
                        # Extract the first transcript consequence (or you can loop through if needed)
                        transcript_consequence = h
                        
                        # Extract and print the relevant data
                        amino_acids = transcript_consequence.get('amino_acids', 'N/A')
                        transcript_id = transcript_consequence.get('transcript_id', 'N/A')
                        pos = transcript_consequence.get('protein_end', 'N/A')
                        
                            
                    
                        change.append(amino_acids)
                        ids.append(transcript_id)
                        poses.append(pos)
                        break
                    
                    
                    
        else:
            print('response issue')
            sys.exit()
            
      
    return ids, poses, change
def fetch_protein_sequence_pdb(mrna):
    # Search for the protein using the UniProt API
    search_url = f"https://rest.uniprot.org/uniprotkb/search?query={mrna}&format=txt"
    response = requests.get(search_url)

    if response.status_code == 200:
        # Parse the response to get the UniProt ID
        text = response.text
        seq_start = text.find("SQ")
        sequence_lines = text[seq_start:].splitlines()[1:]  # Skip the 'SQ' line
        sequence = ''.join(line.strip() for line in sequence_lines if line and not line.startswith('//'))
        sequence = sequence.replace(' ', '')
        return sequence
        




def get_seq(mrna):
# Define the URL for the Entrez efetch API
    succes = False
    while not succes:
        params_0 = { 'term' : mrna,
                  'retmode' : 'json'
        }
        response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein', params_0)
        if response.status_code != 200:
            print('issue with protein seq response')
            continue
        prot_list = response.json()['esearchresult']['idlist']
        
        seq_list = []
        succes = True
    
        
          
            
    if len(prot_list) == 0:
        seq_list.append(fetch_protein_sequence_pdb(mrna))
        
    
        
    
    else:
        for prot in prot_list:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
            
            # Set up parameters for the request
            params_seq = {
                
                'id': prot,
                'retmode': 'fasta',  
                'rettype': 'fasta'
            }
            
            
            succes = False
            while not succes:
                response = requests.get(url, params_seq)
                if response.status_code == 200:
                    
                    succes = True
            
            fasta = response.text
            seq = ''.join(fasta.split('\n')[1:]).strip()
                
            seq_list.append(seq)
        
        
    return seq_list
        





def create_rows(dataframe):
    mismatch_count = 0
    shorter_seq_count = 0
    mrna, poses, amino_acids = get_ids()
    print(len(mrna))
    print(len(poses))
    print(len(amino_acids))
    
    
    for num, i in enumerate(mrna):
        
        print(num)
        
        
        if num%10 == 0:
            print(f"Current number added is {num}")
            print(f"Mismatches between seq and clinvar data : {mismatch_count}")
            print(f"Sequences that were shorter than expected: {shorter_seq_count}")
        
        
        
        seq_list = get_seq(i)
        print(len(seq_list))

        for seq in seq_list:
            
            if seq_list == 0:
                print('Missing sequence or error with request')
                continue
            
            length = len(seq)
            if len(seq) > 0:
                
                aa_ratios = []
                for aa in ['R','H','K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']:
                    aa_ratios.append(seq.count(aa)/length)
                    
            
            
            
            # Transforming the index used in the tsvs for the course
            id = hgvs_df.iloc[num,0]
            
            
            #getting the reference aa, mutated aa and position of the longest protein associated with the mutation
            
            
            
            aa_pos = int(poses[num])
            ref_aa = amino_acids[num][0]
            mut_aa = amino_acids[num][-1]
            
                
            
            if aa_pos == 0:
                print('Not a mismatch mutation')
                continue
            
            if length < aa_pos :
    
            
                
                shorter_seq_count +=1 
                continue
            
            
            
            if ref_aa != seq[aa_pos-1]:
                
            
                continue
            print('Found a matching seq')
            
            break
        print('adding neighbors')
        neighbors = []
        for x in range(-5,5):
            if aa_pos + x < length:
                neighbors.append(seq[aa_pos + x])
            else:
                neighbors.append(0)
                
                    
        mut_features = pd.Series([id, ref_aa, mut_aa, *neighbors, *aa_ratios, len(seq)], index = ["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", "neighbor_minus_2", "neighbor_minus_1", "self", "neighbor_plus_1","neighbor_plus_2", "neighbor_plus_3", "neighbor_plus_4", "neighbor_plus_5", 'R_ratio','H_ratio','K_ratio', 'D_ratio', 'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 'G_ratio', 'P_ratio', 'A_ratio', 'V_ratio', 'I_ratio', 'L_ratio', 'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 'length' ])
        
        dataframe = pd.concat([dataframe, mut_features.to_frame().T]) 
        
         
    return dataframe
        

# Example usage:
# clinvar_id = "12345"  # Replace with an actual ClinVar ID
# print(get_protein_sequence_from_clinvar(clinvar_id))


def main():
    
    #parser = argparse.ArgumentParser()
    #parser.add_argument('-o', dest='out_path')
    #output = parser.parse_args().out_path
    #creating the dataframe with the correct columns
    dataframe = pd.DataFrame(columns= ["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", "neighbor_minus_2", "neighbor_minus_1", "self", "neighbor_plus_1","neighbor_plus_2", "neighbor_plus_3", "neighbor_plus_4", "neighbor_plus_5", 'R_ratio','H_ratio','K_ratio', 'D_ratio', 'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 'G_ratio', 'P_ratio', 'A_ratio', 'V_ratio', 'I_ratio', 'L_ratio', 'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 'length' ])
    
    dataframe = create_rows(dataframe)
    output = 'test.tsv'
            
            
            
    dataframe.to_csv(output, sep='\t', index=False)
        

if __name__ == "__main__":
    main()