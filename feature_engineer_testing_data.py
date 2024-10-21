# -*- coding: utf-8 -*-

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

# Load the HGVS dataset into a pandas DataFrame
hgvs_df = pd.read_csv('data/vep/HGVS_2014_VEP_baseline.tsv', delimiter= '\t', header= 0)

def get_ids():
    # Function to query Ensembl VEP API and retrieve transcript details
    ids = []  # List to store transcript IDs
    change = []  # List to store amino acid changes
    poses = []  # List to store protein positions

    # Loop to process batches of 10 HGVS IDs at a time (40 iterations total)
    for i in range(0, 40):
        print(i)

        # Define the Ensembl VEP API endpoint for batch querying
        VEP_API = "https://rest.ensembl.org/vep/human/hgvs"
    
        # Get a batch of 10 HGVS IDs from the dataset
        hgvs_ids = list(hgvs_df.iloc[i*10:i*10+10, 0])

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
            vep_data = response.json()  # Parse the API response
            print(len(vep_data))

            # Iterate through the response for each HGVS ID
            for k, result in enumerate(vep_data):
                for h in result['transcript_consequences']:
                    # Ensure necessary keys are present before processing
                    if 'protein_end' in h and 'amino_acids' in h and 'transcript_id' in h:
                        transcript_consequence = h

                        # Extract relevant data (amino acids, transcript ID, and protein position)
                        amino_acids = transcript_consequence.get('amino_acids', 'N/A')
                        transcript_id = transcript_consequence.get('transcript_id', 'N/A')
                        pos = transcript_consequence.get('protein_end', 'N/A')

                        # Append data to the respective lists
                        change.append(amino_acids)
                        ids.append(transcript_id)
                        poses.append(pos)
                        break  # Break after finding the first valid transcript consequence
        else:
            print('response issue')
            sys.exit()  # Exit the program if the API response fails
    
    # Return the lists of transcript IDs, positions, and amino acid changes
    return ids, poses, change

def fetch_protein_sequence_pdb(mrna):
    # Function to fetch the protein sequence from UniProt using a given mRNA
    search_url = f"https://rest.uniprot.org/uniprotkb/search?query={mrna}&format=txt"
    response = requests.get(search_url)

    # Check if the request was successful
    if response.status_code == 200:
        text = response.text
        seq_start = text.find("SQ")  # Find the start of the sequence in the response
        sequence_lines = text[seq_start:].splitlines()[1:]  # Skip the 'SQ' line
        sequence = ''.join(line.strip() for line in sequence_lines if line and not line.startswith('//'))
        sequence = sequence.replace(' ', '')  # Remove spaces from the sequence
        return sequence

def get_seq(mrna):
    # Function to fetch a list of protein sequences for a given mRNA using NCBI Entrez and UniProt APIs
    succes = False  # Track if the request is successful
    while not succes:
        # Set search parameters to query NCBI Entrez API for protein sequences
        params_0 = {
            'term': mrna,
            'retmode': 'json'
        }
        response = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein', params_0)

        # If the response is not successful, retry the request
        if response.status_code != 200:
            print('issue with protein seq response')
            continue

        # Parse the response to get the list of protein IDs
        prot_list = response.json()['esearchresult']['idlist']
        seq_list = []  # List to store the sequences
        succes = True  # Mark the request as successful

    # If no protein IDs are found, try fetching the sequence from UniProt
    if len(prot_list) == 0:
        seq_list.append(fetch_protein_sequence_pdb(mrna))
    else:
        # Loop through each protein ID and fetch the sequence using the NCBI Entrez efetch API
        for prot in prot_list:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
            params_seq = {
                'id': prot,
                'retmode': 'fasta',
                'rettype': 'fasta'
            }

            succes = False
            while not succes:
                response = requests.get(url, params_seq)
                if response.status_code == 200:
                    succes = True  # Mark the request as successful
            
            fasta = response.text
            seq = ''.join(fasta.split('\n')[1:]).strip()  # Extract and clean the FASTA sequence
            seq_list.append(seq)
    
    # Return the list of protein sequences
    return seq_list

def create_rows(dataframe):
    # Function to process data and create new rows in the dataframe with sequence and mutation information
    mismatch_count = 0  # Counter for mismatched sequences
    shorter_seq_count = 0  # Counter for sequences shorter than expected
    mrna, poses, amino_acids = get_ids()  # Retrieve the mRNA, positions, and amino acid changes
    print(len(mrna))
    print(len(poses))
    print(len(amino_acids))
    
    # Loop through each mRNA and process its associated sequences
    for num, i in enumerate(mrna):
        print(num)

        if num % 10 == 0:
            # Print progress and summary every 10 iterations
            print(f"Current number added is {num}")
            print(f"Mismatches between seq and clinvar data: {mismatch_count}")
            print(f"Sequences that were shorter than expected: {shorter_seq_count}")
        
        # Get the list of protein sequences for the current mRNA
        seq_list = get_seq(i)
        print(len(seq_list))

        # Loop through each protein sequence
        for seq in seq_list:
            if seq_list == 0:
                print('Missing sequence or error with request')
                continue

            length = len(seq)
            if len(seq) > 0:
                # Calculate the ratio of each amino acid in the sequence
                aa_ratios = []
                for aa in ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']:
                    aa_ratios.append(seq.count(aa) / length)

            # Get the reference amino acid, mutated amino acid, and position from the data
            id = hgvs_df.iloc[num, 0]
            aa_pos = int(poses[num])
            ref_aa = amino_acids[num][0]
            mut_aa = amino_acids[num][-1]

            # If the position is invalid (0), skip this sequence
            if aa_pos == 0:
                print('Not a mismatch mutation')
                continue

            # If the sequence length is shorter than the position, skip it
            if length < aa_pos:
                shorter_seq_count += 1
                continue

            # If the reference amino acid doesn't match the sequence, skip it
            if ref_aa != seq[aa_pos - 1]:
                continue
            print('Found a matching seq')

            # Break after finding the matching sequence
            break
        
        print('adding neighbors')
        neighbors = []
        for x in range(-5, 5):
            if aa_pos + x < length:
                neighbors.append(seq[aa_pos + x])
            else:
                neighbors.append(0)
        
        # Create a new row with mutation and sequence features
        mut_features = pd.Series(
            [id, ref_aa, mut_aa, *neighbors, *aa_ratios, len(seq)],
            index=["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", "neighbor_minus_2", "neighbor_minus_1", 
                   "self", "neighbor_plus_1", "neighbor_plus_2", "neighbor_plus_3", "neighbor_plus_4", "neighbor_plus_5", 
                   'R_ratio', 'H_ratio', 'K_ratio', 'D_ratio', 'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 
                   'G_ratio', 'P_ratio', 'A_ratio, 'V_ratio', 'I_ratio', 'L_ratio', 'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 'length' ])

        # concatenate the new row to the rest of the df
        dataframe = pd.concat([dataframe, mut_features.to_frame().T]) 
        
         
    return dataframe
        

def main():
    
    # Initializing our dataframe
    dataframe = pd.DataFrame(columns= ["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", "neighbor_minus_2", "neighbor_minus_1", "self", "neighbor_plus_1","neighbor_plus_2", "neighbor_plus_3", "neighbor_plus_4", "neighbor_plus_5", 'R_ratio','H_ratio','K_ratio', 'D_ratio', 'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 'G_ratio', 'P_ratio', 'A_ratio', 'V_ratio', 'I_ratio', 'L_ratio', 'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 'length' ])
    # Adding the row to it
    dataframe = create_rows(dataframe)
    output = 'test.tsv'
            
            
    # Saving it to the path specified above       
    dataframe.to_csv(output, sep='\t', index=False)
        

if __name__ == "__main__":
    main()
