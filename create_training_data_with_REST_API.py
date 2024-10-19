import pandas as pd 
import requests
import xmltodict
import argparse
import xml.etree.ElementTree as ET

# Function to get ClinVar IDs for benign and pathogenic variants based on search parameters
def get_ids(amount):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"

    # Parameters for pathogenic variants
    params_pat = {
        "db": "clinvar",
        'term': "clinvar_snp[Filter] AND clinsig_pathogenic[Properties] AND missense_variant[molecular consequence] NOT moi_mitochondrial[prop] AND single_gene[prop]",
        "retmode": "json",     # Get response in JSON format,
        "retmax" : amount      # Limit the number of records returned
    }

    # Parameters for benign variants
    params_ben = {
        "db": "clinvar",
        'term': "clinvar_snp[Filter] AND clinsig_benign[Properties] AND missense_variant[molecular consequence] NOT moi_mitochondrial[prop] AND single_gene[prop]",
        "retmode": "json",     # Get response in JSON format,
        "retmax" : amount
    }

    # Get benign variants
    response = requests.get(url, params_ben)
    if response.status_code == 200:
        data = response.json() # Parse the response as JSON
    else:
        print(f'Error {response.status_code}')
    
    ids_ben = data['esearchresult']['idlist'] # List of benign variant IDs

    # Get pathogenic variants
    response = requests.get(url, params_pat)
    if response.status_code == 200:
        data = response.json() # Parse the response as JSON
    else:
        print('Error')

    ids_pat = data['esearchresult']['idlist'] # List of pathogenic variant IDs
    
    # Return IDs (skipping first 23000 to get relevant records)
    return ids_ben[23000:], ids_pat[23000:]

# Function to get summary information and accession (VCV) for a variant by its ID
def get_vcv(variant_id):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar"
    
    # Parameters to request summary information
    params = {
        'db': 'clinvar',
        'id': variant_id,
        'retmode' : 'json'    
    }
    
    # Try to get the summary and accession for the variant
    try:
        summary = requests.get(url, params).json()
        return summary, summary['result'][variant_id]['accession']
    except:
        return 0, 0  # Return 0 if there's an error

# Function to fetch VCV data for a given variant
def fetch_vcv(vcv):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar"
    
    # Parameters for fetching data
    params = {
        'db': 'clinvar',
        'id': vcv,
        'rettype' : 'vcv'  # Requesting XML format
    }
    
    # Try to fetch protein information from VCV
    try:
        root = ET.fromstring(requests.get(url, params).text)
        protein_expression = root.find(".//ProteinExpression")
        if protein_expression:
            protein = protein_expression.get('sequenceAccessionVersion')
            change = protein_expression.get('change')
        return protein, change
    except:
        return 0, 0  # Return 0 if there's an error

# Function to fetch protein sequence from UniProt using the mRNA sequence
def fetch_protein_sequence_uniprot(mrna):
    search_url = f"https://rest.uniprot.org/uniprotkb/search?query={mrna}&format=txt"
    response = requests.get(search_url)

    if response.status_code == 200:
        text = response.text
        seq_start = text.find("SQ")
        # Extracting protein sequence from UniProt response
        sequence_lines = text[seq_start:].splitlines()[1:]  # Skip the 'SQ' line
        sequence = ''.join(line.strip() for line in sequence_lines if line and not line.startswith('//'))
        sequence = sequence.replace(' ', '')
        return sequence

# Function to get protein sequence in FASTA format using the accession number
def get_seq(prot):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
    
    # Parameters for fetching the protein sequence
    params = {
        'id': prot,
        'retmode': 'fasta',  
        'rettype': 'fasta'
    }
    
    # Try to fetch the sequence from Entrez API
    try:
        response = requests.get(url, params)
        fasta = response.text
        seq = ''.join(fasta.split('\n')[1:]).strip()  # Get sequence without the FASTA header
        return seq
    except:
        return 0  # Return 0 if there's an error

# Function to create rows of data for each variant in the dataframe
def create_rows(dataframe, ids, label):
    mismatch_count = 0
    shorter_seq_count = 0
    
    # Loop through variant IDs
    for num, i in enumerate(ids):
        try:
            # Print progress every 20 variants
            if num % 20 == 0:
                print(f"Current number of {label} added is {num}")
                print(f"Mismatches between seq and clinvar data: {mismatch_count}")
                print(f"Sequences that were shorter than expected: {shorter_seq_count}")
            
            # Get VCV information for the variant
            summary, vcv = get_vcv(i)
            if summary == 0:
                continue
            
            # Fetch protein and mutation information
            protein, protein_change = fetch_vcv(vcv)
            if protein == 0:
                continue
            
            # Fetch protein sequence from Entrez
            seq = get_seq(protein)
            if seq == 0:
                print('Missing sequence or error with request')
                continue

            # Extract and process mutation information
            id = summary['result'][i]['variation_set'][0]['canonical_spdi']
            length = len(seq)
            aa_ratios = []
            for aa in ['R','H','K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']:
                aa_ratios.append(seq.count(aa) / length)  # Calculate amino acid ratios
            
            # Modify ID format to match the course data format
            index = id.find(':')
            id = id[:index + 1] + 'g.' + id[index + 1:]
            index = id.rfind(':')
            id = id[:index] + '>' + id[index + 1:]
            index = id.rfind(':')
            id = id[:index] + id[index + 1:]

            # Extract mutation details
            protein_change = protein_change.split('.')[-1]
            string = protein_change.strip()
            try:
                aa_pos = int(string[3:-3])  # Get mutation position
            except:
                aa_pos = 0
            ref_aa = string[0:3]  # Get reference amino acid
            mut_aa = string[-3:]  # Get mutated amino acid
            aa_dict = {
                'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Sec': 'U',
                'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L',
                'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T',
                'Trp': 'W', 'Tyr': 'Y', 'Val': 'V'
            }
            try:
                ref_aa = aa_dict[ref_aa]
                mut_aa = aa_dict[mut_aa]
            except: 
                continue

            if aa_pos == 0:
                print('Not a mismatch mutation')
                continue

            # Check if the sequence is long enough
            if length < aa_pos:
                shorter_seq_count += 1
                continue

            # Validate if the reference amino acid matches the sequence
            if ref_aa != seq[aa_pos - 1]:
                print(f'Mismatch between sequence and protein change indication')
                mismatch_count += 1
                continue

            # Get neighboring amino acids for analysis
            neighbors = []
            for j in range(-5, 5):
                if aa_pos + j < length:
                    neighbors.append(seq[aa_pos + j])
                else:
                    neighbors.append(0)

            # Create a row with all relevant mutation features and append to the dataframe
            mut_features = pd.Series([id, ref_aa, mut_aa, *neighbors, *aa_ratios, length, label], 
                                     index=["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", 
                                            "neighbor_minus_2", "neighbor_minus_1", "self", "neighbor_plus_1", 
                                            "neighbor_plus_2", "neighbor_plus_3", "neighbor_plus_4", 
                                            "neighbor_plus_5", 'R_ratio', 'H_ratio', 'K_ratio', 'D_ratio', 
                                            'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 
                                            'G_ratio', 'P_ratio', 'A_ratio', 'V_ratio', 'I_ratio', 'L_ratio', 
                                            'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 'length', 'label'])
            
            dataframe = pd.concat([dataframe, mut_features.to_frame().T])  # Append row to dataframe
        except:
            continue
         
    return dataframe  # Return updated dataframe

# Main function to handle data fetching and saving to a TSV file
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='out_path')  # Output file path argument
    output = parser.parse_args().out_path
    
    # Create an empty dataframe with the appropriate columns
    empty = pd.DataFrame(columns=["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", "neighbor_minus_2", 
                                  "neighbor_minus_1", "self", "neighbor_plus_1", "neighbor_plus_2", "neighbor_plus_3", 
                                  "neighbor_plus_4", "neighbor_plus_5", 'R_ratio', 'H_ratio', 'K_ratio', 'D_ratio', 
                                  'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 'G_ratio', 'P_ratio', 
                                  'A_ratio', 'V_ratio', 'I_ratio', 'L_ratio', 'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 
                                  'length', 'label"])
    dataframe = empty.copy()
    
    # Get lists of benign and pathogenic variant IDs
    ids_ben, ids_pat = get_ids(25000)
    
    # Process and add rows to the dataframe in batches of 500
    for i in range(len(ids_ben)):
        if (i + 1) % 500 == 0:
            if i > 500:
                dataframe = create_rows(empty, ids_ben[i-500:i], "Benign")
                dataframe = create_rows(dataframe, ids_pat[i-500:i], "Pathogenic")
            else:
                dataframe = create_rows(empty, ids_ben[0:i], "Benign")
                dataframe = create_rows(dataframe, ids_pat[0:i], "Pathogenic")
            
            # Save the dataframe to a TSV file
            dataframe.to_csv(output, sep='\t', index=False, mode='a', header=False)

# Run the main function if this script is executed
if __name__ == "__main__":
    main()
