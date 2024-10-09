import pandas as pd 
import requests
import xmltodict
import argparse
import xml.etree.ElementTree as ET

def get_ids(amount):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"

    params_pat = {
        "db": "clinvar",
        'term': "clinvar_snp[Filter] AND clinsig_pathogenic[Properties] AND missense_variant[molecular consequence] NOT moi_mitochondrial[prop] AND single_gene[prop]",
        "retmode": "json",     # Get response in JSON format,
        "retmax" : amount
        
    }

    params_ben = {
        "db": "clinvar",
        'term': "clinvar_snp[Filter] AND clinsig_benign[Properties] AND missense_variant[molecular consequence] NOT moi_mitochondrial[prop] AND single_gene[prop]",
        "retmode": "json",     # Get response in JSON format,
        "retmax" : amount
        
    }

    response = requests.get(url, params_ben)
    if response.status_code == 200:
        data = response.json() # Parse the response as JSON
        
    else:
        print(f'Error{response.status_code}')

    

    ids_ben = data['esearchresult']['idlist']
    #ids_ben = ",".join(ids_ben)



    response = requests.get(url, params_pat)
    if response.status_code == 200:
        data = response.json() # Parse the response as JSON
    else:
        print('Error')

    ids_pat = data['esearchresult']['idlist']
    
    return ids_ben[23000:], ids_pat[23000:]




def get_vcv(variant_id):
    # Define the URL for the Entrez efetch API
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar"
    
    # Set up parameters for the request
    params = {
        'db': 'clinvar',
        'id': variant_id,
        'retmode' : 'json'    
    }
    try:
        summary = requests.get(url, params).json()
        return summary, summary['result'][variant_id]['accession']
    except:
        return 0, 0
    

# Fetch the variant


def fetch_vcv(vcv):
    # Define the URL for the Entrez efetch API
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar"
    
    # Set up parameters for the request
    params = {
        'db': 'clinvar',
        'id': vcv,
        'rettype' : 'vcv'    # Requesting XML format
    }
    try:
        root = ET.fromstring(requests.get(url, params).text)
        protein_expression = root.find(".//ProteinExpression")
        if protein_expression:
            protein = protein_expression.get('sequenceAccessionVersion')
            change = protein_expression.get('change')
        return protein, change
    except:
        return 0,0




def fetch_protein_sequence_uniprot(mrna):
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
        

def get_seq(prot):
# Define the URL for the Entrez efetch API
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein"
    
    # Set up parameters for the request
    params = {
        
        'id': prot,
        'retmode': 'fasta',  
        'rettype': 'fasta'
    }
    try:
        response = requests.get(url, params)
    
        fasta = response.text
        seq = ''.join(fasta.split('\n')[1:]).strip()
            
        
        return seq
    except:
        return 0




def create_rows(dataframe, ids, label):
    mismatch_count = 0
    shorter_seq_count = 0
    for num, i in enumerate(ids):
        try:
            
            
            if num%20 == 0:
                print(f"Current number of {label} added is {num}")
                print(f"Mismatches between seq and clinvar data : {mismatch_count}")
                print(f"Sequences that were shorter than expected: {shorter_seq_count}")
            
            summary, vcv = get_vcv(i)
            if summary == 0:
                continue
            protein, protein_change = fetch_vcv(vcv)
            if protein == 0:
                continue
            
            seq = get_seq(protein)
            
            if seq == 0:
                print('Missing sequence or error with request')
                continue
            
            
            id = summary['result'][i]['variation_set'][0]['canonical_spdi']
            length = len(seq)
            aa_ratios = []
            for aa in ['R','H','K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']:
                aa_ratios.append(seq.count(aa)/length)
                
            
            
            
            # Transforming the index used in the tsvs for the course
            index = id.find(':')
            id = id[:index + 1] + 'g.' + id[index + 1:]
            index = id.rfind(':')
            id = id[:index] + '>' + id[index +1:]
            index = id.rfind(':')
            id = id[:index] + id[index + 1:]
            
            
            #getting the reference aa, mutated aa and position of the longest protein associated with the mutation
            
            
            protein_change = protein_change.split('.')[-1]
            string = protein_change.strip()
            try:
                aa_pos = int(string[3:-3])
            except:
                aa_pos = 0
            ref_aa = string[0:3]
            mut_aa = string[-3:]
            aa_dict = {
            'Ala': 'A',  # Alanine
            'Arg': 'R',  # Arginine
            'Asn': 'N',  # Asparagine
            'Asp': 'D',  # Aspartic acid
            'Cys': 'C',  # Cysteine
            'Sec': 'U',  # Selenocysteine
            'Gln': 'Q',  # Glutamine
            'Glu': 'E',  # Glutamic acid
            'Gly': 'G',  # Glycine
            'His': 'H',  # Histidine
            'Ile': 'I',  # Isoleucine
            'Leu': 'L',  # Leucine
            'Lys': 'K',  # Lysine
            'Met': 'M',  # Methionine
            'Phe': 'F',  # Phenylalanine
            'Pro': 'P',  # Proline
            'Ser': 'S',  # Serine
            'Thr': 'T',  # Threonine
            'Trp': 'W',  # Tryptophan
            'Tyr': 'Y',  # Tyrosine
            'Val': 'V'   # Valine
            }
            try:
                ref_aa = aa_dict[ref_aa]
                mut_aa = aa_dict[mut_aa]
            except: 
                continue
                
            
            if aa_pos == 0:
                print('Not a mismatch mutation')
                continue
            
            if length < aa_pos :
                continue
            
                for i in seq2:
                    if len(i) > aa_pos + 1:
                        seq = i
                        
                        break
                
                shorter_seq_count +=1 
                continue
            
            if ref_aa != seq[aa_pos-1]:
                print(f'Mismatch between sequence and protein change indication')
                
                
                
                
                mismatch_count +=1
                continue
            
            neighbors = []
            for i in range(-5,5):
                if aa_pos + i < length:
                    neighbors.append(seq[aa_pos + i])
                else:
                    neighbors.append(0)
                    
                        
            mut_features = pd.Series([id, ref_aa, mut_aa, *neighbors, *aa_ratios, length, label], index = ["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", "neighbor_minus_2", "neighbor_minus_1", "self", "neighbor_plus_1","neighbor_plus_2", "neighbor_plus_3", "neighbor_plus_4", "neighbor_plus_5", 'R_ratio','H_ratio','K_ratio', 'D_ratio', 'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 'G_ratio', 'P_ratio', 'A_ratio', 'V_ratio', 'I_ratio', 'L_ratio', 'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 'length', 'label' ])
            
            dataframe = pd.concat([dataframe, mut_features.to_frame().T])  
        except:
            continue
         
    return dataframe
        

# Example usage:
# clinvar_id = "12345"  # Replace with an actual ClinVar ID
# print(get_protein_sequence_from_clinvar(clinvar_id))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='out_path')
    output = parser.parse_args().out_path
    #creating the dataframe with the correct columns
    empty = pd.DataFrame(columns= ["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", "neighbor_minus_2", "neighbor_minus_1", "self", "neighbor_plus_1","neighbor_plus_2", "neighbor_plus_3", "neighbor_plus_4", "neighbor_plus_5", 'R_ratio','H_ratio','K_ratio', 'D_ratio', 'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 'G_ratio', 'P_ratio', 'A_ratio', 'V_ratio', 'I_ratio', 'L_ratio', 'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 'length', 'label' ])
    dataframe = pd.DataFrame(columns= ["id", "ref_aa", "mut_aa", "neighbor_minus_4", "neighbor_minus_3", "neighbor_minus_2", "neighbor_minus_1", "self", "neighbor_plus_1","neighbor_plus_2", "neighbor_plus_3", "neighbor_plus_4", "neighbor_plus_5", 'R_ratio','H_ratio','K_ratio', 'D_ratio', 'E_ratio', 'S_ratio', 'T_ratio', 'N_ratio', 'Q_ratio', 'C_ratio', 'G_ratio', 'P_ratio', 'A_ratio', 'V_ratio', 'I_ratio', 'L_ratio', 'M_ratio', 'F_ratio', 'Y_ratio', 'W_ratio', 'length', 'label' ])
    ids_ben, ids_pat = get_ids(25000)  
    
    for i in range(len(ids_ben)):
        
        if (i+1)%500 == 0:
            if i > 500:
                dataframe = create_rows(empty, ids_ben[i-500:i], "Benign")
                dataframe = create_rows(dataframe, ids_pat[i-500:i], "Pathogenic")
            else:
                dataframe = create_rows(empty, ids_ben[0:i], "Benign")
                dataframe = create_rows(dataframe, ids_pat[0:i], "Pathogenic")
            
            
            
            
            
            dataframe.to_csv(output, sep='\t', index=False, mode='a', header= False)
        

if __name__ == "__main__":
    main()