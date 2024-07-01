import requests
import logging
import pandas as pd
import re

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

UNIPROT_BASE_URL = 'https://rest.uniprot.org/uniprotkb/'
UNIPROT_VARIANTS_URL = 'https://www.ebi.ac.uk/proteins/api/variation/'
UNIPROT_QUERY_URL = 'https://rest.uniprot.org/uniprotkb/search?query=accession:'

def modify_fasta_sequence(fasta_sequence: str, alt_sequence: str, begin: int, end: int) -> str:
    # Parse the FASTA sequence to get only the amino acid sequence
    lines = fasta_sequence.split('\n')
    amino_acid_sequence = ''.join(line for line in lines if not line.startswith('>'))
    name = f'>MUT({alt_sequence}-{str(begin)})'

    # Modificar la secuencia
    nueva_cadena = amino_acid_sequence[:begin-1] + alt_sequence + amino_acid_sequence[end:]
    
    return name + '\n' + '\n'.join([nueva_cadena[i:i+60] for i in range(0, len(nueva_cadena), 60)])

def get_fasta_from_uniprot_id(uniprot_id: str) -> str:
    url = f'{UNIPROT_BASE_URL}{uniprot_id}.fasta'
    fasta_seq = requests.get(url)
    if fasta_seq.status_code == 200:
        return fasta_seq.text
    else:
        return f'Error: Unable to fetch data. Status code: {fasta_seq.status_code}'

def get_alternative_fasta_from_uniprot_id(uniprot_id: str) -> list:
    url = f'{UNIPROT_VARIANTS_URL}{uniprot_id}?format=json'
    response = requests.get(url)
    results = []
    if response.status_code == 200:
        for variant in response.json()['features']:
            if variant['type'] == 'VARIANT':
                try:
                    alt_sequence = variant['alternativeSequence']
                    begin = int(variant['begin'])
                    end = int(variant['end'])
                    modified_fasta = modify_fasta_sequence(get_fasta_from_uniprot_id(uniprot_id), alt_sequence, begin, end)
                    result = {
                        'id': variant['xrefs'][0]['id'] if 'xrefs' in variant and len(variant['xrefs']) > 0 else None,
                        'prediction_type': variant['predictions'][0]['predAlgorithmNameType'] if 'predictions' in variant and len(variant['predictions']) > 0 else None,
                        'prediction_score': variant['predictions'][0]['score'] if 'predictions' in variant and len(variant['predictions']) > 0 else None,
                        'modified_fasta': modified_fasta
                    }
                    if result['prediction_score'] is not None and result['prediction_score'] > 0.9:
                        results.append(result)
                        if len(results) > 30: break
                except KeyError as e:
                    logging.warning(f"KeyError: {e} in variant {variant}")
        return results
    else:
        return f'Error: Unable to fetch data. Status code: {response.status_code}'
    
def query_variant_uniprot(uniprot_id: str, mutation: str):
    try:
        wt_aa = mutation[0]
        mut_aa = mutation[-1]
        pos = re.search(r'\d+', mutation).group()
        results = {}
        url = UNIPROT_VARIANTS_URL + f'{uniprot_id}?sourcetype=uniprot&wildtype={wt_aa}&alternativesequence={mut_aa}&location={pos}'
        response = requests.get(url)
        
        if response.status_code == 200:
            response_json = response.json()
            if 'sequence' in response_json:
                wt_seq = response_json['sequence']
                wt_list = list(wt_seq)
                wt_list[int(pos) - 1] = mut_aa
                mut_seq = ''.join(wt_list)
                
                results['response'] = response_json
                results['fasta'] = {
                    'wild_type': f">{uniprot_id}\n{wt_seq}",
                    'mutated': f">{uniprot_id}-{mutation}\n{mut_seq}"
                }
                return results
            else:
                return f'Error: La respuesta no contiene la secuencia.'
        else:
            return f'Error: Unable to fetch data. Status code: {response.status_code}'
    
    except Exception as e:
        return f'Error: {str(e)}'