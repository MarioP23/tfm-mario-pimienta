import requests
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

ESMFOLD_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

def fold_sequence(fasta_seq:str, output_pdb_filepath:str) -> dict:
    results = {}
    if len(fasta_seq) >= 400:
        logging.error("Error: la longitud de data debe ser menor a 400 caracteres.")
        exit
    else:
        try:
            response = requests.post(ESMFOLD_URL, data=fasta_seq, verify=False)
            response.raise_for_status()  # Esto lanzará una excepción si la solicitud no fue exitosa
            results['status_code'] = response.status_code
            with open(output_pdb_filepath, "w") as pdb_file:
                pdb_file.write(response.text)
            results['pdb_file_path'] = output_pdb_filepath
            return results
        except requests.exceptions.RequestException as e:
            logging.error(f"Ocurrió un error al hacer la solicitud: {e}")