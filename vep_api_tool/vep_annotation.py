import requests
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

VEP_BASE_URL = 'https://rest.ensembl.org/variation/human/'

def get_variant_annotation(rsid:str) -> dict:
  url = f'{VEP_BASE_URL}{rsid}/?content-type=application/json;phenotypes=1'
  result_json = {}
  response = requests.get(url)
  if response.status_code == 200:
      annotation_data = response.json()
      result_json['clinical_significance'] = annotation_data['clinical_significance'][0]
      result_json['trait'] = annotation_data['phenotypes'][0]['trait']
      result_json['gene'] = annotation_data['phenotypes'][0]['genes']
      result_json['variant'] = rsid
      return result_json
  else:
      logging.error(f"Error en la solicitud: {response.status_code}")