import requests

RSCB_BASE_URL = "https://files.rcsb.org/download/"

def pdb_content(pdb_id:str):
    pdb_url = f"{RSCB_BASE_URL}{pdb_id}.pdb"
    response = requests.get(pdb_url)
    return response.content

def download_pdb_file(pdb_id, file_format = "pdb", file_path = None):
    pdb_url = f"{RSCB_BASE_URL}{pdb_id}.{file_format}"
    response = requests.get(pdb_url)
    if response.status_code == 200:
        if not file_path:
            file_path = f"{pdb_id}.{file_format}"
        with open(file_path, 'wb') as file:
            file.write(response.content)
        print(f"Archivo {file_path} descargado con éxito.")
    else:
        print(f"Fallo al descargar el archivo. Código de estado: {response.status_code}")
