import streamlit as st
import sys
import os
import tempfile
from io import StringIO
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
from parsers import fasta_parser, pdb_parser
from esmfold_api_tool import esmfold
from show_pdb import display_protein

def main():
    st.title("Predicción y visualización de proteínas con ESMFold2")

    # Carga del archivo FASTA
    fasta_file = st.file_uploader("Sube tu fichero FASTA (máx 400 aa)", type=["fasta"])
    if fasta_file is not None:
        if st.button("Ejecutar", use_container_width=True):
            with st.status("Ejecutando ESMFold2..."):
                content = fasta_file.read().decode("utf-8")
                st.code(f'Parseando fichero FASTA: {fasta_file.name}')

                fasta_content = StringIO(content)
                fasta_protein_dict = fasta_parser.parse_fasta_file(fasta_content=fasta_content)

                st.code(f'Secuencias incluidas: {list(fasta_protein_dict.keys())}')
                
                pdb_contents = {}
                tmp_files = {}
                for protein_name in fasta_protein_dict.keys():
                    st.code(f'Ejecutando ESMFold2 en {protein_name} via API')
                    
                    tmpfile = tempfile.NamedTemporaryFile(delete=False)
                    try:
                        esmfold.fold_sequence(fasta_seq=fasta_protein_dict[protein_name], output_pdb_filepath=tmpfile.name)
                        tmpfile.seek(0)
                        pdb_data = tmpfile.read().decode("utf-8")
                        pdb_contents[protein_name] = pdb_data
                        tmp_files[protein_name] = tmpfile.name
                    finally:
                        tmpfile.close()
                
                st.code(f'Todas las secuencias se han plegado exitosamente.')
                
                st.session_state.pdb_contents = pdb_contents
                st.session_state.tmp_files = tmp_files
                st.session_state.fasta_protein_dict = fasta_protein_dict

    # Visualización y descarga de la proteína
    if 'pdb_contents' in st.session_state:
        protein_selected = st.selectbox("Selecciona qué proteína quieres visualizar", options=list(st.session_state.pdb_contents.keys()))
        if protein_selected:
            pdb_data = st.session_state.pdb_contents[protein_selected]
            st.write(f"### Estructura 3D de {protein_selected}")
            display_protein.render_protein_3d(pdb_data)

            # Botón para descargar el fichero PDB
            st.download_button(
                label="Descargar fichero PDB",
                data=pdb_data,
                file_name=f"{protein_selected}.pdb",
                mime="chemical/x-pdb"
            )

            # Eliminar el fichero temporal después de la descarga
            temp_file_path = st.session_state.tmp_files[protein_selected]
            if os.path.exists(temp_file_path):
                os.remove(temp_file_path)

if __name__ == "__main__":
    main()