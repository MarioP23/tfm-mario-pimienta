import streamlit as st
import sys
import os
from io import StringIO
import tempfile
import requests
import re

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from parsers import fasta_parser
from esmfold_api_tool import esmfold
from show_pdb import display_protein, display_quality_data
from uniprot_api_funcs.uniprot import query_variant_uniprot

def main():
    st.title("Herramienta de comparación de estructuras II")

    st.write(
        """
        - **Paso 1:** Introduce el UniProt ID de la primera proteína. Se abrirá un pop-up con su información de Uniprot. *Es recomendable expandir la barra lateral para una mejor visualización*
        - **Paso 2:** Introduce la mutación con el formato A123B.
        - **Paso 3:** Compara las dos proteínas para observar las diferencias causadas por la mutación.
        """
    )

    # Entrada del código UniProt
    uniprot_query_code = st.text_input("Introduce el código UniProt", placeholder="P09936")
    
    if st.button("Buscar UniProt ID", use_container_width=True):
        st.session_state.uniprot_query_code = uniprot_query_code

    if 'uniprot_query_code' in st.session_state:
        with st.sidebar:
            url = f"https://www.uniprot.org/uniprotkb/{st.session_state.uniprot_query_code}/variant-viewer"
            st.components.v1.iframe(url, height=600, scrolling=True)
        
        mutation = st.text_input("Introduce la mutación a buscar", placeholder="Q2P")

        if st.button("Buscar mutación", use_container_width=True):
            with st.spinner(f"Buscando secuencias {st.session_state.uniprot_query_code}-WT y {mutation} en UniProt..."):
                try:
                    query_uniprot = query_variant_uniprot(uniprot_id=st.session_state.uniprot_query_code, mutation=mutation)
                    if 'Error' in query_uniprot:
                        st.error(query_uniprot)
                    else:
                        fasta_wt = query_uniprot['fasta']['wild_type']
                        fasta_mutated = query_uniprot['fasta']['mutated']
                        st.success(f"Secuencias FASTA obtenidas para {st.session_state.uniprot_query_code}-WT y {st.session_state.uniprot_query_code}-{mutation}")
                        
                        st.session_state.fasta_wt = fasta_wt
                        st.session_state.fasta_mutated = fasta_mutated

                except Exception as e:
                    st.error(f"Error al obtener o procesar las mutaciones de UniProt: {e}")

    # Plegamiento y visualización de proteínas
    if 'fasta_wt' in st.session_state and 'fasta_mutated' in st.session_state:
        if st.button("Plegar proteínas", use_container_width=True):
            try:
                with st.spinner("Plegando la secuencia Wild Type..."):
                    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
                        esmfold.fold_sequence(fasta_seq=st.session_state.fasta_wt, output_pdb_filepath=tmpfile.name)
                        tmpfile.seek(0)
                        wt_pdb_data = tmpfile.read().decode("utf-8")
                        st.session_state.wt_pdb_data = wt_pdb_data
                
                with st.spinner("Plegando la secuencia Mutated..."):
                    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
                        esmfold.fold_sequence(fasta_seq=st.session_state.fasta_mutated, output_pdb_filepath=tmpfile.name)
                        tmpfile.seek(0)
                        mutated_pdb_data = tmpfile.read().decode("utf-8")
                        st.session_state.mutated_pdb_data = mutated_pdb_data
                
                # Comparación de las proteínas
                display_protein.show_protein_grid(pdb_data_1=st.session_state.wt_pdb_data, pdb_data_2=st.session_state.mutated_pdb_data)
                try:
                    rmsd_value = display_quality_data.calculate_rmsd_from_strings(
                        pdb_str_1=st.session_state.wt_pdb_data, 
                        pdb_str_2=st.session_state.mutated_pdb_data
                    )
                    st.metric(f"RMSD", f"{rmsd_value:.3f} Å")
                except Exception as e:
                    st.error(f"Error al calcular RMSD: {e}")
                
                # Descarga de ficheros PDB
                st.download_button(
                    label="Descargar Wild Type PDB",
                    data=st.session_state.wt_pdb_data,
                    file_name="wild_type.pdb",
                    mime="chemical/x-pdb"
                )
                
                st.download_button(
                    label="Descargar Mutated PDB",
                    data=st.session_state.mutated_pdb_data,
                    file_name="mutated.pdb",
                    mime="chemical/x-pdb"
                )
                
                # Crear fichero multiFASTA
                multifasta = f"{st.session_state.fasta_wt}\n{st.session_state.fasta_mutated}"
                st.download_button(
                    label="Descargar multiFASTA",
                    data=multifasta,
                    file_name="multifasta.fasta",
                    mime="text/plain"
                )

                # Crear fichero PDB combinado
                combined_pdb = st.session_state.wt_pdb_data + "\n" + st.session_state.mutated_pdb_data
                st.download_button(
                    label="Descargar PDB combinado",
                    data=combined_pdb,
                    file_name="combined.pdb",
                    mime="chemical/x-pdb"
                )
                
            except Exception as e:
                st.error(f"Error al plegar las secuencias: {e}")

if __name__ == "__main__":
    main()
