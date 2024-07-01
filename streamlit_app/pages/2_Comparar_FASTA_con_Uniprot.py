import streamlit as st
import sys
import os
from io import StringIO
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from parsers import fasta_parser
from esmfold_api_tool import esmfold
from show_pdb import display_protein, display_quality_data
import tempfile
import pandas as pd
from uniprot_api_funcs.uniprot import get_alternative_fasta_from_uniprot_id

def main():
    st.title("Herramienta de comparación de estructuras I")

    st.write(
        """
        - **Paso 1:** Sube tu fichero FASTA local.
        - **Paso 2:** Introduce el accesion ID de la proteína de UniProt.
        - **Paso 3:** Se expandirá una tabla con 30 variantes patogénicas relacionadas y podrás seleccionar una para comparar. *Puede tardar varios minutos*
        """
    )

    # Carga del archivo FASTA
    fasta_file = st.file_uploader("Sube tu fichero FASTA", type=["fasta", "fa"])
    if fasta_file:
        if st.button("Plegar FASTA"):
            with st.spinner("Plegando la secuencia del archivo FASTA..."):
                try:
                    fasta_content = fasta_file.read().decode("utf-8")
                    fasta_protein_dict = fasta_parser.parse_fasta_file(fasta_content=StringIO(fasta_content))

                    st.code(f'Secuencias incluidas: {list(fasta_protein_dict.keys())}')

                    pdb_contents = {}
                    tmp_files = {}
                    for protein_name in fasta_protein_dict.keys():
                        if len(fasta_protein_dict[protein_name]) > 400: 
                            st.error("El fichero FASTA contiene más de 400 aminoácidos.")
                            break
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

                    st.session_state.pdb_contents_uploaded = pdb_contents
                    st.session_state.tmp_files_uploaded = tmp_files
                    st.session_state.fasta_protein_dict_uploaded = fasta_protein_dict

                    # Visualizar la proteína plegada inmediatamente
                    protein_name = list(fasta_protein_dict.keys())[0]
                    pdb_data = st.session_state.pdb_contents_uploaded[protein_name]
                    st.session_state.pdb_data_2 = pdb_data
                    st.write(f"### Estructura 3D de {protein_name}")
                    display_protein.render_protein_3d(pdb_data)

                    # Botón para descargar el fichero PDB
                    st.download_button(
                        label="Descargar fichero PDB",
                        data=pdb_data,
                        file_name=f"{protein_name}.pdb",
                        mime="chemical/x-pdb"
                    )

                    # Eliminar el fichero temporal después de la descarga
                    temp_file_path = st.session_state.tmp_files_uploaded[protein_name]
                    if os.path.exists(temp_file_path):
                        os.remove(temp_file_path)

                except Exception as e:
                    st.error(f"Error al plegar la secuencia del archivo FASTA: {e}")

    # Entrada del código UniProt
    uniprot_query_code = st.text_input("Introduce el código UniProt")
    if st.button("Buscar", use_container_width=True):
        with st.spinner("Buscando las 30 primeras variantes desde UniProt..."):
            try:
                alternative_variants = get_alternative_fasta_from_uniprot_id(uniprot_query_code)
                if isinstance(alternative_variants, str):
                    st.error(alternative_variants)
                else:
                    df_variants = pd.DataFrame(alternative_variants).drop(columns=['modified_fasta'])
                    df_variants.columns = ['ID', 'Algorythm', 'Patogenicity score']
                    df_variants = df_variants.sort_values(by='Patogenicity score', ascending=False)
                    st.session_state.df_variants = df_variants
                    st.session_state.alternative_variants = alternative_variants
            except Exception as e:
                st.error(f"Error al obtener o procesar las variantes de UniProt: {e}")

    # Mostrar el dataframe si está en session_state
    if 'df_variants' in st.session_state:
        st.subheader("Variantes encontradas")
        st.dataframe(st.session_state.df_variants, use_container_width=True, hide_index=True)

    # Plegamiento de la variante seleccionada
    if 'df_variants' in st.session_state:
        st.subheader("Selecciona una variante para plegar y visualizar")
        #selected_variant = st.selectbox("Selecciona la variante", options=st.session_state.df_variants.index)
        selected_variant = st.selectbox("Selecciona la variante", options=st.session_state.df_variants.index, format_func=lambda x: st.session_state.df_variants.at[x, 'ID'])

        if st.button("Plegar Mutated FASTA"):
            with st.spinner("Plegando la secuencia seleccionada..."):
                try:
                    selected_variant_data = st.session_state.df_variants.loc[selected_variant]
                    modified_fasta = st.session_state.alternative_variants[selected_variant]['modified_fasta']
                    fasta_protein_dict = fasta_parser.parse_fasta_file(fasta_content=StringIO(modified_fasta))

                    st.code(f'Secuencias incluidas: {list(fasta_protein_dict.keys())}')

                    pdb_contents = {}
                    tmp_files = {}
                    for protein_name in fasta_protein_dict.keys():
                        st.code(f'Ejecutando ESMFold2 en {protein_name} via API')

                        tmpfile = tempfile.NamedTemporaryFile(delete=False)
                        try:
                            esmfold.fold_sequence(fasta_seq=fasta_protein_dict[protein_name], output_pdb_filepath=tmpfile.name)
                            st.code(f'Todas las secuencias se han plegado exitosamente.')
                            tmpfile.seek(0)
                            pdb_data = tmpfile.read().decode("utf-8")
                            pdb_contents[protein_name] = pdb_data
                            tmp_files[protein_name] = tmpfile.name
                        finally:
                            tmpfile.close()

                    st.session_state.pdb_contents = pdb_contents
                    st.session_state.tmp_files = tmp_files
                    st.session_state.fasta_protein_dict = fasta_protein_dict

                    # Visualizar la proteína plegada inmediatamente
                    protein_name = list(fasta_protein_dict.keys())[0]
                    pdb_data = st.session_state.pdb_contents[protein_name]
                    st.session_state.pdb_data_1 = pdb_data
                    st.write(f"### Estructura 3D de {protein_name}")
                    display_protein.render_protein_3d(pdb_data)

                    # Botón para descargar el fichero PDB
                    st.download_button(
                        label="Descargar fichero PDB",
                        data=pdb_data,
                        file_name=f"{protein_name}.pdb",
                        mime="chemical/x-pdb"
                    )

                    # Eliminar el fichero temporal después de la descarga
                    temp_file_path = st.session_state.tmp_files[protein_name]
                    if os.path.exists(temp_file_path):
                        os.remove(temp_file_path)

                except Exception as e:
                    st.error(f"Error al plegar la secuencia seleccionada: {e}")

    # Comparación de las proteínas
    if st.button("Comparar", use_container_width=True):
        if 'pdb_data_1' in st.session_state and 'pdb_data_2' in st.session_state:
            display_protein.show_protein_grid(pdb_data_1=st.session_state.pdb_data_1, pdb_data_2=st.session_state.pdb_data_2)
            try:
                st.metric(f"RMSD", f"{display_quality_data.calculate_rmsd_from_strings(pdb_str_1=st.session_state.pdb_data_1, pdb_str_2=st.session_state.pdb_data_2):.3f} Å")
            except Exception as e:
                st.error(f"Error al calcular RSMD: {e}")

        else:
            if 'pdb_data_1' not in st.session_state:
                st.error("Por favor, proporciona un código UniProt válido y selecciona una variante.")
            if 'pdb_data_2' not in st.session_state:
                st.error("Por favor, sube un archivo FASTA válido.")

if __name__ == "__main__":
    if 'pdb_data_1' not in st.session_state:
        st.session_state.pdb_data_1 = None
    if 'pdb_data_2' not in st.session_state:
        st.session_state.pdb_data_2 = None
    main()
