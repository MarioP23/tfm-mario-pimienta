import streamlit as st
import sys
import os
from io import StringIO
import tempfile
import pandas as pd
import zipfile
import io

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from parsers import fasta_parser
from esmfold_api_tool import esmfold
from show_pdb import display_protein, display_quality_data
from uniprot_api_funcs.uniprot import get_alternative_fasta_from_uniprot_id
from eval.eval_sequences import align_sequences, ramachandran_plot

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
        try:
            fasta_content = fasta_file.read().decode("utf-8")
            fasta_protein_dict = fasta_parser.parse_fasta_file(fasta_content=StringIO(fasta_content))
            st.success("Fichero FASTA cargado exitosamente.")
            st.session_state.fasta_protein_dict_uploaded = fasta_protein_dict
        except Exception as e:
            st.error(f"Error al cargar el archivo FASTA: {e}")

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
                    df_variants.columns = ['ID', 'Algorithm', 'Pathogenicity score']
                    df_variants = df_variants.sort_values(by='Pathogenicity score', ascending=False)
                    st.session_state.df_variants = df_variants
                    st.session_state.alternative_variants = alternative_variants
            except Exception as e:
                st.error(f"Error al obtener o procesar las variantes de UniProt: {e}")

    # Mostrar el dataframe si está en session_state
    if 'df_variants' in st.session_state:
        st.subheader("Variantes encontradas")
        st.dataframe(st.session_state.df_variants, use_container_width=True, hide_index=True)

    # Comparación de las proteínas
    if 'df_variants' in st.session_state and 'fasta_protein_dict_uploaded' in st.session_state:
        st.subheader("Selecciona una variante para comparar")
        selected_variant = st.selectbox("Selecciona la variante", options=st.session_state.df_variants.index, format_func=lambda x: st.session_state.df_variants.at[x, 'ID'])

        if st.button("Comparar", use_container_width=True):
            with st.spinner("Plegando y comparando las secuencias..."):
                try:
                    # Plegar el FASTA cargado
                    protein_name = list(st.session_state.fasta_protein_dict_uploaded.keys())[0]
                    fasta_seq = st.session_state.fasta_protein_dict_uploaded[protein_name]

                    with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as tmpfile_wt:
                        esmfold.fold_sequence(fasta_seq=fasta_seq, output_pdb_filepath=tmpfile_wt.name)
                        tmpfile_wt.seek(0)
                        wt_pdb_data = tmpfile_wt.read().decode("utf-8")
                        st.session_state.wt_pdb_filepath = tmpfile_wt.name

                    # Plegar la variante seleccionada
                    selected_variant_data = st.session_state.df_variants.loc[selected_variant]
                    modified_fasta = st.session_state.alternative_variants[selected_variant]['modified_fasta']
                    fasta_protein_dict = fasta_parser.parse_fasta_file(fasta_content=StringIO(modified_fasta))
                    protein_name_mutated = list(fasta_protein_dict.keys())[0]
                    mutated_fasta_seq = fasta_protein_dict[protein_name_mutated]

                    with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as tmpfile_mutated:
                        esmfold.fold_sequence(fasta_seq=mutated_fasta_seq, output_pdb_filepath=tmpfile_mutated.name)
                        tmpfile_mutated.seek(0)
                        mutated_pdb_data = tmpfile_mutated.read().decode("utf-8")
                        st.session_state.mutated_pdb_filepath = tmpfile_mutated.name

                    # Comparación de las proteínas
                    display_protein.show_protein_grid(pdb_data_1=wt_pdb_data, pdb_data_2=mutated_pdb_data)
                    
                    # Calculo del RMSD
                    try:
                        rmsd_value = display_quality_data.calculate_rmsd_from_strings(
                            pdb_str_1=wt_pdb_data, 
                            pdb_str_2=mutated_pdb_data
                        )
                        colu1, colu2 = st.columns(2)
                        colu1.metric("RMSD", f"{rmsd_value:.3f} Å")
                    except Exception as e:
                        st.error(f"Error al calcular RMSD: {e}")
                    
                    # Calculo de otros valores
                    try:
                        alignment_score, alignment_img_buf = align_sequences(fasta_seq, mutated_fasta_seq)
                        ramachandran_img_both_buf = ramachandran_plot(st.session_state.wt_pdb_filepath, st.session_state.mutated_pdb_filepath)
                        
                        colu2.metric("Alignment Score", f"{alignment_score}")

                        with colu1:
                            st.image(alignment_img_buf, caption="Alineamiento de secuencias")
                        with colu2:
                            st.image(ramachandran_img_both_buf, caption="Gráfico de Ramachandran - WT")

                    except Exception as e:
                        st.error(f"Error al calcular las métricas: {e}")

                    # Crear un archivo ZIP en memoria
                    buffer = io.BytesIO()
                    multifasta = f"{st.session_state.fasta_wt}\n{st.session_state.fasta_mutated}"
                    with zipfile.ZipFile(buffer, "w") as zip_file:
                        zip_file.writestr(f"{uniprot_query_code}-WildType.pdb", wt_pdb_data)
                        zip_file.writestr(f"{uniprot_query_code}-{mutation}.pdb", mutated_pdb_data)
                        zip_file.writestr(f"{uniprot_query_code}-WT_{mutation}.fasta", multifasta)

                    # Mover el puntero al inicio del archivo
                    buffer.seek(0)

                    # Botón para descargar el archivo ZIP
                    st.download_button(
                        label=f"Descargar ficheros PDB y FASTA de {uniprot_query_code}",
                        data=buffer,
                        file_name=f"{uniprot_query_code}_files.zip",
                        mime="application/zip",
                        use_container_width=True
                    )

                except Exception as e:
                    st.error(f"Error al plegar y comparar las secuencias: {e}")
                finally:
                    if 'wt_pdb_filepath' in st.session_state:
                        os.unlink(st.session_state.wt_pdb_filepath)
                    if 'mutated_pdb_filepath' in st.session_state:
                        os.unlink(st.session_state.mutated_pdb_filepath)

if __name__ == "__main__":
    if 'pdb_data_1' not in st.session_state:
        st.session_state.pdb_data_1 = None
    if 'pdb_data_2' not in st.session_state:
        st.session_state.pdb_data_2 = None
    main()
