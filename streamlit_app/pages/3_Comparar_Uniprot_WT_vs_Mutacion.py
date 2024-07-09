import streamlit as st
import sys
import os
import logging
import tempfile
import zipfile
import io

# Setup logging
logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(message)s')

# Import custom modules
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from parsers import fasta_parser
from esmfold_api_tool import esmfold
from show_pdb import display_protein, display_quality_data
from uniprot_api_funcs.uniprot import query_variant_uniprot
from eval.eval_sequences import align_sequences, ramachandran_plot

def main():
    st.title("Herramienta de comparación de estructuras II")

    st.write(
        """
        - **Paso 1:** Introduce el UniProt ID de la primera proteína. Se abrirá un pop-up con su información de Uniprot. *Es recomendable expandir la barra lateral para una mejor visualización*
        - **Paso 2:** Introduce la mutación con el formato A123B. *Si el cambio es Q>P y la posición es 2, hay que introducir Q2P*
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
                with st.spinner(f"Plegando la secuencia {uniprot_query_code}-WT ..."):
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as tmpfile_wt:
                        esmfold.fold_sequence(fasta_seq=st.session_state.fasta_wt, output_pdb_filepath=tmpfile_wt.name)
                        tmpfile_wt.seek(0)
                        wt_pdb_data = tmpfile_wt.read().decode("utf-8")
                        st.session_state.wt_pdb_filepath = tmpfile_wt.name
                
                with st.spinner(f"Plegando la secuencia {uniprot_query_code}-{mutation} ..."):
                    with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as tmpfile_mutated:
                        esmfold.fold_sequence(fasta_seq=st.session_state.fasta_mutated, output_pdb_filepath=tmpfile_mutated.name)
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
                    alignment_score, alignment_img_buf = align_sequences(st.session_state.fasta_wt, st.session_state.fasta_mutated)
                    ramachandran_img_both_buf = ramachandran_plot(st.session_state.wt_pdb_filepath, st.session_state.mutated_pdb_filepath)
                    
                    colu2.metric("Alignment Score", f"{alignment_score}")

                    with colu1:
                        st.image(alignment_img_buf, caption="Alineamiento de secuencias")
                    with colu2:
                        st.image(ramachandran_img_both_buf, caption="Gráfico de Ramachandran")

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
                st.error(f"Error al plegar las secuencias: {e}")
            finally:
                if 'wt_pdb_filepath' in st.session_state:
                    os.unlink(st.session_state.wt_pdb_filepath)
                if 'mutated_pdb_filepath' in st.session_state:
                    os.unlink(st.session_state.mutated_pdb_filepath)

if __name__ == "__main__":
    main()
