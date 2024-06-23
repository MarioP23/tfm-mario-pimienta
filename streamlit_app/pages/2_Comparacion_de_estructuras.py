import streamlit as st
from data import download_pdb_file_rscb
from show_pdb.display_protein import show_protein_grid

def main():
    st.title("Herramienta de comparación de estructuras")

    # Entrada del código RSCB PDB
    pdb_query = st.text_input("Introduce el código RSCB PDB", placeholder="2PZD")
    if st.button("Buscar por código RSCB PDB", use_container_width=True):
        pdb_data_2_bytes = download_pdb_file_rscb.pdb_content(pdb_query)
        if pdb_data_2_bytes:
            st.session_state.pdb_data_2 = pdb_data_2_bytes.decode('utf-8')
            st.success(f"Proteína {pdb_query} cargada con éxito")
        else:
            st.error("No se pudo descargar los datos PDB. Verifica el código.")

    # Carga del archivo PDB
    pdb_file = st.file_uploader("Sube tu fichero PDB", type=["pdb"])
    if pdb_file:
        st.session_state.pdb_data = pdb_file.read().decode("utf-8")
        st.success("Archivo PDB cargado exitosamente.")

    # Comparación de las proteínas
    if st.button("Visualizar proteínas", use_container_width=True):
        if 'pdb_data' in st.session_state and 'pdb_data_2' in st.session_state:
            show_protein_grid(pdb_data_1=st.session_state.pdb_data, pdb_data_2=st.session_state.pdb_data_2)
        else:
            if 'pdb_data' not in st.session_state:
                st.error("Por favor, sube un archivo PDB.")
            if 'pdb_data_2' not in st.session_state:
                st.error("Por favor, proporciona un código RSCB PDB válido.")

if __name__ == "__main__":
    if 'pdb_data' not in st.session_state:
        st.session_state.pdb_data = None
    if 'pdb_data_2' not in st.session_state:
        st.session_state.pdb_data_2 = None
    main()
