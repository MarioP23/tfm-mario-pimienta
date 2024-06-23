import streamlit as st

def main():
    st.set_page_config(
    page_title="Home",
    layout="wide"
    )

    st.title("Protein Structure Prediction and Visualization")
    
    st.header("Descripción del Proyecto")
    st.write("""
    Este proyecto tiene como objetivo proporcionar herramientas para la comparación y visualización de estructuras proteicas utilizando datos PDB. 
    Se utiliza la biblioteca Py3Dmol para la visualización en 3D y Streamlit para la creación de la interfaz de usuario.
    Además, se ofrece la capacidad de predecir estructuras proteicas a partir de secuencias FASTA utilizando la API ESMFold.
    """)

    st.header("Funcionalidades")
    st.subheader("Comparación de Estructuras PDB")
    st.write("""
    - **Entrada por código RSCB PDB**: Permite al usuario introducir un código PDB para descargar y visualizar una estructura.
    - **Carga de archivos PDB**: Permite al usuario cargar su propio archivo PDB para la visualización.
    - **Personalización de la Visualización**: Ofrece opciones para personalizar la visualización, incluyendo estilo, opacidad de la superficie, color de la proteína y color de fondo.
    """)

    st.subheader("Predicción y Visualización de Estructuras Proteicas")
    st.write("""
    - **Carga de archivos FASTA**: Permite al usuario cargar un archivo FASTA para la predicción de la estructura proteica.
    - **Integración con ESMFold**: Utiliza la API ESMFold para predecir la estructura proteica a partir de la secuencia FASTA.
    - **Visualización 3D**: Muestra la estructura 3D predicha de la proteína.
    - **Descarga de archivos PDB**: Permite descargar la estructura PDB predicha.
    """)

    st.header("Ejemplo de Uso")
    st.subheader("Comparación de Estructuras PDB")
    st.write("""
    1. Introduce un código RSCB PDB y presiona "Buscar por código RSCB PDB".
    2. Sube un archivo PDB desde tu máquina.
    3. Visualiza ambas proteínas
    """)

    st.subheader("Predicción y Visualización de Estructuras Proteicas")
    st.write("""
    1. Sube un archivo FASTA.
    2. Presiona "Ejecutar" para predecir la estructura proteica utilizando ESMFold.
    3. Selecciona la proteína que deseas visualizar y ajusta las opciones de visualización según tus preferencias.
    4. Descarga el archivo PDB si lo deseas.
    """)

    st.header("Bibliografía")
    st.write("""
    EJEMPLOS:
    
    1. Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H., ... & Bourne, P. E. (2000). The Protein Data Bank. Nucleic acids research, 28(1), 235-242.
    2. Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., ... & Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596(7873), 583-589.
    """)

if __name__ == "__main__":
    main()
