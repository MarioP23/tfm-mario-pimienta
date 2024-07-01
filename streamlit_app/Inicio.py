import streamlit as st

def main():
    st.set_page_config(
    page_title="Home",
    layout="wide"
    )

    st.title("ANÁLISIS DE ESTRUCTURAS DE PROTEÍNAS CON ESMFold2 y UniProt")
    st.subheader("Mario Pimienta Calderón - Máster en Bioinformática - Universidad de Murcia")

    # Descripción general
    st.write("Esta aplicación te permite analizar y comparar proteínas de diversas maneras. A continuación, se explican cada una de las funcionalidades disponibles:")

    # Funcionalidad 1: Plegar proteínas a partir de un fichero FASTA o MultiFASTA
    st.header("📄 Plegar proteínas desde FASTA o multiFASTA")
    st.write(
        """
        Con esta funcionalidad, puedes cargar un fichero FASTA o MultiFASTA y obtener la estructura plegada de una o varias proteínas.
        - **Paso 1:** Sube tu fichero FASTA o MultiFASTA.
        - **Paso 2:** Procesa el archivo para obtener las estructuras plegadas de las proteínas.
        """
    )

    # Funcionalidad 2: Comparar una proteína local con una de UniProt
    st.header("🔍 Comparar proteína local con UniProt")
    st.write(
        """
        Aquí puedes comparar una proteína de un fichero FASTA local con una proteína de UniProt utilizando su accesion ID.
        - **Paso 1:** Sube tu fichero FASTA local.
        - **Paso 2:** Introduce el accesion ID de la proteína de UniProt.
        - **Paso 3:** Se expandirá una tabla con 30 variantes patogénicas relacionadas (puede demorarse varios minutos) y podrás seleccionar una para comparar.
        """
    )

    # Funcionalidad 3: Comparar dos proteínas a partir de UniProt ID y una mutación
    st.header("⚖️ Comparar proteínas a partir de UniProt ID y mutación")
    st.write(
        """
        Esta funcionalidad te permite comparar dos proteínas utilizando el ID de UniProt y una mutación introducida.
        - **Paso 1:** Introduce el UniProt ID de la proteína. Se abrirá un pop-up con su información de Uniprot.
        - **Paso 2:** Introduce la mutación con el formato A123B.
        - **Paso 3:** Compara las dos proteínas para observar las diferencias causadas por la mutación.
        """
    )

    st.header("Bibliografía")
    st.write("""
    EJEMPLOS:
    
    1. Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H., ... & Bourne, P. E. (2000). The Protein Data Bank. Nucleic acids research, 28(1), 235-242.
    2. Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., ... & Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596(7873), 583-589.
    """)

if __name__ == "__main__":
    main()
