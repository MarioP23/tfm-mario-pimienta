
# Protein Structure Prediction and Visualization

Este proyecto permite la comparación y visualización de estructuras proteicas utilizando datos PDB. Utiliza la biblioteca Py3Dmol para la visualización en 3D y Streamlit para la creación de la interfaz de usuario.

## Tabla de Contenidos
1. [Instalación](#instalación)
2. [Uso](#uso)
3. [Estructura del Proyecto](#estructura-del-proyecto)
4. [Funcionalidades](#funcionalidades)
5. [Ejemplo de Uso](#ejemplo-de-uso)
6. [Contribuciones](#contribuciones)
7. [Licencia](#licencia)

## Instalación

### Prerrequisitos
- Python 3.7 o superior
- Pip (administrador de paquetes de Python)

### Instalación de Dependencias

1. Clona este repositorio:
    \`\`\`bash
    git clone https://github.com/tu-usuario/protein-structure-visualization.git
    cd protein-structure-visualization
    \`\`\`

2. Crea un entorno virtual (opcional pero recomendado):
    \`\`\`bash
    python -m venv env
    source env/bin/activate  # En Windows usa `env\Scriptsctivate`
    \`\`\`

3. Instala las dependencias necesarias:
    \`\`\`bash
    pip install -r requirements.txt
    \`\`\`

### Instalación de Py3Dmol y Streamlit

Si no están incluidos en \`requirements.txt\`, instala Py3Dmol y Streamlit:
\`\`\`bash
pip install py3Dmol streamlit
\`\`\`

## Uso

Para iniciar la aplicación de Streamlit, navega al directorio del proyecto y ejecuta:
\`\`\`bash
streamlit run streamlit_app/main.py
\`\`\`

## Estructura del Proyecto

\`\`\`
protein-structure-visualization/
│
├── data/
│   ├── download_pdb_file_rscb.py      # Función para descargar archivos PDB desde el RSCB
│
├── parsers/
│   ├── fasta_parser.py                # Funciones para parsear archivos FASTA
│
├── show_pdb/
│   ├── display_protein.py             # Funciones para visualizar proteínas con Py3Dmol
│
├── esmfold_api_tool/
│   ├── esmfold.py                     # Funciones para interactuar con la API ESMFold
│
├── streamlit_app/
│   ├── main.py                        # Archivo principal de Streamlit
│   ├── pages/
│       ├── 2_RSCB_PDB_viewer.py       # Página de comparación de estructuras PDB
│       ├── 1_Protein_Structure_Prediction.py # Página de predicción y visualización de estructuras proteicas
│
├── requirements.txt                   # Dependencias del proyecto
├── README.md                          # Documentación del proyecto
\`\`\`

## Funcionalidades

### Comparación de Estructuras PDB
- **Entrada por código RSCB PDB**: Permite al usuario introducir un código PDB para descargar y visualizar una estructura.
- **Carga de archivos PDB**: Permite al usuario cargar su propio archivo PDB para la visualización.
- **Personalización de la Visualización**: Ofrece opciones para personalizar la visualización, incluyendo estilo, opacidad de la superficie, color de la proteína y color de fondo.

### Predicción y Visualización de Estructuras Proteicas
- **Carga de archivos FASTA**: Permite al usuario cargar un archivo FASTA para la predicción de la estructura proteica.
- **Integración con ESMFold**: Utiliza la API ESMFold para predecir la estructura proteica a partir de la secuencia FASTA.
- **Visualización 3D**: Muestra la estructura 3D predicha de la proteína.
- **Descarga de archivos PDB**: Permite descargar la estructura PDB predicha.

## Ejemplo de Uso

### Comparación de Estructuras PDB

1. Inicia la aplicación Streamlit:
    \`\`\`bash
    streamlit run streamlit_app/pages/2_RSCB_PDB_viewer.py
    \`\`\`

2. Introduce un código RSCB PDB y presiona "Buscar por código RSCB PDB".
3. Sube un archivo PDB desde tu máquina.
4. Presiona "Personalizar display" y ajusta las opciones de visualización según tus preferencias.

### Predicción y Visualización de Estructuras Proteicas

1. Inicia la aplicación Streamlit:
    \`\`\`bash
    streamlit run streamlit_app/pages/1_Protein_Structure_Prediction.py
    \`\`\`

2. Sube un archivo FASTA.
3. Presiona "Ejecutar" para predecir la estructura proteica utilizando ESMFold.
4. Selecciona la proteína que deseas visualizar y ajusta las opciones de visualización según tus preferencias.
5. Descarga el archivo PDB si lo deseas.
