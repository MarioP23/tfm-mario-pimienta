import py3Dmol
from stmol import showmol

def show_protein_grid(pdb_data_1: str, pdb_data_2):
    """Muestra dos proteínas en una cuadrícula usando py3Dmol.
    """
    pdbview = py3Dmol.view(width=1250, height=600, viewergrid=(1, 2), linked=True)
    
    # Configurar la primera proteína
    pdbview.addModel(pdb_data_1, 'pdb', viewer=(0, 0))
    pdbview.setStyle({'model': -1, 'cartoon': {'color': 'red'}}, viewer=(0, 0))
    pdbview.zoomTo(viewer=(0, 0))
    pdbview.addSurface(py3Dmol.VDW, {'opacity': 0.25}, viewer=(0, 0))

    # Configurar la segunda proteína
    pdbview.addModel(pdb_data_2, 'pdb', viewer=(0, 1))
    pdbview.setStyle({'model': -1, 'cartoon': {'color': 'cyan'}}, viewer=(0, 1))
    pdbview.zoomTo(viewer=(0, 1))
    pdbview.addSurface(py3Dmol.VDW, {'opacity': 0.25}, viewer=(0, 1))
    
    pdbview.setBackgroundColor("white")
    pdbview.spin(False)
    
    showmol(pdbview, height=600, width=1250)


def render_protein_3d(pdb_data):
    """Renders a 3D visualization of a protein given its PDB data."""
    pdbview = py3Dmol.view(width=1250, height=600)
    pdbview.addModel(pdb_data, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(1, 800)
    pdbview.spin(False)
    pdbview.addSurface(py3Dmol.VDW, {'opacity': 0.5})
    showmol(pdbview, height=950, width=1300)
