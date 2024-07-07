import biotite.sequence as seq
import biotite.sequence.align as align
import numpy as np
import matplotlib.pyplot as plt
import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.sequence.graphics as graphics
from io import BytesIO

def align_sequences(query_seq_str, hit_seq_str):
    query_seq_str = query_seq_str.split('\n', 1)[-1].replace(' ', '').replace('\n', '')
    hit_seq_str = hit_seq_str.split('\n', 1)[-1].replace(' ', '').replace('\n', '')

    query_seq = seq.ProteinSequence(query_seq_str)
    hit_seq = seq.ProteinSequence(hit_seq_str)
    matrix = align.SubstitutionMatrix.std_protein_matrix()
    GAP_PENALTY = (-12, -1)
    alignment = align.align_optimal(
        query_seq, hit_seq, matrix,
        local=True, gap_penalty=GAP_PENALTY, max_number=1
    )[0]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    graphics.plot_alignment_similarity_based(
        ax, alignment, matrix=matrix, labels=["Query", "Hit"],
        show_numbers=True, show_line_position=True
    )
    fig.tight_layout()
    img_buf = BytesIO()
    plt.savefig(img_buf, format='png')
    plt.close(fig)
    img_buf.seek(0)
    
    return alignment.score, img_buf

def ramachandran_plot(pdb_filepath):
    try:
        array = strucio.load_structure(pdb_filepath)
      
        if array.coord.shape[0] < 3:
            raise ValueError("La estructura tiene menos de tres átomos, no se pueden calcular todas las medidas geométricas.")
      
        phi, psi, omega = struc.dihedral_backbone(array)
      
        plt.figure()
        plt.plot(phi * 360/(2*np.pi), psi * 360/(2*np.pi), marker="o", linestyle="None")
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        plt.xlabel("$\phi$")
        plt.ylabel("$\psi$")
      
        img_buf = BytesIO()
        plt.savefig(img_buf, format='png')
        plt.close()
        img_buf.seek(0)
      
        return img_buf
    except Exception as e:
        raise ValueError(f"Error en ramachandran_plot: {e}")