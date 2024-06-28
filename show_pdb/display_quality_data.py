import numpy as np
from Bio.PDB import PDBParser, Superimposer
from io import StringIO
from prody import *

def calculate_rmsd_from_strings(pdb_str_1, pdb_str_2):
    parser = PDBParser(QUIET=True)

    # Leer las estructuras desde los strings
    structure1 = parser.get_structure('structure1', StringIO(pdb_str_1))
    structure2 = parser.get_structure('structure2', StringIO(pdb_str_2))

    # Seleccionar el primer modelo de cada estructura
    model1 = structure1[0]
    model2 = structure2[0]

    # Extraer los átomos CA (Carbono Alpha) de ambas estructuras
    atoms1 = [atom for atom in model1.get_atoms() if atom.get_id() == 'CA']
    atoms2 = [atom for atom in model2.get_atoms() if atom.get_id() == 'CA']

    if len(atoms1) != len(atoms2):
        raise ValueError("Las estructuras no tienen el mismo número de átomos CA.")

    # Crear un Superimposer para alinear las estructuras
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)

    # Calcular el RMSD
    rmsd = super_imposer.rms

    return rmsd

def align_and_compare(pdb_str_1, pdb_str_2):
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure('structure1', StringIO(pdb_str_1))
    structure2 = parser.get_structure('structure2', StringIO(pdb_str_2))

    atoms1 = [atom for atom in structure1.get_atoms() if atom.get_id() == 'CA']
    atoms2 = [atom for atom in structure2.get_atoms() if atom.get_id() == 'CA']

    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(structure2.get_atoms())

    return super_imposer.rms, prody.calcTMscore(atoms1, atoms2)

def get_prody_matches(pdb1, pdb2):
    pdb1_io = StringIO(pdb1)
    pdb2_io = StringIO(pdb2)
    atom_group_1 = parsePDBStream(pdb1_io)
    atom_group_2 = parsePDBStream(pdb2_io)
    return matchChains(atom_group_1, atom_group_2)


