def parse_plddt_from_pdb(pdb_data):
    """Parses pLDDT scores from PDB data and returns a dictionary of residue index to pLDDT score."""
    plddt_scores = {}
    for line in pdb_data.splitlines():
        if line.startswith("ATOM"):
            resi = int(line[22:26].strip())
            plddt = float(line[60:66].strip())
            plddt_scores[resi] = plddt
    return plddt_scores