from Bio import SeqIO
from io import StringIO

def parse_fasta_file(fasta_content) -> dict:
    protein_dict = {}
    for record in SeqIO.parse(fasta_content, "fasta"):
        protein_name = record.id
        protein_sequence = str(record.seq)
        protein_dict[protein_name] = protein_sequence
    return protein_dict
