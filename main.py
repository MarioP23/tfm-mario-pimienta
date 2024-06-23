from parsers import fasta_parser
from esmfold_api_tool import esmfold
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

fasta_filepath = fasta_file = "test/data/HTRA2_WT.fasta"

logging.info(msg=f'Parsing FASTA file: {fasta_filepath}')

fasta_protein_dict = fasta_parser.parse_fasta_file(fasta_file_path=fasta_filepath)

logging.info(msg=f'Parsed file: {fasta_filepath}')
logging.info(msg=f'Sequences included: {list(fasta_protein_dict.keys())}')

for protein_name in list(fasta_protein_dict.keys()):
  logging.info(msg=f'Folding protein sequence {protein_name} with ESMFold2 via API')
  esmfold.fold_sequence(fasta_seq=fasta_protein_dict[protein_name], 
                        output_pdb_filepath=f'data/output/{protein_name}.pdb')

logging.info('All sequences folded succesfully')

pathogenic_variants_htra = ['rs1407675367', 'rs767006508', 'rs1675486182'] #Se han seleccionado los 3 SNPs considerados patogenicos en clinvar. 



