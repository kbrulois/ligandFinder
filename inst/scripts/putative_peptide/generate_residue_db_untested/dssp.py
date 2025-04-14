from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os
import datatable as dt
import numpy as np
import re



input = '/Users/kbrulois/peptide_alg/UP000005640_9606_HUMAN_v4'
output = /Users/kbrulois/peptide_alg/alphafold_dssp'

pdb_files = [f for f in os.listdir(input) if f.endswith('.pdb')]

for f in pdb_files[0:10]:
	file_name = "".join([input, '/', f])
	out_file = "".join([output, '/', re.sub(r".pdb$", ".csv", f)])
	p = PDBParser()
	structure = p.get_structure("1MOT", file_name)
	model = structure[0]
	dssp = DSSP(model, file_name)
	all_keys = list(dssp.keys())
	dat = [dssp[x] for x in all_keys]
	dat = np.vstack(dat)
	dat = dt.Frame(dat)
	dat.to_csv(out_file)