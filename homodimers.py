# Copyright 2007 Peter Cock, all rights reserved.
# Licenced under the GPL v2 (or later at your choice)
#
# Please see this website for details:
# http://www.warwick.ac.uk/go/peter_cock/python/protein_superposition/
from Bio.PDB import *
import argparse
import numpy
import os
import sys
import re

parser = argparse.ArgumentParser(description="""This program reconstructs the macrocomplex from protein interactions pdb files""")

parser.add_argument('-i', '--input',
					dest = "infile",
					action = "store",
					default = os.getcwd(),
					help = "Enter a list of interaction files, a directory or the current directory will be selected")

options = parser.parse_args()

#Extract input files

motif = re.compile('.pdb$')
path = options.infile

	#check if input is a directory
if os.path.isdir(path):
	pdb_files = [f for f in os.listdir(path) if motif.search(f) is not None]
	os.chdir(path)
else:
	pdb_files = [options.infile]


n = 0
structures = Structure.Structure("homodimer")
p = PDBParser()
pdb_out_filename = "test.pdb"


for f in pdb_files:	
	str = p.get_structure(n, f)

	for model in str:
		if n != 0:
			model.id = n
		for chain in model:
			for residue in chain:
				if residue.id[0] != ' ':
					chain.detach_child(residue.id) #delete heteroatoms
			if len(chain) <= 10: #select chains which forms sequences
				str.detach_child(chain)
		structures.add(model)
	n += 1

#print(list(structures.get_models()))
print("Everything aligned to first model...")
ref_model = structures[0]
#print(ref_model)
	
for alt_model in structures:
	#print(alt_model)
	#Build paired lists of c-alpha atoms:
	ref_atoms = []
	alt_atoms = []
	for (ref_chain, alt_chain) in zip(ref_model, alt_model):
		#print(ref_chain)
		for ref_res, alt_res in zip(ref_chain, alt_chain):
			if ref_res.resname == alt_res.resname and ref_res.id == alt_res.id:
				ref_atoms.append(alt_res['CA'])
				alt_atoms.append(alt_res['CA'])

		#Align these paired atom lists:
	super_imposer = Superimposer()
	super_imposer.set_atoms(ref_atoms, alt_atoms)

	if ref_model.id == alt_model.id :
        #Check for self/self get zero RMS, zero translation
        #and identity matrix for the rotation.
		assert numpy.abs(super_imposer.rms) < 0.0000001
		assert numpy.max(numpy.abs(super_imposer.rotran[1])) < 0.000001
		assert numpy.max(numpy.abs(super_imposer.rotran[0]) - numpy.identity(3)) < 0.000001
	else :
        #Update the structure by moving all the atoms in
        #this model (not just the ones used for the alignment)
		super_imposer.apply(alt_model.get_atoms())

	print("RMS(first model, model %i) = %0.2f" % (alt_model.id, super_imposer.rms))


print("Saving aligned structure as PDB file %s" % pdb_out_filename)

io=PDBIO()
io.set_structure(structures)
io.save(pdb_out_filename)

#print("Done")

