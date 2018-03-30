#####################################################
#################### HOMODIMERS #####################
#####################################################

from Bio.PDB import *
import numpy
import string
import copy

#Function to build the complex for homodimers from all files
def get_structure_homodimer_from_files(structures):
	"""
	Predict the structure of a homodimer of the same chains and heterodimer of two chains
	"""
	print("Aligned all to the first structure")
	ref_model = structures[0].copy()

	#Iterate over all interactions
	for alt_model in structures:

		#Build paired lists of c-alpha atoms:
		ref_atoms = []
		alt_atoms = []
		for (ref_chain, alt_chain) in zip(ref_model, alt_model):

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

	return structures





#Function to build the complex for homodimers from n interactions
def get_structure_homodimer_from_number_interaction(main_inter, n):
	"""
	Predict the structure of a homodimer of the same chains and heterodimer of two chains
	"""

	super_imposer = Superimposer()
	inter = main_inter[0].copy()

	new_structure = Structure.Structure("New")
	new_structure.add(Model.Model('ref'))
	chain1 = list(inter.get_chains())[0]
	new_structure['ref'].add(chain1.copy())

	for i in range(n):
		new_structure.add(Model.Model(i))
		alt_chain = list(new_structure.get_chains())[-1]
		print(alt_chain.get_full_id())
		#print("Superposing %s with %s" %(ref_model.get_full_id(), alt_model.get_full_id()))
		super_imposer.set_atoms(list(alt_chain.get_atoms()), list(list(inter.get_chains())[0].get_atoms()))
		#Align these paired atom lists:
		print(super_imposer.rotran)
		super_imposer.apply(list(inter.get_atoms()))
		chain2 = list(inter.get_chains())[1]
		new_structure[i].add(chain2.copy())

	return new_structure

