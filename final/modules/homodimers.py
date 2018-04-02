#####################################################
#################### HOMODIMERS #####################
#####################################################

from Bio.PDB import *
from Bio import pairwise2
from common_functions import *
import numpy
import string
import copy

#Function to build the complex for homodimers from all files
def get_structure_homodimer(structure):
	"""
	Builds the structure when the input are homodimers.
	Returns the complete predicted structure.
	"""
	sup = Superimposer()
	new_structure = structure.copy()
	ref_model = new_structure[0]

	#Iterate over all interactions
	for alt_model in new_structure:
		#Build paired lists of c-alpha atoms:
		ref_atoms = []
		alt_atoms = []
		for (ref_chain, alt_chain) in zip(ref_model, alt_model):
			for ref_res, alt_res in zip(ref_chain, alt_chain):
				if ref_res.resname == alt_res.resname and ref_res.id == alt_res.id:
					ref_atoms.append(alt_res['CA'])
					alt_atoms.append(alt_res['CA'])

		#Align these paired atom lists:
		sup.set_atoms(ref_atoms, alt_atoms)

		if ref_model.id == alt_model.id :
	        #Check for self/self get zero RMS, zero translation
	        #and identity matrix for the rotation.
			assert numpy.abs(sup.rms) < 0.0000001
			assert numpy.max(numpy.abs(sup.rotran[1])) < 0.000001
			assert numpy.max(numpy.abs(sup.rotran[0]) - numpy.identity(3)) < 0.000001

		else :
	        #Update the structure by moving all the atoms in
	        #this model (not just the ones used for the alignment)
			sup.apply(alt_model.get_atoms())
			alt_model.detach_child(list(alt_model.get_chains())[0].id)

		print("RMS(first model, model %i) = %0.2f" % (alt_model.id, sup.rms))

	return new_structure


#Function to build the complex for homodimers from n interactions
def get_structure_mix(all_interactions, homodimers, heterodimers):
	"""
	Predict the macrocomplex fromed by heterodimers and homodimers.
	Superposes the homodimer interactions to get the right coordinates of the chains.
	Adds the heterodimers interactions within the complex
	Builds and returns the new structure.
	 input: two dictionaries with the structures from both homodimers and heterodimers interactins respectively.
	 ouput: new structure formed by those interactions.
	"""

	super_imposer = Superimposer()
	sup_chains = {}

	#Get homodimer chains
	n = 0
	for str_id, str in homodimers.items():
		chains = all_interactions[str_id]
		ref_model = chains[0]

		sup_chains[n] = [ref_model[0],ref_model[1]]

		for alt_model in chains[1:]:
			ref_atoms = []
			alt_atoms = []
			print("Comparing %s with %s" %(ref_model, alt_model))
			for (ref_chain, alt_chain) in zip(ref_model, alt_model):
				for ref_res, alt_res in zip(ref_chain, alt_chain):
					if ref_res.resname == alt_res.resname and ref_res.id == alt_res.id:
						ref_atoms.append(alt_res['CA'])
						alt_atoms.append(alt_res['CA'])

			#Align these paired atom lists:
			super_imposer.set_atoms(ref_atoms, alt_atoms)

		    #Update the structure by moving all the atoms in
		    #this model (not just the ones used for the alignment)
			super_imposer.apply(list(alt_model[0].get_atoms()) + list(alt_model[1].get_atoms()))
			sup_chains[n].append(ref_model[1])

			print("RMS(ref %s, model %s) = %0.2f" % (ref_model, alt_model, super_imposer.rms))
		n += 1
	
	#Add heterodimer interactions
	for str_id, str in heterodimers.items():
		n += 1
		chains = all_interactions[str_id]
		alt_model = chains[0]
		seq1 = get_sequence(alt_model[0])
		seq2 = get_sequence(alt_model[1])
		for n_inter in sup_chains:
			ref_chain = sup_chains[n_inter][0]
			ref_seq = get_sequence(ref_chain)
			if (seq_comparison(ref_seq, seq1)):
				ref_atoms = list(alt_model[0].get_atoms())
				alt_atoms = list(ref_chain.get_atoms())
				super_imposer.set_atoms(ref_atoms, alt_atoms)
				super_imposer.apply(list(alt_model[0].get_atoms()) + list(alt_model[1].get_atoms()))
				sup_chains[n] = [alt_model[1]]
				print("RMS(ref %s, model %s) = %0.2f" % (chain1, chain, super_imposer.rms))
				break

			elif seq_comparison(ref_seq, seq2):
				ref_atoms = list(alt_model[1].get_atoms())
				alt_atoms = list(ref_chain.get_atoms())
				super_imposer.set_atoms(ref_atoms, alt_atoms)
				super_imposer.apply(list(alt_model[0].get_atoms()) + list(alt_model[1].get_atoms()))
				sup_chains[n] = [alt_model[0]] 
				print("RMS(ref %s, model %s) = %0.2f" % (alt_model, ref_chain, super_imposer.rms))
				break
			else:
				continue


	#Create the new structure
	new_structure = Structure.Structure("macrocomplex")
	new_structure.add(Model.Model(0))

	n = 0
	for chains in sup_chains.values():
		for chain in chains:
			try:
				new_structure[0].add(chain)
			except:
				n += 1
				new_structure.add(Model.Model(n))
				new_structure[n].add(chain)

	return new_structure



