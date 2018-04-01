#########################################################
################### COMMON FUNCTIONS ####################
#########################################################

from Bio.PDB import *
from Bio import pairwise2

#Get sequence from chain taking into account alpha carbons
def get_sequence(chain):
	"""
	Get sequence from chain with its alpha carbons
	""" 
	ppb = CaPPBuilder()
	residues_to_remove = []
	for residue in chain:
		if residue.id[0] != ' ':
			residues_to_remove.append(residue.id) #select heteroatoms
	
	for residue in residues_to_remove:
		chain.detach_child(residue) #remove heteroatoms

	if len(chain) <= 30:
		return None	#select chains are suitable to be protein chains
	else:
		pp = ppb.build_peptides(chain)
		seq = pp[0].get_sequence()
		return seq


#Compare chain sequences
def seq_comparison(seq1, seq2):
	"""
	Use Pairwise alignment to align to sequences.
	Return True if similar and False if not.
	"""
	alignment = pairwise2.align.globalxx(seq1, seq2)
	score = alignment[0][2]
	ident_perc = score / max(len(seq1), len(seq2))
	#print("%s and %s has a identity of %f" %(seq1, seq2, ident_perc))

	if ident_perc > 0.99:
		return True
	else:
		return False

#Clash identifier
def clash_identifier(structure):
	"""
	Check if there is a clash within a structure
	It returns 1 if there is no clash and 0 if there is a clash.
	 inputs: structure we want to evaluate
	 output: number of clashes
	"""

	clash = 0
	ns = NeighborSearch(list(structure.get_atoms()))
	for chain in structure.get_chains():
		for atom in chain.get_atoms():
			for chain2 in ns.search(atom.get_coord(), 2, level='C'):
				if chain2 == chain:
					clash += 1
	if clash == 0:
		print("THERE is no clash")
	return clash