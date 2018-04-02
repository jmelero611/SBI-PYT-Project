#########################################################
################### COMMON FUNCTIONS ####################
#########################################################

from Bio.PDB import *
from Bio import pairwise2

#Get sequence from chain taking into account alpha carbons
def get_sequence(chain):
	"""
	Deletes heteroatoms and returns the sequence of proteins with more than 30 aminoacids from their alpha carbons. 
	It uses CaPPBuilder package from Bio.PDB.
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
	Uses pariwise alignment to align a pair of sequences to find if two sequences are similar or not.
	Returns True if percentage of identity is greater than 95%, and False if it is not.
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
	Returns the number of clashes and if there is no clash, prints it.
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
		print("There is no clash")
	return clash
