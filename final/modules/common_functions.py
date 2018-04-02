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
	 input: chain from which we want to get the sequence
	 output: sequence object of the chain
	""" 
	ppb = CaPPBuilder()
	residues_to_remove = []
	for residue in chain:
		if residue.id[0] != ' ':
			residues_to_remove.append(residue.id) #select heteroatoms
	
	for residue in residues_to_remove:
		chain.detach_child(residue) #remove heteroatoms

	if len(chain) <= 30:
		return None	#select chains that are suitable to be protein chains
	else:
		pp = ppb.build_peptides(chain)
		seq = pp[0].get_sequence()
		return seq


#Compare chain sequences
def seq_comparison(seq1, seq2):
	"""
	Uses pariwise alignment to align a pair of sequences to find if two sequences are similar or not.
	Returns True if percentage of identity is greater than 95%, and False if it is not.
	 input: both sequences we want to align
	 output True or False dependend of the percentage od identity
	"""
	alignment = pairwise2.align.globalxx(seq1, seq2)
	score = alignment[0][2]
	ident_perc = score / max(len(seq1), len(seq2))

	if ident_perc > 0.99:
		return True
	else:
		return False


#Compare chains from two structures
def chains_comparison(str1, str2):
	"""
	Compares both chains from the input structures to know if they are similar or not
	 input: both structures from which we want to perform the chains_comparison
	 output: True if both structures have the same chains and False if they are differents
	"""
	rs = 0
	rs_tot = [[0,0],[0,0]]
	for chain1 in str1.get_chains():
		seq1 = get_sequence(chain1)
		for chain2 in str2.get_chains():
			seq2 = get_sequence(chain2)

			rs2 = 0
			if seq_comparison(seq1, seq2) is True:
				rs_tot[rs][rs2] = 1
			else:
				rs_tot[rs][rs2] = 0
		rs += 1

	if 1 in rs_tot[0] and 1 in rs_tot[1]:
		return True
	else:
		return False

#Create temporal structure
def temp_structure(interaction_dic, name):
	"""
	Builds a temporary structure and fill it with the interactions.
	 input: dictionary with all pairwise interaction chains and the name of the structure
	 output: new structure
	"""
	
	new_structure = Structure.Structure(name)

	i = 0
	for interaction in interaction_dic.values():
		for inter in interaction:
			new_structure.add(Model.Model(i))
			new_structure[i].add(inter[0])
			new_structure[i].add(inter[1])
			i += 1

	return new_structure

#Clash identifier
def clash_identifier(structure, add_chain):
	"""
	Check if there is a clash within a structure
	Returns the number of clashes and if there is no clash, prints it.
	 inputs: structure we want to evaluate, chain we want to add
	 output: number of clashes
	"""

	ns = NeighborSearch(list(structure.get_atoms()))

	for atom in add_chain.get_atoms():
		for chain in ns.search(atom.get_coord(), 1, level='C'):
			return True

	return False

#Save structure into a new pdb file
def save_complex(structure, outfile):
	"""Saves the structure into a PDB file. Uses Bio.PBIO package.
	 input: structure from which we want to create a pdb file and the name of the file.
	 ouput: new pdb file.
	"""
	io = PDBIO()
	io.set_structure(structure)
	io.save(outfile)