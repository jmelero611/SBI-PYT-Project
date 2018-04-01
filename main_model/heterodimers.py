##PYTHON SBI PROJECT

########################################################
######################HETERODIMERS #####################
########################################################

#Simple input file with iterations between heterodimers
from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import string
from common_functions import *

###################### FUNCTIONS ########################

##Function to align sequences from chains
def align_sequences_heterodimers(structure):
	"""
	Align all sequences generated from pdb files and get the best alignments.
	 input: structure from which we want to perform the alignments
	 output: the ids of the best alignments
	"""

	chains = list(structure.get_chains())
	best_aln = []
	n = len(chains)

	#Select sequences to align
	for i in range(n-1):
		id1 = chains[i].get_full_id()[1:]
		seq1 = get_sequence(chains[i])
		
		for j in range(i+1,n):
			id2 = chains[j].get_full_id()[1:]
			seq2 = get_sequence(chains[j])	
			#print("align %s with %s" %(id1, id2))
			#Align sequence
			alns = pairwise2.align.globalxx(seq1, seq2)
			score = alns[0][2]

			#Calculate sequence identity
			ident_perc = score / max(len(seq1), len(seq2))
			#print(ident_perc)
			if ident_perc > 0.99:
				print("Sequence of %s has aligned with %s with a identity of %d\n" %(id1, id2, ident_perc))
				best_aln.append((id1,id2))
	
	return best_aln

## Function to superimpose sequences aligned
def superimpose_structures_heterodimers(structure, best_aln):
	"""
	Superimpose structures from the best alginments.
	 input: structure with all the chains and a list with the chains we want to superimpose
	 output: structure with the correct coordinates of the macrocomplex
	"""
	sup = Superimposer()
	new_structure = structure.copy()

	#get fixed and moving atoms from best alginments
	for aln1, aln2 in best_aln:
		ref_atoms = list(new_structure[aln1[0]][aln1[1]].get_atoms())
		alt_atoms = list(new_structure[aln2[0]][aln2[1]].get_atoms())

		#Align these paired atoms
		sup.set_atoms(ref_atoms, alt_atoms) 
		sup.apply(new_structure[aln2[0]].get_atoms())
		new_structure[aln1[0]].detach_child(aln1[1])
		print("RMS(model %s, model %s) = %0.2f" % (aln1, aln2, sup.rms))
		
	return new_structure
