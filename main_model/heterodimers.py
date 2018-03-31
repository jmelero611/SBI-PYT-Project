##PYTHON SBI PROJECT

########################################################
######################HETERODIMERS #####################
########################################################

#Simple input file with iterations between heterodimers
from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import string

###################### FUNCTIONS ########################
#Get sequence (CA) from chain
def get_sequence(chain):
	ppb = CaPPBuilder()

	for residue in chain:
		if residue.id[0] != ' ':
			chain.detach_child(residue.id) #delete heteroatoms
		if len(chain) <= 40:
			return None	#select chains which forms sequences from more than 
		else:
			pp = ppb.build_peptides(chain)
			seq = pp[0].get_sequence() #get the
			return seq

##Function to align sequences from chains
def align_sequences_heterodimers(structure):
	"""Align all sequences generated from pdb files and get the best alignments """

	chains = list(structure.get_chains())
	best_aln = []
	n = len(chains)

	#Select sequences to align
	for i in range(n-1):
		id1 = chains[i].get_full_id()[1:]
		print(id1)
		seq1 = get_sequence(chains[i])
		
		for j in range(i+1,n):
			id2 = chains[j].get_full_id()[1:]
			seq2 = get_sequence(chains[j])	

			#Align sequence
			alns = pairwise2.align.globalxx(seq1, seq2)
			score = alns[0][2]

			#Calculate sequence identity
			ident_perc = score / max(len(seq1), len(seq2))
			#print(ident_perc)
			if ident_perc > 0.95:
				print("Sequence of %s has aligned with %s with a identity of %d\n" %(id1, id2, ident_perc))
				best_aln.append((id1,id2))

	return best_aln

## Function to superimpose sequences aligned
def superimpose_structures_heterodimers(structures, best_aln):
	"""
	Superimpose structures.
	Get the chains to create the new structure.
	Save the structure in a new pdb file.
	"""
	print("Perform the superimposition")
	sup = Superimposer() #initialize superimpose
	for aln1, aln2 in best_aln:
		print(aln1)
		#print("first aln %s, correspond to chain %s" %(aln1, structures[0][aln1]))
		#print("second aln %s, correspond to chain %s" %(aln2, structures[0][aln2]))
		#get fixed and moving atoms for superimposition
		fixed = list(structures[aln1[0]][aln1[1]].get_atoms())
		moving = list(structures[aln2[0]][aln2[1]].get_atoms())

		#apply superimposition
		sup.set_atoms(fixed, moving) # set coords
		sup.apply(structures[aln2[0]][aln2[1]]) # apply to structure-chain objecy
		
		print("RMS(model %s, model %s) = %0.2f" % (aln1, aln2, sup.rms))
		#Delete the fixed chain from structure
		structures[aln1[0]].detach_child(aln1[1])

	return structures
