##PYTHON SBI PROJECT

########################################################
######################HETERODIMERS #####################
########################################################

#Simple input file with iterations between heterodimers
from Bio.PDB import *
from Bio import pairwise2
import sys
from common_functions import *

###################### FUNCTIONS ########################

##Function to align sequences from chains
def align_sequences_heterodimers(structure, option):
	"""
	Align all sequences generated from PDB files and returns the best alignments.
	 input: structure from which we want to perform the alignments
	 output: the ids of the best alignments
	"""

	chains = list(structure.get_chains())
	best_aln = []
	n = len(chains)

	if option:
		sys.stderr.write("Align chain sequences: \n")
	#Select sequences to align
	for i in range(n-1):
		id1 = chains[i].get_full_id()[1:]
		seq1 = get_sequence(chains[i])
		
		for j in range(i+1,n):
			id2 = chains[j].get_full_id()[1:]
			seq2 = get_sequence(chains[j])	

			#Align sequence
			alns = pairwise2.align.globalxx(seq1, seq2)
			score = alns[0][2]

			#Calculate sequence identity
			ident_perc = score / max(len(seq1), len(seq2))

			if ident_perc > 0.99:
				best_aln.append((id1,id2))
				if option:
					sys.stderr.write("Sequence of %s has aligned with %s with a identity of %d\n" %(id1, id2, ident_perc))
	
	return best_aln

## Function to superimpose sequences aligned
def superimpose_structures_heterodimers(structure, best_aln, option):
	"""
	Superimposes the best alignments from the input structures, gets the chains to create the new structure and returns the new structure superimposed.
	 input: structure with all the chains and a list with the chains we want to superimpose
	 output: structure with the correct coordinates of the macrocomplex
	"""
	sup = Superimposer()
	structure = structure.copy()
	new_structure = Structure.Structure("macrocomplex")
	new_structure.add(Model.Model(0))
	n = 0

	if option:
		sys.stderr.write("Performing the superimposition with heterodimer structure with best alignments: \n")

	#get fixed and moving atoms from best alginments
	for aln1, aln2 in best_aln:
		ref_atoms = []
		alt_atoms = []
		for ref_res, alt_res in zip(structure[aln1[0]][aln1[1]], structure[aln2[0]][aln2[1]]):
			if ref_res.resname == alt_res.resname and ref_res.id == alt_res.id:
				ref_atoms.append(ref_res['CA'])
				alt_atoms.append(alt_res['CA'])

		#Align these paired atoms
		sup.set_atoms(ref_atoms, alt_atoms) 
		sup.apply(structure[aln2[0]].get_atoms())

		#Add chain if there is no clashes
		add_chain = structure[aln2[0]][aln2[1]]
		
		if bool(list(new_structure.get_chains())):
				if not clash_identifier(new_structure, add_chain):
					if option:
						sys.stderr.write("\tRMSD of the model %s with respect model %s is %0.2f" % (aln1, aln2, sup.rms))
					try:
						new_structure[n].add(add_chain)
					except:
						n += 1
						new_structure[n].add(add_chain)
		else:
			new_structure[n].add(add_chain)
			if option:
				sys.stderr.write("RMSD of the model %s with respect model %s is %0.2f" % (aln1, aln2, sup.rms))

	return new_structure