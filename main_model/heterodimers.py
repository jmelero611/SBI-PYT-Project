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
def align_sequences(structures):
	"""Align all sequences generated from pdb files and get the best alignments """

	##Function to calculate the percentaje of identity inside the alignment function
	def calculate_identity(seqA, seqB):
		"""
		Returns the precentage of identical characters between two sequences.
		Assumes the sequences are aligned.
		"""
		sa, sb, sl = seqA, seqB, len(seqA)
		matches = [sa[i] == sb[i] for i in range(sl)]
		seq_id = (100 * sum(matches)) / sl

		gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
		gap_id = (100 * sum(matches)) / gapless_sl
        
		return (seq_id, gap_id)

    #Return to main function
	print("Align sequences")
	
	matrix = matlist.blosum62
	gap_open = -10.0
	gap_extend = -0.5
	best_aln = []
	all_aln = []

	chains = list(structures.get_chains())

	n = len(chains)

	#Select sequences to align
	for i in range(n-1):
		id1 = chains[i].get_id()
		seq1 = get_sequence(chains[i])
		
		for j in range(i+1,n):
			id2 = chains[j].get_id()
			seq2 = get_sequence(chains[j])	

			#Align sequence
			alns = pairwise2.align.globalds(seq1, seq2, 
											matrix, gap_open, gap_extend,
											penalize_end_gaps = (False, False))
			top_aln = alns[0]
			al1, al2, score, begin, end = top_aln

			#Calculate sequence identity
			seq_id, g_seq_id = calculate_identity(al1, al2)

			all_aln.append((id1, id2, seq_id))

			if seq_id > 99:
				#print("Sequence of %s has aligned with %s with a identity of %d\n" %(id1, id2, seq_id))
				best_aln.append((id1,id2))

	return best_aln

## Function to superimpose sequences aligned
def superimpose_structures(structures, best_aln):
	"""
	Superimpose structures.
	Get the chains to create the new structure.
	Save the structure in a new pdb file.
	"""
	print("Perform the superimposition")
	sup = Superimposer() #initialize superimpose
	for aln1, aln2 in best_aln:
		#print("first aln %s, correspond to chain %s" %(aln1, structures[0][aln1]))
		#print("second aln %s, correspond to chain %s" %(aln2, structures[0][aln2]))
		#get fixed and moving atoms for superimposition
		fixed = list(structures[0][aln1].get_atoms())
		moving = list(structures[0][aln2].get_atoms())

		#apply superimposition
		sup.set_atoms(fixed, moving) # set coords
		sup.apply(structures[0][aln2]) # apply to structure-chain objecy
		
		print("RMS(model %s, model %s) = %0.2f" % (aln1, aln2, sup.rms))
		#Delete the fixed chain from structure
		structures[0].detach_child(aln1)


