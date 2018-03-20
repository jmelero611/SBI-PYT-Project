##PYTHON SBI PROJECT

########################################################
######################TOY MODEL ########################
########################################################

#Simple input file with iterations between heterodimers
from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import sys

###################### FUNCTIONS ########################

##Function to extract chains from pdb files
def get_pdb_info(pdb_files):
	"""
	Extract all chains of the structures from a list of pdb files
	"""
	p = PDBParser(PERMISSIVE=1)#initialize pdb parser

	structures = []
	chains = []
	n = 0
	#get structures for each pdb file
	for f in pdb_files:
		n += 1
		str = p.get_structure(n, f)
		structures.append(str)

	for str in structures:
		residue_to_remove = []
		chain_to_remove = []
		for model in str:
			for chain in model:
				for residue in chain:
					if residue.id[0] != ' ':
						chain.detach_child(residue.id)
				if len(chain) >= 10:
					chains.append(chain)
	return chains

##Function to align sequences from chains
def align_sequences(dic):
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
	ppb = CaPPBuilder() #create sequence with Ca-Ca 
	n = len(dic)
	matrix = matlist.blosum62
	gap_open = -10.0
	gap_extend = -0.5
	best_aln = []
	all_aln = []

	#Select sequences to align
	for i in range(n-1):
		ch1 = dic[i].get_id()
		pp1 = ppb.build_peptides(dic[i])
		
		if pp1 != []:
			id1 = "%s-%d" %(ch1, i)
			seq1 = pp1[0].get_sequence()
		else:
			continue
		
		for j in range(i+1,n):
			ch2 = dic[j].get_id()
			pp2 = ppb.build_peptides(dic[j])
			
			if pp2 != []:
				id2 = "%s-%d" %(ch2, j)
				seq2 = pp2[0].get_sequence()
			else:
				continue


			#Align sequence
			alns = pairwise2.align.globalds(seq1, seq2, 
											matrix, gap_open, gap_extend,
											penalize_end_gaps = (False, False))
			top_aln = alns[0]
			al1, al2, score, begin,end = top_aln

			#Calculate sequence identity
			seq_id, g_seq_id = calculate_identity(al1, al2)

			all_aln.append((id1, id2, seq_id))

			if seq_id > 99:
				best_aln.append((id1,id2))

	return best_aln, all_aln


## Function to get the atoms to superimpose
def get_superimpose_atoms(chains, best_aln):
	"""
	Get the atoms list for the chains that had aligned.
	Get the chains object for the new structure
	"""
	n = 0
	fixed = []
	moving = []
	id_moving = []
	id_fixed = []

	for aln in best_aln:
		fix_id = aln[0]
		mov_id = aln[1]
		#print((fix_id, mov_id))
		for chain in chains:
			id_c = "%s-%d" %(chain.get_id(), n)
			
			if id_c == fix_id:
				fixed.append(list(chain.get_atoms()))
				id_fixed.append(n)

			if id_c == mov_id:
				moving.append(list(chain.get_atoms()))
				id_moving.append(n)

			n += 1

	return fixed, moving, id_moving, id_fixed


## Function to superimpose sequences aligned
def superimpose_atoms(fixed_atoms, moving_atoms, ids_moving, ids_fixed, chains):
	"""
	Superimpose structures.
	Get the chains to create the new structure.
	"""
	sup = Superimposer() #initialize superimpose
	
	for i in range(len(fixed_atoms)):
		sup.set_atoms(fixed_atoms[i], moving_atoms[i]) # set coords
		sup.apply(chains[ids_moving[i]]) # apply to chains objects
		del chains[ids_fixed[i]]

## Function to create the new structure

def get_structure(chains):
	"""
	Create a new structure object from the selected chains.
	"""
	io = PDBIO()
	for chain in chains:
		ch_id = chain.get_id()
		io.set_structure(chain)
		io.save("chain_%s.pdb" %(ch_id))


###################### MAIN SCRIPT ########################

#get files
files = []

for f in sys.argv[1:]:
	files.append(f)

#get chains
chains = get_pdb_info(files)

#align seqs
alignment = align_sequences(chains)
best_aln = alignment[0]
all_aln = alignment[1]
	
#get superimpose atoms
atoms = get_superimpose_atoms(chains, best_aln)
fixed = atoms[0]
moving = atoms[1]
id_moving = atoms[2]
id_fixed = atoms[3]

#superimposition
superimpose_atoms(fixed, moving, id_moving, id_fixed, chains)
	
#create output
get_structure(chains)
