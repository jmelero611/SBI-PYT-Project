#####################################################
################### MAIN PROGRAM ####################
#####################################################

#import modules
import re
import sys
import os
import argparse
from Bio.PDB import *
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

#Define command-line argumets

parser = argparse.ArgumentParser(description="""This program reconstructs the macrocomplex from protein interactions pdb files""")

parser.add_argument('-i', '--input',
					dest = "infile",
					action = "store",
					default = os.getcwd(),
					help = "Enter a list of interaction files, a directory or the current directory will be selected")

parser.add_argument('-o', '--output',
					dest = "outfile",
					action = "store",
					default = sys.stdout,
					help = "Enter the name of output file")

parser.add_argument('-s', '--sequence',
					dest = "sequence",
					action = "store",
					default = None,
					help = "Enter a sequence of the macrocomplex if avaiable")

parser.add_argument('-vz', '--visualize',
					dest = "visualize",
					action = "store_true",
					default = False,
					help = "Visualize obtained macrocomplex in Chimera")

parser.add_argument('-v', '--verbose',
					dest = 'verbose',
					action = 'store_true',
					default = False, 
					help = "Print log in stderr")

options = parser.parse_args()

###################### FUNCTIONS ########################

#Extract input files
def get_input(inputs):
	"""Handling with different kind of input: a list of pdb files or a given path with files """
	motif = re.compile('.pdb$')
	path = inputs

	#check if input is a directory
	if os.path.isdir(path):
		pdb_files = [f for f in os.listdir(path) if motif.search(f) is not None]
		os.chdir(path)
	else:
		pdb_files = [input]

	return pdb_files


#Get unique chain interactions from pdb files
def get_pdb_info(pdb_files):
	"""
	Extract all chains of the structures from a list of pdb files
	Consider only unique interactions
	"""
	p = PDBParser(PERMISSIVE=1)#initialize pdb parser
	ppb = CaPPBuilder() #initialize sequence builder by Ca-Ca
	
	structures = {}
	unique_inter = []
	n = 0

	#get structures for each pdb file
	for f in pdb_files:
		n += 1
		str = p.get_structure(n, f)
		str_id = str.get_id()
		structures[str_id] = []
		int = []

		for model in str:
			for chain in model:
				for residue in chain:
					if residue.id[0] != ' ':
						chain.detach_child(residue.id) #delete heteroatoms
				if len(chain) >= 10: #select chains which forms sequences
					pp = ppb.build_peptides(chain)
					seq = pp[0].get_sequence()
					int.append((chain,seq))

		structures[str_id].extend(int)

	#Store first interaction
	ch1 = structures[1][0][0]
	ch2 = structures[1][1][0]
	sq1 = structures[1][0][1]
	sq2 = structures[1][1][1]
	
	unique_chains = [(ch1, ch2)]
	unique_seqs = [(sq1, sq2)]
	
	#print(structures)
	#print(unique_chains)
	if sq1 == sq2:
		print("Homodimers")

	#Iterate the structures to get all unique interactions
	n = 0
	n_iter = 0
	#print(structures)
	for str, inter in structures.items():
		#print(str)
		n +=1
		if n == 1:	continue
		else:

			ch3 = inter[0][0]
			ch4 = inter[1][0]
			sq3 = inter[0][1]
			sq4 = inter[1][1]

			if sq1 != sq3 or sq2 != sq4:
				unique_chains.append((ch3, ch4))
				unique_seqs.append((sq3, sq4))
			else:
				n_iter += 1 #number of interactions of the same type
	
	return unique_chains, unique_seqs, n_iter

## Construct the macrocomplex

def get_macrocomplex(chains, n_iteration):
	"""
	Construct homodimers of ceratins number of proteins.
	"""
	io = PDBIO()
	sup = Superimposer()
	file = "output_homodimer.pdb"
	base_struct = Structure.Structure('0')
	base_struct.add(Model.Model(0))

	main_chain = chains[0][0]
	fix = main_chain.copy()
	complex = []
	complex.append(fix)
	main_id = fix.id
	#base_struct[0].add(fix)
	#print(fix.get_full_id())
	#print(list(base_struct.get_chains()))
	#fix2 = list(inter[1].get_atoms())

	dic = ['B', 'C', 'D', 'E', 'F', 'G', 'H', 'I','J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
	
	for i in range(1, n_iteration):
		#base_struct.id = str(n)
		mov = main_chain.copy()
		mov.id = dic[i]
		#if base_struct[0].has_id(mov.id):
		#	print("mierda")
		#else:
		#	print("ye")
		#print(mov.get_full_id())
		sup.set_atoms(list(fix.get_atoms()), list(mov.get_atoms()))
		sup.apply(mov)
		#print(complex)
		#base_struct.add(mov)
		complex.append(mov.copy())
		#print(list(base_struct[0].get_chains()))

		#fix = mov
		#print(mov.id)
		#print("\n")

		#mov.id = ""

	#i = 0
	#print(complex)
	
	for chain in complex:
		#id_c = chain.get_id()
		#print(i)
		#chain.id = "%s%d" %(id_c, i)
		#print(chain.id)
		base_struct[0].add(chain)
		#i += 1

	print(list(base_struct.get_chains()))
	#for r in residues:
		#print(type(r.get_id()[1]))
	#for a in atoms:
	#	print(a.get_fullname())
	
	io.set_structure(base_struct)
	#ch = list(base_struct.get_chains())
	#for c in ch:
		#print(c.get_id())
	io.save(file)


files = get_input(options.infile)
info = get_pdb_info(files)
chains = info[0]
n = info[2]
get_macrocomplex(chains,n)

#Extract sequence
def get_sequences(seqs):
	"""Get the sequence of the interactions if avaiable"""
	motif = re.compile('.fa$|.fasta$')
	seq_path = seqs

	if os.path.isdir(seq_path):
		fasta_files = [f for f in os.listdir(seq_path) if motif.search(f) is not None]
	else:
		fasta_files = [seqs]

	return fasta_files

