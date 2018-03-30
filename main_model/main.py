#####################################################
################### MAIN PROGRAM ####################
#####################################################

#import modules
import re
import sys
import os
import subprocess
import argparse
import string
from Bio.PDB import *
from Bio import pairwise2
import homodimers
import heterodimers

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
					default = None, #"interacion_macrocomplex.pdb",
					help = "Enter the name of output file")

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
		pdb_files = [inputs]

	return pdb_files

#Get the name of the structure
def get_name_structure(filename):
	p = re.compile('(.*).pdb')
	m = p.match(filename)

	return m.group(1)

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

#Compare chain sequences
def seq_comparison(seq1, seq2):
	
	alignment = pairwise2.align.globalxx(seq1, seq2)
	score = alignment[0][2]
	ident_perc = score / max(len(seq1), len(seq2))

	if ident_perc > 0.95:
		return True
	else:
		return False

#Compare chains from two structures
def structure_chain_comparison(str1, str2):
	rs = 0
	rs_tot = [[0,0], [0,0]]
	#print((list(str2.get_chains())))
	for chain1 in str1.get_chains():
		seq1 = get_sequence(chain1)
		rs2 = 0
		for chain2 in str2.get_chains():
			seq2 = get_sequence(chain2)

			print("Comparint %s with %s" %(chain1.get_id(), chain2.get_id()))

			if seq_comparison(seq1, seq2) is True:
				print("The same chain")
				rs_tot[rs][rs2] = 1
			else:
				rs_tot[rs][rs2] = 0
				print("Not The same chain")
			rs2 += 1
		rs += 1

	return rs_tot

#Save structure

def save_complex(structure, outfile):
	print("Saving aligned structure as PDB file %s" % outfile)
	io = PDBIO()
	io.set_structure(structure)
	io.save(outfile)

#Extract information from pdb files
def get_pdb_info(pdb_files):
	"""
	Extract structures from pdb
	Check if we are dealing with homodimers or heterodimers
	"""
	p = PDBParser(PERMISSIVE=1) #initialize pdb parser
	
	pairwise_inter = {}
	uniques_sequences = []
	main_interaction = []

	#get structures for each pdb file
	for f in pdb_files:
		str_id = get_name_structure(f)
		str = p.get_structure(str_id, f)

		chains = []
		sequences = []
		for chain in str.get_chains():
			seq = get_sequence(chain)
			if seq is None:
				str[0].detach_child(chain.id)
			else:
				sequences.append(seq)
				chains.append(chain)
			
		if seq_comparison(sequences[0], sequences[1]): #same chains in the interaction files
			if uniques_sequences and sequences[0] in uniques_sequences: #if we have that sequences in other interaction				
				pairwise_inter["homodimer"].append([chains[0], chains[1]])
					
			else:
				uniques_sequences.append(sequences[0])
				pairwise_inter["homodimer"] = [[chains[0], chains[1]]]
				main_interaction.append(str)

			
		else:
			uniques_sequences.append(sequences[0])
			uniques_sequences.append(sequences[1])
			if "heterodimer" in pairwise_inter.keys():
				pairwise_inter["heterodimer"].extend([chains[0], chains[1]])
			else:
				pairwise_inter["heterodimer"] = ([chains[0], chains[1]])

	if "homodimer" in pairwise_inter.keys():
		new_structure = Structure.Structure('Homodimers')
		n = 0
		for inter in pairwise_inter.values():
			for chains in inter:
				new_structure.add(Model.Model(n))
				new_structure[n].add(chains[0])
				new_structure[n].add(chains[1])
				n += 1

	
	if "heterodimer" in pairwise_inter.keys():
		new_structure = Structure.Structure('Heterodimers')
		new_structure.add(Model.Model(0))
		n = 0
		for chains in pairwise_inter.values():
			for chain in chains:
				try:
					chain.id = alphabet[n]
					new_structure[0].add(chain)
				except:
					new_structure[0].add(chain)	
				n +=1

	return new_structure, main_interaction


files = get_input(options.infile)
structure = get_pdb_info(files)[0]
main_structure = get_pdb_info(files)[1]

if structure.id == "Homodimers":
	print("You have enter %d chains to interact" %(len(list(structure.get_chains()))))
	n_inter = input("How much interactions do you want (enter for current files): ")
	
	if n_inter == "":
		print("Extracting interaction for files")
		new_structure = homodimers.get_structure_homodimer_from_files(structure)
		
	else:
		n = int(n_inter)
		print("Extracting interaction for n interactions")
		#chain_to_add = main_structure[0].copy()
		new_structure = homodimers.get_structure_homodimer_from_number_interaction(main_structure, n)


if structure.id == "Heterodimers":
	best_aln = heterodimers.align_sequences(structure)
	new_structure = heterodimers.superimpose_structures(structure, best_aln)

save_complex(new_structure, options.outfile)




