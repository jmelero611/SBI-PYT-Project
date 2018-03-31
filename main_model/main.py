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
from homodimers import *
from heterodimers import *

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
					default = None, 
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
	"""
	Handling with different kind of input:
	 a list of pdb files
	 a given path with files
	 default input which is the current directory.
	 """
	motif = re.compile('.pdb$')
	path = inputs

	#check if input is a directory
	if os.path.isdir(path):
		pdb_files = [f for f in os.listdir(path) if motif.search(f) is not None]
		os.chdir(path)
	else:
		#check if input is a list
		if ',' in inputs:
			pdb_files = inputs.split(',')
		else:
			pdb_files = [inputs]

	return pdb_files

#Get the name of the structure
def get_name_structure(filename):
	p = re.compile('(.*).pdb')
	m = p.match(filename)

	return m.group(1)

#Get sequence from chain taking into account alpha carbons
def get_sequence(chain):
	"""
	Get sequence from chain with its alpha carbons
	""" 
	ppb = CaPPBuilder()

	for residue in chain:
		if residue.id[0] != ' ':
			chain.detach_child(residue.id) #delete heteroatoms
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

	if ident_perc > 0.95:
		return True
	else:
		return False

#Compare chains from two structures
def chains_comparison(str1, str2):
	"""
	Compare both chains from each structure.
	Return True if both sequences from str1 are similar to one or both of sequences from str
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


#Extract information from pdb files
def get_pdb_info(pdb_files):
	"""
	Extract structures from pdb
	Check if we are dealing with homodimers or heterodimers
	"""
	p = PDBParser(PERMISSIVE=1, QUIET=True) #initialize pdb parser without showing warnings
	
	pairwise_inter = {}
	homodimer_interactions = {}
	heterodimer_interactions = {}
	homodimers = False

	#get structures for each pdb file
	for f in pdb_files:
		str_id = get_name_structure(f)
		str = p.get_structure(str_id, f)

		chains = []
		sequences = []
		for chain in str.get_chains():
			seq = get_sequence(chain)
			if seq is None:
				str[0].detach_child(chain.id) #delte chain from structure
			else:
				sequences.append(seq)
				chains.append(chain)
		
		if seq_comparison(sequences[0], sequences[1]): #same chains in the interaction files
			homodimers = True
			if homodimer_interactions: #if we have that sequences in other interaction			
				for inter_id, inter in homodimer_interactions.items():
					if chains_comparison(inter, str): #check if it the same interaction as befor
						#here we need to check if interactions are the same....
						pairwise_inter[inter_id].append([chains[0], chains[1]])
					else:
						pairwise_inter[str_id] = [[chains[0], chains[1]]]
						homodimer_interactions[str_id] = str
					
			else:
				pairwise_inter[str_id] = [[chains[0], chains[1]]] #first interaction
				homodimer_interactions[str_id] = str
		
		else:
			if heterodimer_interactions:
				for inter_id, inter in heterodimer_interactions.items():
					if chains_comparison(inter, str): #check same interaction from heterodimer chains
						homodimers = True
						pairwise_inter[inter_id].append([chains[0], chains[1]])
					else:
						pairwise_inter[str_id] = [chains[0], chains[1]]
						
			else:
				pairwise_inter[str_id] = [chains[0], chains[1]]
				heterodimer_interactions[str_id] = str
	
	return pairwise_inter, homodimers

#Create temporal structure

def temp_structure(interaction_dic, name):
	"""
	Build a structure with interactions 
	"""
	new_structure = Structure.Structure(name)
	n = 0
	i = 0
	for inter in interaction_dic.values():
		if len(inter) > 2:
			for chains in inter:
				new_structure.add(Model.Model(i))
				new_structure[i].add(chains[0])
				new_structure[i].add(chains[1])	
				i += 1
		else:
			new_structure.add(Model.Model(n))
			new_structure[n].add(inter[0])
			new_structure[n].add(inter[1])	
		n += 1

	return new_structure

#Save structure

def save_complex(structure, outfile):
	"""
	Save the new structure in output pdb file.
	"""
	print("Printing results in %s" %(outfile))
	io = PDBIO()
	io.set_structure(structure)
	io.save(outfile)

files = get_input(options.infile)
interactions = get_pdb_info(files)[0]
homodimers = get_pdb_info(files)[1]


if homodimers:
	#print("You have enter %d chains to interact" %(len(list(structure.get_chains()))))
	n_inter = ""#input("How much interactions do you want (enter for current files): ")
	tmp_str = temp_structure(interactions, "homodimers")
	
	if n_inter == "":
		print("Extracting interaction for files")
		new_structure = get_structure_homodimer_from_files(tmp_str)
		
	else:
		n = int(n_inter)
		print("Extracting interaction for n interactions")
		#chain_to_add = main_structure[0].copy()
		new_structure = get_structure_homodimer_from_number_interaction(main_structure, n)

	if options.outfile is None:
		options.outfile = "homodimer.pdb"

else:
	tmp_str = temp_structure(interactions, "heterodimers")
	best_aln = align_sequences_heterodimers(tmp_str)
	print(best_aln)
	new_structure = superimpose_structures_heterodimers(tmp_str, best_aln)

	if options.outfile is None:
		options.outfile = "heterodimer.pdb"


save_complex(new_structure, options.outfile)

if visualize:
	subprocess.Popen('/usr/bin/chimera "options.outfile"', shell= True, executable="/bin/bash")




