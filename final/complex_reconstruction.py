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
from common_functions import *

#Define command-line argumets
if __name__ = "__main__":
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
						help = "Visualize the obtained macrocomplex in Chimera")

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
		Handles with different kind of inputs:
		 a list of pdb files
		 a given path with files
		 default input which is the current directory.

		Returns a list of PDB files.
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
		"""Parses the names of the input files and returns the name of the files without .pdb extension"""
		p = re.compile('(.*).pdb')
		m = p.match(filename)

		return m.group(1)


	#Compare chains from two structures
	def chains_comparison(str1, str2):
		"""
		Compares both chains from the input structures.
		Return True if both sequences from one structure are similar to one or both sequences from the other structure.
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
		Extracts the structures and the sequences from PDB input files.
		The function returns three elements:
		 Returns a dictionary with chains of all pairwise interactions
		 Returns a dictionary with homodimers structures
		 Returns a dictionary with heterodimers structures
		"""
		p = PDBParser(PERMISSIVE=1, QUIET=True) #initialize pdb parser without showing warnings

		pairwise_inter = {}
		homodimer_interactions = {}
		heterodimer_interactions = {}

		#get structures for each pdb file
		for f in pdb_files:
			str_id = get_name_structure(f)
			str = p.get_structure(str_id, f)

			#Save chains and sequences for each interaction file
			chains = []
			sequences = []
			for chain in str.get_chains():
				seq = get_sequence(chain)
				if seq is None:
					str[0].detach_child(chain.id) #delte chain from structure
				else:
					sequences.append(seq)
					chains.append(chain)

			#check type of interaction
			if seq_comparison(sequences[0], sequences[1]): #same chains in the interaction files
				#print("Sequences from %s are %s and %s are identical" %(str_id, sequences[0], sequences[1]))
				if homodimer_interactions: #if we have that sequences in other interaction			
					for inter_id, inter in homodimer_interactions.items():
						if chains_comparison(inter, str): #check if the interaction is already save
							#here we need to check if interactions are the same....
							pairwise_inter[inter_id].append((chains[0], chains[1]))
						else:
							pairwise_inter[str_id] = [(chains[0], chains[1])]
							homodimer_interactions[str_id] = str
							break

				else:
					pairwise_inter[str_id] = [(chains[0], chains[1])]
					homodimer_interactions[str_id] = str

			else:
				if heterodimer_interactions:
					for inter_id, inter in heterodimer_interactions.items():
						if chains_comparison(inter, str): #check same interaction from heterodimer chains
							pairwise_inter[inter_id].append((chains[0], chains[1]))
						else:
							pairwise_inter[str_id] = [(chains[0], chains[1])]
							heterodimer_interactions[str_id] = str
							break

				else:
					pairwise_inter[str_id] = [(chains[0], chains[1])]
					heterodimer_interactions[str_id] = str

		return pairwise_inter, homodimer_interactions, heterodimer_interactions

	#Create temporal structure
	def temp_structure(interaction_dic, name):
		"""Builds a new structure and fill it with the interactions. Returns the new built structure"""

		new_structure = Structure.Structure(name)

		i = 0
		for interaction in interaction_dic.values():
			for inter in interaction:
				new_structure.add(Model.Model(i))
				new_structure[i].add(inter[0])
				new_structure[i].add(inter[1])
				i += 1


		return new_structure

	#Save structure
	def save_complex(structure, outfile):
		"""Saves the structure into a PDB file. Uses Bio.PBIO package."""

		print("Printing results in %s" %(outfile))
		io = PDBIO()
		io.set_structure(structure)
		io.save(outfile)

	files = get_input(options.infile)
	print(files)
	pdb_information = get_pdb_info(files)
	interactions = pdb_information[0]
	print(interactions)
	homodimer_interactions = pdb_information[1]
	heterodimer_interactions = pdb_information[2]
	print(homodimer_interactions)
	print(heterodimer_interactions)

	#Control different types of interactions
	if bool(homodimer_interactions) and not bool(heterodimer_interactions):
		print("All interactins are homodimers")
		tmp_str = temp_structure(interactions, "homodimers")
		print(tmp_str)
		new_structure = get_structure_homodimer(tmp_str)
		print(list(new_structure.get_chains()))

		if options.outfile is None:
			options.outfile = "homodimer_complex.pdb"

	elif bool(homodimer_interactions) and bool(heterodimer_interactions):
		print("The pdf file contains both homodimer and heterodimer interactions")
		new_structure = get_structure_homodimer_heterodimer(interactions, homodimer_interactions, heterodimer_interactions)

		if options.outfile is None:
			options.outfile = "homodimer_heterodimer_complex.pdb"

	else:
		print("All interactins are heterodimers")
		tmp_str = temp_structure(interactions, "heterodimers")
		#print(list(tmp_str.get_chains()))
		best_aln = align_sequences_heterodimers(tmp_str)
		print(best_aln)
		new_structure = superimpose_structures_heterodimers(tmp_str, best_aln)
		print(list(new_structure.get_chains()))
		if options.outfile is None:
			options.outfile = "heterodimer_complex.pdb"


	save_complex(new_structure, options.outfile)

	if options.visualize:
		subprocess.Popen('/usr/bin/chimera "options.outfile"', shell= True, executable="/bin/bash")


