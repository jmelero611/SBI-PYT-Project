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

###################### COMMAND-LINE ARGUMENTS ########################

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

########################## FUNCTIONS #############################

#Extract input files
def get_input(inputs):
	"""
	Handles with different kind of inputs:
	 a list of pdb files, separed by a comma
	 a given path with pdb files
	 default input which is the current directory.

	 output: a list of PDB files.
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
			files = inputs.split(',')
			pdb_files = [f for f in files]
		else:
			pdb_files = [inputs]

	return pdb_files

#Get the name of the structure
def get_name_structure(filename):
	"""Parses the names of the input files and returns the name of the files without .pdb extension"""
	p = re.compile('(.*).pdb')
	m = p.match(filename)

	return m.group(1)


#Extract information from pdb files
def get_pdb_info(pdb_files):
	"""
	Extracts the structures and the chains from PDB input files.
	The function returns three outputs:
	 a dictionary with chains of all pairwise interactions
	 a dictionary with homodimers structures
	 a dictionary with heterodimers structures
	 a list with files with errors.
	"""
	p = PDBParser(PERMISSIVE=1, QUIET=True)
	
	pairwise_inter = {}
	homodimer_interactions = {}
	heterodimer_interactions = {}
	error_files = []
	#get structures for each pdb file
	for f in pdb_files:
		try:						
			str_id = get_name_structure(f)
			str = p.get_structure(str_id, f)	
		except:
			error_files.append(f)
			continue

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

		#Check input file is a pairwise interaction
		if len(chains) > 2 or len(chains) == 1 or len(list(str.get_models())) > 1:
			error_files.append(f)
			continue

		#check type of interaction
		if seq_comparison(sequences[0], sequences[1]): #same chains in the interaction file
			if homodimer_interactions:			
				for inter_id, inter in homodimer_interactions.items():
					if chains_comparison(inter, str): #check if we have already that homodimer interaction
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
					if chains_comparison(inter, str): #check if we have already that heterodimer interaction
						pairwise_inter[inter_id].append((chains[0], chains[1]))
					else:
						pairwise_inter[str_id] = [(chains[0], chains[1])]
						heterodimer_interactions[str_id] = str
						break
						
			else:
				pairwise_inter[str_id] = [(chains[0], chains[1])]
				heterodimer_interactions[str_id] = str
	
	return pairwise_inter, homodimer_interactions, heterodimer_interactions, error_files

########################## MAIN PROGRAM #############################

if __name__ == "__main__":

	files = get_input(options.infile)
	pdb_information = get_pdb_info(files)

	interactions = pdb_information[0]
	homodimer_interactions = pdb_information[1]
	heterodimer_interactions = pdb_information[2]
	wrong_files = pdb_information[3]

	if options.verbose:
		if wrong_files:
			sys.stderr.write("The following %i files were not analizing: \n" %(len(wrong_files)))
			for file in wrong_files:
				if "pdb" in file:
					sys.stderr.write("\tThe file %s because it is not a pairwise interaction file.\n" %(file))
				else:
					sys.stderr.write("\tThe file %s because it is not a pdb file\n" %(file))
		
		sys.stderr.write("The following %i files are going to be analyzed: \n" %(len(files) - len(wrong_files)))
		for file in files:
			if file not in wrong_files:
				sys.stderr.write("\t- %s\n" %(file))

	#Control different types of interactions
	#All homodimers
	if bool(homodimer_interactions) and not bool(heterodimer_interactions):
		if options.verbose:
			sys.stderr.write("All pdb files contains %i homodimer interactions" %(len(homodimer_interactions)))

		tmp_str = temp_structure(interactions, "homodimers")
		new_structure = get_structure_homodimer(tmp_str)

		if options.outfile is None:
			options.outfile = "homodimer_complex.pdb"

	#All heterodimer
	elif bool(heterodimer_interactions) and not bool(homodimer_interactions):
		if options.verbose:
			sys.stderr.write("All pdb files contains %i heterodimer interactions" %(len(heterodimer_interactions)))
	
		tmp_str = temp_structure(interactions, "heterodimers")
		best_aln = align_sequences_heterodimers(tmp_str, options.verbose)
		new_structure = superimpose_structures_heterodimers(tmp_str, best_aln, options.verbose)
		
		if options.outfile is None:
			options.outfile = "heterodimer_complex.pdb"

	#Both, homodimer and heterodimer
	elif bool(homodimer_interactions) and bool(heterodimer_interactions):
		if options.verbose:
			sys.stderr.write("The pdb file contains %i homodimer and %i heterodimer interactions" %(len(homodimer_interactions), len(heterodimer_interactions)))
	
		new_structure = get_structure_homodimer_heterodimer(interactions, homodimer_interactions, heterodimer_interactions)

		if options.outfile is None:
			options.outfile = "homodimer_heterodimer_complex.pdb"

	else:
		sys.stderr.write("All the files introduced are not valid. Please check if there are pairwise interaction PDB files")
		exit()

	#Save the structure in a pdb file
	if options.verbose:
		sys.stderr.write("Printing the results in the PDB file: %s" %(options.outfile))

	save_complex(new_structure, options.outfile)

	if options.visualize:
		subprocess.Popen('/usr/bin/chimera "options.outfile"', shell= True, executable="/bin/bash")


		
