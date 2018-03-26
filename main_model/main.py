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
		pdb_files = [inputs]

	return pdb_files


#Extract information from pdb files
def get_pdb_info(pdb_files):
	"""
	Extract structures from pdb
	Check if we are dealing with homodimers or heterodimers
	"""
	p = PDBParser(PERMISSIVE=1) #initialize pdb parser
	ppb = CaPPBuilder() #initialize sequence builder by Ca-Ca
	
	structures = Structure.Structure('Macrocomplex')
	uniques_sequences = []
	n = 0

	#get structures for each pdb file
	for f in pdb_files:
		str = p.get_structure(n, f)
		
		for model in str:
			if n != 0:
				model.id = n
			for chain in model:
				for residue in chain:
					if residue.id[0] != ' ':
						chain.detach_child(residue.id) #delete heteroatoms
				#print((chain.get_full_id(), len(chain)))
				if len(chain) <= 20:
					#select chains which forms sequences
					model.detach_child(chain.id)
				else:
					pp = ppb.build_peptides(chain)
					seq = pp[0].get_sequence() #get the sequence of the chain
					if seq not in uniques_sequences:
						uniques_sequences.append(seq) #add unique sequences to know the type of complex
			
			structures.add(model)
		n += 1

	print(uniques_sequences)
	#Analyse type of dimer
	if len(uniques_sequences) == 1:	
		type_inter = "Homodimer"
		print("Handling with homodimers")
	elif len(uniques_sequences) == 2:
		type_inter = "Homodimer"
		print("Handling with heterodimers of two chains")
	else:
		type_inter = "Heterodimer"

	return structures, type_inter


files = get_input(options.infile)
info = get_pdb_info(files)
structures = info[0]
interaction = info[1]

if interaction == "Homodimer":
	homodimers.get_structure_homodimer(structures, options.outfile)

if interaction == "Heterodimer":
	best_aln = heterodimers.align_sequences(structures)
	heterodimers.superimpose_structures(structures, best_aln, options.outfile)


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

