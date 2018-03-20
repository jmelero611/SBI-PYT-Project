#####################################################
################### MAIN PROGRAM ####################
#####################################################

#import modules
import re
import sys
import os
import argparse

#Define command-line argumets

parser = argparse.ArgumentParser(description="""This program reconstructs the macrocomplex from protein interactions pdb files""")

parser.add_argument('-i', '--input',
					dest = "input",
					action = "store",
					default = os.getcwd(),
					help = "Enter a list of interaction files, a directory or the current directory will be selected")

parser.add_argument('-o', '--output',
					dest = outfile,
					action = store,
					default = sys.stdout,
					help = "Enter the name of output file")

parser.add_argument('-s', '--sequence',
					dest = sequence,
					action = store,
					default = None,
					help = "Enter a sequence of the macrocomplex if avaiable")

parser.add_argument('-v', '--verbose',
					dest = 'verbose',
					action = 'store_true',
					default = False, 
					help = "Print log in stderr")


#Function definition

#Extract input files
def get_input(input):
	"""Handling with different kind of input: a list of pdb files or a given path with files """
	motif = re.compile('.pdb$')
	path = input

	#check if input is a directory
	if os.path.isdir(path):
		pdb_files = [f for f in os.listdir(path) if motif.search(f) is not None]
		os.chdir(path)
	else:
		pdb_files = [input]

	return pdb_files

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

