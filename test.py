#import modules
import re
import sys
import os
import subprocess
import argparse
import string
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

#Get minime distance between two chains
def get_min_distance(chain1, chain2):
	distances = []
	for res1 in chain1.get_residues():
		for res2 in chain1.get_residues():
			distances.append(res1['CA'] - res2['CA'])

	min_distance = min(distances)
	return min_distance


#Compare chain sequence
def seq_comparison(seq1, seq2):
	
	alignment = pairwise2.align.globalxx(seq1, seq2)
	score = alignment[0][2]
	ident_perc = score / max(len(seq1), len(seq2))

	if ident_perc > 0.95:
		return True
	else:
		return False

def get_sequence(chain):
	ppb = CaPPBuilder()
	for residue in chain:
		if residue.id[0] != ' ':
			chain.detach_child(residue.id) #delete heteroatoms
		if len(chain) <= 40:
			return None	#select chains which forms sequences
		else:
			pp = ppb.build_peptides(chain)
			seq = pp[0].get_sequence() #get the
			return seq


#Compare chains from two structures
def structure_chain_comparison(str1, str2):
	rs = 0
	rs_tot = [[0,0], [0,0]]
	#print((list(str2.get_chains())))
	for chain1 in str1.get_chains():
		seq1 = get_sequence(chain1)
		#esto me ralla tambien, porque los resultados son rs_tot[0] para las comparaciones cadena 2 y rs_tot[1] para la cadena1
		#if rs == 0:
		#	rs += 1
		#else:
		#	rs = 0
		rs2 = 0
		for chain2 in str2.get_chains():
			seq2 = get_sequence(chain2)

			print("Comparint %s with %s" %(chain1.get_id(), chain2.get_id()))

			if seq_comparison(seq1, seq2) is True:
				print("The same chain")
				rs_tot[rs][rs2] = 1 #donde hay un 1 es que las cadenas coinciden
			else:
				rs_tot[rs][rs2] = 0
				print("Not The same chain")
			rs2 += 1
		rs += 1 #creo que es esto lo que hay que hacer
	print(rs_tot)
	return rs_tot


#Function to create a dictionary with pdb unique information

def extract_pdb_information(pdb_files, pdb_interaction):

	parser = PDBParser(PERMISSIVE=1)

	for file in pdb_files:
		n = 0
		m = 0
		counter = 0
		pdb_id = file.split('.')[0]
		structure = parser.get_structure(pdb_id, file)
		
		for chain in structure.get_chains(): #Delete secuences with low aa
			if get_sequence(chain) is None:
				structure[0].detach_child(chain.id)

		if pdb_interaction:
			for str in pdb_interaction:
				print("we compare %s with %s" %(pdb_interaction[str], structure))
				tmp_chain_id = []
				results = structure_chain_comparison(pdb_interaction[str], structure)
				print("the results of comparison are %s" %results)

				if 1 in results[0] and 1 in results[1]: #All chains are the same
					counter = 0 #qqui habría que mirar si las interacciones son iguales o no
					#no se si es necesario
				
				#aqui no entiendo bien lo que hace, porque no se que es lo que quiere conseguir
				elif 1 in results[0]:
					if str[0] not in tmp_chain_id:
						tmp_chain_id.append(str[0])
					counter += 1
				elif 1 in results[1]:
					if str[0] not in tmp_chain_id:
						tmp_chain_id.append(str[1])
					counter += 1
				else:
					counter += 1 

				print(tmp_chain_id, counter)	

		else:
			seq = []
			for chain in structure.get_chains():
				seq.append(get_sequence(chain))
			
			if seq_comparison(seq[0], seq[1]) is False:
				m += 1
				chains_id = (alphabet[n], alphabet[m], 0)
			else:
				chains_id = (alphabet[n], alphabet[m], 0)


			pdb_interaction[chains_id] = structure
			
			if chains_id[0] in alphabet:
				alphabet.remove(chains_id[0])
			
			if chains_id[1] in alphabet:
				alphabet.remove(chains_id[1])

#		Esta parte no entiendo bien que ees lo que quiere hacer, si coger aquellas que son iguales o las diferentes
		print(counter, len(pdb_interaction))
		if counter == len(pdb_interaction):
			seq = []
			for chain in structure.get_chains():
				print(chain.id)
				seq.append(get_sequence(chain))

			print(seq)
			print(len(tmp_chain_id))

			if len(tmp_chain_id) == 0:
				if seq_comparison(seq[0], seq[1]) is False:
					m += 1
					chain_id = (alphabet[n], alphabet[m], 0)
				else:
					chain_id = (alphabet[n], alphabet[m], 0)
				print("no tmp chain id")
				print(chain_id)

			if len(tmp_chain_id) == 1:
				if seq_comparison(seq[0], seq[1]) is False:
					chain_id = (tmp_chain_id[0], alphabet[n], 0)
				else:
					chain_id = (tmp_chain_id[0], tmp_chain_id[0], 0)
				print("tmp chain id")
				print(chain_id)

			if len(tmp_chain_id) == 2:
				tmp_chain_id.append(0)
				chain_id = tuple(tmp_chain_id)
			
			if len(tmp_chain_id) == 3:
				chain_id = tuple(tmp_chain_id)

			pdb_interaction[chain_id] = structure

			if tmp_chain_id[0] in alphabet:
				alphabet.remove(tmp_chain_id[0])
            
			#if tmp_chain_id and tmp_chain_id[1] in alphabet:
			#	alphabet.remove(tmp_chain_id[1])

	print(pdb_interaction)
alphabet = list(string.ascii_uppercase) + list(string.ascii_lowercase)
dic = {}
files = get_input(options.infile)
print(files)
print(extract_pdb_information(files, dic))
