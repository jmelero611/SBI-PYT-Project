#import modules
import re
import sys
import os
import subprocess
import argparse
import string
import numpy as np
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


#Compare chain sequences
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

#get numeric array from the alignment
def get_numeric_array(seq_alignment):
	array = []
	n = 0
	for char in seq_alignment:
		if char == "-":
			array.append('-')
		else:
			array.append(n)
		n += 1
		
	return array

#Compare chains from two structures
def structure_chain_comparison(str1, str2):
	rs = 0
	rs_tot = [[0,0], [0,0]]
	#print((list(str2.get_chains())))
	for chain1 in str1.get_chains():
		seq1 = get_sequence(chain1)

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
				rs_tot[rs][rs2] = 1
			else:
				rs_tot[rs][rs2] = 0
				print("Not The same chain")
			rs2 += 1
		rs += 1

	return rs_tot

#Compare two interactions
def compare_interactions(inter1, inter2):
	str1 = Structure.Structure('1')
	str2 = Structure.Structure('2')

	str1.add(Model.Model(0))
	str2.add(Model.Model(0))

	homodimer = False
	for chain in inter1:
		if len(list(str1[0].get_chains())) == 1 and seq_comparison(chain, list(str1[0].get_chains()[0])):
			homodimer = True

		str1[0].add(Chain.Chain(chain.get_id()))
		res_count = 0
		for residue in chain:
			for atom in residue:
				if atom.get_id() == 'CA':
					str1[0][chain.get_id()].add(Residue.Residue(('',res_count,''),residue.get_resname(),residue.get_segis()))

					str1[0][chain.get_id()][('',res_count,'')].add(atom.copy())
					res_count += 1

	for chain in inter2:
		str2[0].add(Chain.Chain(chain.get_id()))
		res_count = 0
		for residue in chain:
			for atom in residue:
				if atom.get_id() == 'CA':
					str2[0][chain.get_id()].add(Residue.Residue(('',res_count,''),residue.get_resname(),residue.get_segis()))

					str2[0][chain.get_id()][('',res_count,'')].add(atom.copy())
					res_count += 1

	print(list(str1.get_chains()))
	print(list(str2.get_chains()))
	print(homodimer)


#Compare if two structures are structurally similar
def compare_structures_superimpose(str1, str2):
	chains1 = [ch for ch in str1.get_chains()]
	chains2 = [ch for ch in str2.get_chains()]

	sup = Superimposer()
	mean_dist = []

	for chain1 in chains1:
		for chain2 in chains2:
			seq1 = get_sequence(chain1)
			seq2 = get_sequence(chain2)

			if not seq_comparison(seq1, seq2):
				continue

			for round in range(10):
				sup.set_atoms(list(str1[0][chain1.get_id()].get_atoms()), list(str2[0][chain1.get_id()].get_atoms()))
				sup.apply(str2)


			distances = []
			chain3 = [ch.get_id() for ch in chains1 if ch != chain1.get_id()][0]
			chain4 = [ch.get_id() for ch in chains2 if ch != chain2.get_id()][0]

			CA3 = [res['CA'] for res in str1[0][chain3].get_residues() if 'CA' in [atm.get_id() for atm in res.get_atoms()]]
			CA4 = [res['CA'] for res in str2[0][chain4].get_residues() if 'CA' in [atm.get_id() for atm in res.get_atoms()]]

			for atoms in zip(CA3, CA4):
				atom1 = atoms[0].get_coord()
				atom2 = atoms[1].get_coord()
				diff = atom1 - atom2
				diff_sq = list(map(lambda x: pow(x,2), diff))
				dist = np.sqrt(sum(diff_sq))
				distances.append(dist)

			mean_dist.append(sum(distances)/len(distances))

			if min(mean_dist) < 9:
				return 0
	return 1
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
				#print("we compare %s with %s" %(pdb_interaction[str], structure))
				tmp_chain_id = []
				results = structure_chain_comparison(pdb_interaction[str], structure)
				print("the results of comparison are %s" %results)

				if 1 in results[0] and 1 in results[1]: #All chains are the same
					str_comp = compare_structures_superimpose(pdb_interaction[str], structure)
					#0 if structures are similar ans 1 if not
					counter += str_comp

					if str_comp == 1: #different structures
						n = 0
						repeat = 0

						for key in [id for id in pdb_interaction if list(id[:2]) == [str[0], str[1]]]:
							repeat.append(compare_structures_superimpose(pdb_interaction[key], structure))
							n.append(key[2])

						if sum(repeat) == len(pdb_interaction):
							tmp_chain_id = [str[0], str[1], max(num) +1]
							counter = len(pdb_interaction)
							break
						else:
							counter = 0
						break
					else:
						break


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

		
		if counter == len(pdb_interaction):
			seq = []

			for chain in structure.get_chains():
				seq.append(get_sequence(chain))

			if len(tmp_chain_id) == 0:
				if seq_comparison(seq[0], seq[1]) is False: #<0,95%
					m += 1
					chain_id = (alphabet[n], alphabet[m], 0)
				else:
					chain_id = (alphabet[n], alphabet[m], 0)

			if len(tmp_chain_id) == 1:

				if seq_comparison(seq[0], seq[1]) is False:
					chain_id = tuple([tmp_chain_id[0], alphabet[n], 0])
					tmp_chain_id.append(alphabet[n])
				else:
					chain_id = tuple([tmp_chain_id[0], tmp_chain_id[0], 0])
					tmp_chain_id.append(tmp_chain_id[0])

			if len(tmp_chain_id) == 2:
				tmp_chain_id.append(0)
				chain_id = tuple(tmp_chain_id)
			
			if len(tmp_chain_id) == 3:
				chain_id = tuple(tmp_chain_id)

			pdb_interaction[chain_id] = structure


			if tmp_chain_id and tmp_chain_id[0] in alphabet:
				alphabet.remove(tmp_chain_id[0])
            
			if tmp_chain_id and tmp_chain_id[1] in alphabet:
				alphabet.remove(tmp_chain_id[1])

	return pdb_interaction

def macrocomplex_builder(str, sup_dic = {}):
	for interact in dic.values():
		sup = Superimposer()
		sup.set_atoms(list(str.get_atoms()), list(interact.get_atoms()))
		if not np.abs(sup.rms) < 10:
			macrocomplex_builder(sup)
		else:
			sup_dic = sup

	print(sup_dic)
alphabet = list(string.ascii_uppercase) + list(string.ascii_lowercase)
dic = {}
files = get_input(options.infile)
extract_pdb_information(files, dic)

print(dic)


