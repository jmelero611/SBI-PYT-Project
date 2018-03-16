##PYTHON SBI PROJECT

###########################
########TOY MODEL #########
###########################

#Simple input file with iterations between heterodimers
from Bio.PDB import *
from Bio import pairwise2
import sys

#get files
files = []

for f in sys.argv[1:]:
	files.append(f)

#extract information for pdb files: atoms and sequences for each chain
p = PDBParser(PERMISSIVE=1)
ppb = CaPPBuilder() #create sequence with Ca-Ca 

structures = []
n = 0

#get the structure for each file
for f in files:
	n += 1
	str = p.get_structure(n, f)
	structures.append(str)

pdb_info = {}
for struct in structures:
	for model in struct:
		for chain in model:
			atom = chain.get_atoms()
			pp = ppb.build_peptides(chain)
			if pp != []:
				seq = pp[0].get_sequence()
				ch = chain.get_id()
				fi = struct.get_id()
				pdb_info[(ch,fi)] = {}
				pdb_info[(ch,fi)]["atoms"] = atom
				pdb_info[(ch,fi)]["sequence"] = seq

#align sequences
n = len(pdb_info)
max_score = 0
best_alignment = None #should be more than one
alignments = {}
ids = list(pdb_info.keys())

for i in range(n-1):
	id1 = ids[i]
	seq1 = pdb_info[id1]['sequence'] 
	for j in range(i+1,n):
		id2 = ids[j]
		seq2 = pdb_info[id2]['sequence']

		alns = pairwise2.align.globalxx(seq1, seq2)
		top_aln = alns[0]
		al1, al2, score, begin,end = top_aln
		alignments[(id1,id2)] = top_aln
		if score > max_score: #we should use a cut-off
			max_score = score
			best_alignment = (id1,id2)
			
#Superimpose aligned sequences
sup = Superimposer()

#get the list of atoms to perform the superimposer
fixed = list(pdb_info[best_alignment[0]]['atoms'])
moving = list(pdb_info[best_alignment[1]]['atoms'])
id_m = best_alignment[1][1] #get id for moving structure

sup.set_atoms(fixed, moving)

sup.apply(structures[id_m - 1].get_atoms()) #change coordinates for moving structure

#get new structure
pdb_out = "toy_model.pdb"

io = PDBIO()

#it only save the last structure
for s in structures:
	io.set_structure(s)
	io.save(pdb_out)
