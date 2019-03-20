from __future__ import print_function
from Bio import SeqIO
from collections import defaultdict
import os, sys, re, glob

#take a set of orthologues, and filter them using the output from a phylter (Damien de Vienne) analysis.
to_delete = []
to_edit = defaultdict(list) #key = file to edit, values = list of seqs (cell outliers) to remove

phylter_out = open("OUTLIERS-REMOVED.txt")
delete_mode = 0
edit_mode = 0
for line in phylter_out:
	if line.startswith("## Complete outlier"):
		delete_mode = 1
		continue
	elif line.startswith("## Cell outliers"):
		edit_mode = 1
		continue
	elif re.match("^$", line.rstrip()):
		delete_mode = 0
		edit_mode = 0
	else:
		if delete_mode == 1:
			name_bits = re.split("\.", line.rstrip())
			to_delete.append(name_bits[0])
		elif edit_mode == 1:
			fields = re.split(" ", line.rstrip())
			name_bits = re.split("\.", fields[0])
			to_edit[name_bits[0]].append(fields[1])

print(to_delete)
print(to_edit)

#now edit the alignment files accordingly.

input_alignments = glob.glob("*aln")
if os.path.exists("edited"):
	pass	
else:
	os.system("mkdir edited")
for align in input_alignments:
	name_bits = re.split("\.", align)
	root = name_bits[0]
	if root in to_delete:
		continue
	if root in to_edit:
		seqs = SeqIO.index(align, "fasta")
		outh = open("edited/" + align, "w")	
		for rec in seqs:
			if seqs[rec].description in to_edit[root]:
				continue
			else:
				outh.write(">" + seqs[rec].description + "\n" + str(seqs[rec].seq) + "\n")
		outh.close()
	else:
		os.system("cp " + align + " edited/")
