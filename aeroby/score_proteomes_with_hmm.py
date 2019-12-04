from __future__ import print_function
import os, re, sys, glob
from collections import defaultdict
#search a HMM against a set of proteomes and print out a table with 0/1 for hits above a certain threshold. Maybe do this for multiple hmms and print out one table.

def map_id_to_species_name(proteome_name):
	id_to_return = ''
	for id in mapping.keys():
		if proteome_name.startswith(id):
			if id_to_return == '':
				id_to_return = mapping[id]
			else:
				print("Uh oh! Duplicates in the mapping..." + id + "\t" + mapping[id])
	if id_to_return == '':
		id_to_return = 'None'
	return id_to_return
	

#create mapping between proteome names and species taxonomy string
mapping = {}
inh = open("arc_taxonomy_r86.tsv")
for line in inh:
	fields = re.split("\t", line.rstrip())
	mapping[fields[0][3:]] = fields[1]
inh.close()

results = defaultdict(dict)
hmms = glob.glob("hmms/*hmm")
proteomes = glob.glob("proteomes/*fas")
for hmm in hmms:
	for proteome in proteomes:
		hmm_outfile = os.path.basename(hmm)[:-4] + "_" + os.path.basename(proteome)[:-4] + "_tbl.out"
		os.system("hmmsearch --tblout " + hmm_outfile + " -E 0.0000001 " + hmm + " " + proteome + " > /dev/null 2>&1")
		#check if there was a good enough hit
		inh = open(hmm_outfile)
		good_hit = 0
		for line in inh:
			if line.startswith("#"):
				continue
			else:
				good_hit = 1
				print(line.rstrip())
		if good_hit == 1:
			results[proteome][hmm] = 1
		else:
			results[proteome][hmm] = 0

#print out results structure
for proteome in results:
	species_string = map_id_to_species_name(os.path.basename(proteome))
	toprint = os.path.basename(proteome) + "\t" + species_string + "\t"
	for hmm in results[proteome]:
			toprint = toprint + "\t" + str(results[proteome][hmm])
	print(toprint)
