#given a HMM and a target proteome, find the ortholog (if any) by the cognitor-ish method
from __future__ import print_function
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
import os, re, sys

hmmfile = sys.argv[1]
target_proteome = sys.argv[2]

#which ortholog group are we looking for?
hmm_bits = re.split("_", os.path.basename(hmmfile))
ortho_group = hmm_bits[0]
reference_mapping_file = "marker_map_ecoli.txt" #this should be domain-appropriate
refmap = {}
inh = open(reference_mapping_file)
for line in inh:
	fields = re.split("\t", line.rstrip())
	refmap[fields[1]] = fields[0]

target_proteins = SeqIO.index(target_proteome, "fasta")
target_protein_ids = []
ortholog = 'None'

for rec in target_proteins:
	target_protein_ids.append(rec)

tag = os.path.basename(hmmfile) + "_" + os.path.basename(target_proteome)
#search the HMM against the target proteome
os.system("hmmsearch --tblout " + tag + ".tbl" + " " + hmmfile + " " + target_proteome + " > /dev/null")
#extract the top hit
inh = open(tag + ".tbl")
first_hit = 1
for line in inh:
	if line.startswith("#"):
		continue
	else:
		if first_hit == 1: #take this one
			first_hit = 0
			fields = re.split("\s+", line.rstrip())
			best_hit = fields[0].rstrip()
			hit_seq = ''
			for id in target_protein_ids:
				if id.startswith(best_hit):
					hit_seq = str(target_proteins[best_hit].seq)
			#do more...	
			target_hit_seqfile = tag + "_back.seq"
			outh = open(target_hit_seqfile, "w")
			outh.write(">" + best_hit + "\n" + hit_seq + "\n")	
			outh.close()
			
			blastx_cline = NcbiblastxCommandline(cmd='blastp', query=target_hit_seqfile, db="GCF_000005845.2_ASM584v2_protein.faa", evalue=0.001, outfmt=5, out=target_hit_seqfile + ".search")
			stdout, stderr = blastx_cline()
			result_handle = open(target_hit_seqfile + ".search")
			blast_record = NCBIXML.read(result_handle)
			if len(blast_record.descriptions) > 0: 
				top_hit = blast_record.descriptions[0]
				th_fields = re.split(" ", str(top_hit))
				genome_name = os.path.basename(target_proteome)
				if th_fields[1] in refmap:
					if refmap[th_fields[1]] == ortho_group:
						print(ortho_group + "\t" + genome_name + "\t" + best_hit)
					else:
						print(ortho_group + "\t" + genome_name + "\t" + "None")
				else:
					print(ortho_group + "\t" + genome_name + "\t" + "None")
