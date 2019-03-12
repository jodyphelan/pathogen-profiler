from .utils import *
import json
import re

def get_missense_codon(x):
	re_obj = re.search("([0-9]+)",x)
	if re_obj:
		return int(re_obj.group(1))
	else:
		log("Error can't find codon number in %s" % x,True)

def get_indel_nucleotide(x):
	re_obj = re.search("([0-9]+)",x)
	if re_obj:
		return int(re_obj.group(1))
	else:
		log("Error can't find nucleotide number in %s" % x,True)

def barcode(mutations,barcode_bed):
	bed_num_col = len(open(barcode_bed).readline().rstrip().split())
	cols = [1,3,4,5,6]+list(range(7,bed_num_col+1)) if bed_num_col>6 else [1,3,4,5,6]
	bed = load_bed(barcode_bed,cols,1,3,intasint=True)
	add_info = load_bed(barcode_bed,cols,4)
	#{'Chromosome':{'4392120': ('Chromosome', '4392120', 'lineage4.4.1.2', 'G', 'A', 'Euro-American', 'T1', 'None')}}

	barcode_support = defaultdict(list)
	for chrom in bed:
		for pos in bed[chrom]:
			marker = bed[chrom][pos]
			tmp = [0,0]
			if chrom in mutations and pos in mutations[chrom]:
				if marker[3] in mutations[chrom][pos]: tmp[0] = mutations[chrom][pos][marker[3]]
				if marker[4] in mutations[chrom][pos]: tmp[1] = mutations[chrom][pos][marker[4]]
			if  tmp==[0,0]: continue
			barcode_support[marker[2]].append(tmp)
	barcode_frac = defaultdict(float)
	for l in barcode_support:
		if stdev([x[1]/(x[0]+x[1]) for x in barcode_support[l]])>0.15: continue
		barcode_pos_reads = sum([x[1] for x in barcode_support[l]])
		barcode_neg_reads = sum([x[0] for x in barcode_support[l]])
		lf = barcode_pos_reads/(barcode_pos_reads+barcode_neg_reads)
		if lf<0.05:continue
		barcode_frac[l] = lf
	final_results = []

	for l in barcode_frac:
		tmp = {"annotation":l,"freq":barcode_frac[l],"info":[]}
		if bed_num_col>6:
			for i in range(5,bed_num_col-1):
				tmp["info"].append(add_info[l][i])
		final_results.append(tmp)
	return final_results

def db_compare(mutations,db_file,bed_file):
	db = json.load(open(db_file))
	annotated_mutations = mutations
	for i in range(len(mutations["variants"])):
		#var = {'genome_pos': 6140, 'gene_id': 'Rv0005', 'chr': 'Chromosome', 'freq': 0.975609756097561, 'type': 'missense', 'change': '301V>301L'}
		var = mutations["variants"][i]
		print(var)
		if var["gene_id"] in db:
			db_var_match = None
			if var["change"] in db[var["gene_id"]]:
				db_var_match = db[var["gene_id"]][var["change"]]
			elif "frameshift" in var["type"] and "frameshift" in db[var["gene_id"]]:
				db_var_match = db[var["gene_id"]]["frameshift"]
			elif "missense" in var["type"] and "any_missense_codon_%s" % get_missense_codon(var["change"]) in db[var["gene_id"]]:
				db_var_match = db[var["gene_id"]]["any_missense_codon_%s" % get_missense_codon(var["change"])]
			elif "frame" in var["type"] and "any_indel_nucleotide_%s" % get_indel_nucleotide(var["change"]) in db[var["gene_id"]]:
				db_var_match = db[var["gene_id"]]["any_indel_nucleotide_%s" % get_indel_nucleotide(var["change"])]
			elif "stop_gained" in var["type"] and "premature_stop" in db[var["gene_id"]]:
				db_var_match = db[var["gene_id"]]["premature_stop"]
			if db_var_match:
				if "annotation" not in annotated_mutations["variants"][i]:
					annotated_mutations["variants"][i]["annotation"] = {}
				for key in db_var_match:
					annotated_mutations["variants"][i]["annotation"][key] = db_var_match[key]
	for var in deletions:
		bed = load_bed(bed_file,[1,2,3,4,5,6],4)
		for deletion in deletions:
			if deletion["region"] in db and "large_deletion" in db[deletion["region"]]:
				tmp = {"genome_pos":deletion["start"],"gene_id":deletion["region"],"chr":deletion["chr"],"freq":1,"type":"large_deletion","change":"%(chr)s:g.%(start)s_%(end)sdel" % deletion,"annotation":db[deletion["region"]]}
				mutations["variants"].append(tmp)

	#1883443: {u'C': 0.8, u'A': 0.2}
	return annotated_mutations

def annotate_deletions(deletions,mutations,bed_file):

	return mutations
