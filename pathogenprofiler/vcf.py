from .utils import *
from collections import defaultdict
import re
re_seq = re.compile("([0-9\-]*)([A-Z\*]+)")
re_I = re.compile("([A-Z\*]+)")
number_re = re.compile("[0-9\-]+")

def parse_mutation(x):
	tmp = x.split(">")
	aa_changed = True if len(tmp)>1 else False
	re_obj = re_seq.search(tmp[0])
	change_num = re_obj.group(1)
	ref_aa = re_obj.group(2)
	alt_aa = re_seq.search(tmp[1]).group(2) if aa_changed else None
	return change_num,ref_aa,alt_aa

class bcf:
	def __init__(self,filename,prefix=None,threads=4):
		self.samples = []
		add_arguments_to_self(self,locals())
		if prefix==None:
			if filename[-4:]==".bcf":
				self.prefix = filename[:-4]
			elif filename[-5:]==".gbcf":
				self.prefix = filename[:-5]
			elif filename[-7:]==".vcf.gz":
				self.prefix = filename[-7:]==".vcf.gz"
			elif filename[-4:]==".vcf":
				self.prefix = filename[-4:]==".vcf"
			else:
				self.prefix = filename
		else:
			self.prefix = prefix
		self.prefix = self.prefix
		self.temp_file = get_random_file()
		index_bcf(filename,self.threads)
		run_cmd("bcftools query -l %(filename)s > %(temp_file)s" % vars(self))
		for l in open(self.temp_file):
			self.samples.append(l.rstrip())
		rm_files([self.temp_file])
	def del_pos2bed(self):
		self.del_bed = "%s.del_pos.bed" % self.prefix
		OUT = open(self.del_bed,"w")
		j = 0
		for l in cmd_out("bcftools view --threads %(threads)s -Ou -v indels %(filename)s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' | awk 'length($3)>1'" % vars(self)):
			j+=1
			row = l.split()
			start_pos = int(row[1])+1
			for i in range(start_pos,start_pos+len(row[2])-1):
				OUT.write("%s\t%s\t%s\n" % (row[0],i-1,i))
		if j==0:
			OUT.write("dummy\t1\t1\n")
		OUT.close()
		return self.del_bed
	def load_csq(self,ann_file=None):
		ann = defaultdict(dict)
		if ann_file:
			for l in open(ann_file):
				#chrom pos gene gene/codon_pos
				row = l.rstrip().split()
				ann[row[0]][int(row[1])] = (row[2],row[3])

		nuc_variants = self.load_variants()
		variants = {s:[] for s in self.samples}
		for line in cmd_out("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE\\t%%TBCSQ\\t%%TGT\\t%%AD]\\n' %s" % self.filename):
			row = line.split()
			chrom = row[0]
			pos = int(row[1])
			ref = row[2]
			alts = row[3].split(",")
			alleles = [ref]+alts
			if chrom in ann and pos in ann[chrom]:
				ann_pos = int(ann[chrom][pos][1])
				ann_gene = ann[chrom][pos][0]
			else:
				ann_pos = None
				ann_gene = None
			if len(row)==4:
				for alt in alts:
					if chrom in ann and pos in ann[chrom]:
						cng = "%s%s>%s" % (ann_pos,ref,alt)
						for sample in self.samples:
							if sample in nuc_variants[chrom][pos] and alt in nuc_variants[chrom][pos][sample]:
								variants[sample].append({"sample":sample,"gene_id":ann_gene,"chr":chrom,"genome_pos":pos,"type":"non_coding","change":cng,"freq":nuc_variants[chrom][pos][sample][alt]})
					else:
						log("ERROR in loading alts",True)
				continue

			for i in range(4,len(row)-4,5):
				sample = row[i]
				info = row[i+1].split("|") if row[i+1]!="." else row[i+2].split("|")
				call1,call2 = row[i+3].split("/")
				ad = [int(x) if x!="." else 0 for x in row[i+4].split(",")]

				adr = {alleles[i]:d/sum(ad) for i,d in enumerate(ad)}
				if row[i+1][0]=="@": continue
				if info[-1]=="pseudogene": continue
				gene = info[1]
				if info[0]=="intron":continue
				if info[0]=="coding_sequence":
					cng = "%s%s>%s" % (ann_pos,call1,call2)
					variants[sample].append({"sample":sample,"gene_id":ann_gene,"chr":chrom,"genome_pos":pos,"type":"non_coding","change":cng,"freq":adr[call2]})
				elif  "missense" in info[0] or "start_lost" in info[0] or "stop_gained" in info[0]:
					variants[sample].append({"sample":sample,"gene_id":gene,"chr":chrom,"genome_pos":pos,"type":info[0],"change":info[5],"freq":adr[call2]})
				elif "synonymous" in info[0] or info[0]=="stop_retained":
					change_num,ref_nuc,alt_nuc =  parse_mutation(info[6])
					change = "%s%s>%s" % (ann_pos,ref_nuc,alt_nuc) if ann_pos else "%s%s>%s" % (pos,ref_nuc,alt_nuc)
					variants[sample].append({"sample":sample,"gene_id":gene,"chr":chrom,"genome_pos":pos,"type":info[0],"change":change,"freq":adr[call2]})
				elif "frame" in info[0] or "stop_lost" in info[0]:
					if len(info)<6:
						if chrom in ann and pos in ann[chrom]:
							change = "%s%s>%s" % (pos,ref,call2)
							variants[sample].append({"sample":sample,"gene_id":gene,"chr":chrom,"genome_pos":pos,"type":info[0],"change":change,"freq":adr[call2]})
					else:
						variants[sample].append({"sample":sample,"gene_id":gene,"chr":chrom,"genome_pos":pos,"type":info[0],"change":info[6],"freq":adr[call2]})
				elif info[0]=="non_coding":
					if chrom in ann and pos in ann[chrom]:
						gene = ann[chrom][pos][0]
						gene_pos = ann[chrom][pos][1]
						change = "%s%s>%s" % (gene_pos,ref,call2)
						variants[sample].append({"sample":sample,"gene_id":gene,"chr":chrom,"genome_pos":pos,"type":info[0],"change":change,"freq":adr[call2]})
				else:
					log(line)
					log(info[0]+"\n")
					log("Unknown variant type...Exiting!\n",True)

		return variants
	def load_variants(self,chrom=None,pos=None):
		variants = defaultdict(lambda:defaultdict(lambda:defaultdict(dict)))
		raw_variants = defaultdict(lambda:defaultdict(lambda:defaultdict(dict)))
		if chrom and pos:
			cmd = "bcftools view --threads %(threads)s %s %s:%s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%TGT:%%AD]\\n'  | sed 's/\.\/\./N\/N/g'" % (self.filename,chrom,pos)
		else:
			cmd = "bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%TGT:%%AD]\\n' %s  | sed 's/\.\/\./N\/N/g'" % self.filename
		for l in cmd_out(cmd):
			row = l.split()
			alts = row[3].split(",")
			alleles = [row[2]]+alts
			for i in range(len(self.samples)):
				calls,ad = row[i+4].split(":")
				call1,call2 = calls.split("/")
				if calls=="N/N":
					raw_variants[row[0]][row[1]][self.samples[i]]["N"] = 1.0
					continue
				elif calls=="%s/%s" % (row[2],row[2]) and ad==".":
					raw_variants[row[0]][row[1]][self.samples[i]][row[2]] = 1.0
					continue
				ad = [int(x) if x!="." else 0 for x in ad.split(",")]
				sum_ad = sum(ad)
				for j in range(1,len(alleles)):
					if ad[j]==0: continue
					raw_variants[row[0]][row[1]][self.samples[i]][alleles[j]] = ad[j]/sum_ad
		for tchrom in raw_variants:
			for tpos in raw_variants[tchrom]:
				variants[tchrom][int(tpos)] = raw_variants[tchrom][tpos]
		if chrom and pos and len(variants)==0:
			log("Variant not found",True)
		if chrom and pos:
			return variants[chrom][int(pos)]
		else:
			return variants
	def get_positions(self):
		results = []
		for l in cmd_out("bcftools query -f '%%CHROM\\t%%POS\\n' %s" % self.filename):
			row = l.split()
			results.append((row[0],int(row[1])))
		return results


class delly_bcf(bcf):
	def __init__(self,filename):
		 bcf.__init__(self,filename)
	def get_robust_calls(self):
		results = []
		for l in cmd_out("bcftools query -f '%%CHROM\\t%%POS\\t[%%END\\t%%GT\\t%%DR\\t%%DV\\t%%RR\\t%%RV]\\n' %(filename)s" % vars(self)):
			row = l.split()
			if row[3]!="1/1":continue
			if int(row[2])-int(row[1])>100000: continue
			results.append(row)
		return results
	def overlap_bed(self,bed_file):
		results = []
		bed = load_bed(bed_file,[1,2,3,4,5],4)
		calls = self.get_robust_calls()
		for call in calls:
			set_call_pos = set(range(int(call[1]),int(call[2])))
			for region in bed:
				if bed[region][0]!=call[0]: continue
				set_region_pos = set(range(int(bed[region][1]),int(bed[region][2])))
				intersect = set_call_pos.intersection(set_region_pos)
				if len(intersect)>1:
					results.append({"chr":call[0],"region":region,"start":min(intersect),"end":max(intersect)})
		return results
