#! /usr/bin/env python3
import pysam
from collections import defaultdict

variants = defaultdict(list)
input_vcf = pysam.VariantFile("-")
for var in input_vcf:
	variants[(var.chrom,var.pos)].append(var)

new_header = input_vcf.header.copy()
##INFO=<ID=SAF,Number=A,Type=Integer,Description="Number of alternate observations on the forward strand">
##INFO=<ID=SAR,Number=A,Type=Integer,Description="Number of alternate observations on the reverse strand">
info_tags = [
	{'ID':'SAF','Number':'A','Type':'Integer','Description':'Number of alternate observations on the forward strand'},
	{'ID':'SAR','Number':'A','Type':'Integer','Description':'Number of alternate observations on the reverse strand'},	
]

format_tags = [
	{'ID':'GT','Number':'1','Type':'String','Description':'Genotype'},
	{'ID':'DP','Number':'1','Type':'Integer','Description':'Read Depth'},
	{'ID':'AD','Number':'R','Type':'Integer','Description':'Allelic depths (high-quality bases)'},
]



if len(new_header.samples)==0:
	for tag in info_tags:
		new_header.add_meta('INFO',items=tag.items())
	for tag in format_tags:
		new_header.add_meta('FORMAT',items=tag.items())
	new_header.add_sample('sample')
output_vcf = pysam.VariantFile("-", 'w', header=new_header)
for vars in variants.values():
	if "SAF" in vars[0].info:
		output_vcf.write(vars[0])
	else:
		new_var = output_vcf.new_record()
		new_var.chrom = vars[0].chrom
		new_var.pos = vars[0].pos
		new_var.ref = vars[0].ref
		new_var.filter.add("PASS")
		new_var.info['DP'] = vars[0].info['DP']
		if len(vars) == 1:
			new_var.alts = vars[0].alts 
		else:
			new_var.alts = vars[0].alts + vars[1].alts
			# copy the first variant record's INFO and FORMAT fields

		new_var.info['SAF'] = [v.info['DP4'][2] for v in vars]
		new_var.info['SAR'] = [v.info['DP4'][3] for v in vars]
	
		new_var.samples[0]['GT'] = [1,1]
		new_var.samples[0]['DP'] = vars[0].info['DP']
		new_var.samples[0]['AD'] = [vars[0].info['DP4'][0] + vars[0].info['DP4'][1]] + [v.info['DP4'][2]+v.info['DP4'][3] for v in vars]

		output_vcf.write(new_var)
output_vcf.close()
