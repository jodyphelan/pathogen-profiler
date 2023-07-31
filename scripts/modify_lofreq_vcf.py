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
tags = [
	{'ID':'SAF','Number':'A','Type':'Integer','Description':'Number of alternate observations on the forward strand'},
	{'ID':'SAR','Number':'A','Type':'Integer','Description':'Number of alternate observations on the reverse strand'},	
]
for tag in tags:
	new_header.add_meta('INFO',items=tag.items())
new_header.add_sample('sample')
output_vcf = pysam.VariantFile("-", 'w', header=new_header)
for vars in variants.values():
	# create new variant record
	new_var = output_vcf.new_record()
	# copy the first variant record
	new_var.chrom = vars[0].chrom
	new_var.pos = vars[0].pos
	new_var.ref = vars[0].ref
	new_var.alts = vars[0].alts + vars[1].alts
	new_var.filter.add("PASS")
	# copy the first variant record's INFO and FORMAT fields
	new_var.info['DP'] = vars[0].info['DP']
	new_var.info['SAF'] = [v.info['DP4'][2] for v in vars]
	new_var.info['SAR'] = [v.info['DP4'][3] for v in vars]
	output_vcf.write(new_var)
output_vcf.close()
