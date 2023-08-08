#! /usr/bin/env python
import pysam


vcf_in = pysam.VariantFile("-")
header = vcf_in.header

# add float INFO field
tag = {'ID':'MAAF','Number':'1','Type':'Float','Description':'Major Alternative Allele Frequency'}
header.add_meta('INFO',items=tag.items())

vcf_out = pysam.VariantFile("-", "w", header=header)
for var in vcf_in:
    # calculate MAF from FMT/AD
    maaf = max(var.samples[0]['AD'][1:])/sum(var.samples[0]['AD'])
    var.info['MAAF'] = maaf
    vcf_out.write(var)
