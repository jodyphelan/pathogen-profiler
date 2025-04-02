#! /usr/bin/env python
import pysam


vcf_in = pysam.VariantFile("-")
header = vcf_in.header

# add float INFO field
tag = {'ID':'MAAF','Number':'1','Type':'Float','Description':'Major Alternative Allele Frequency'}
header.add_meta('INFO',items=tag.items())

vcf_out = pysam.VariantFile("-", "w", header=header)
for var in vcf_in:
    # short variant
    if 'AD' in var.samples[0]:
        if var.samples[0]['AD'] == (None,):
            maaf = 0
        else:
        # calculate MAF from FMT/AD
            maaf = max(var.samples[0]['AD'][1:])/sum(var.samples[0]['AD'])

    elif 'SVTYPE' in var.info:
        # long variant
        # use DR:DV:RR:RV
        if var.samples[0]['RR'] == (None,):
            maaf = 0
        else:
            maaf = (var.samples[0]['DV'] + var.samples[0]['RV'])/(var.samples[0]['RR'] + var.samples[0]['RV'] + var.samples[0]['DR'] + var.samples[0]['DV'])
    else:
        # panic
        raise ValueError("No AD or SVTYPE in INFO field")
    var.info['MAAF'] = maaf
    vcf_out.write(var)
