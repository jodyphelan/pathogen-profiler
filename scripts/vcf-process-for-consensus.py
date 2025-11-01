#! /usr/bin/env python3
import argparse
import pysam

parser = argparse.ArgumentParser(description="Set GT to ./ for SNPs in indel areas.")
parser.add_argument("--vcf", help="Path to the input VCF file.")
parser.add_argument("--output", help="Path to the output BED file.")
parser.add_argument("--af", type=float, default=0.8, help="Allele frequency threshold.")
parser.add_argument("--min-dp", type=int, default=10, help="Minimum depth threshold.")

args = parser.parse_args()



if not args.vcf:
    args.vcf = '/dev/stdin'
if not args.output:
    args.output = '/dev/stdout'

vcf = pysam.VariantFile(args.vcf)  # Open the input VCF file
vcf_lines = []
for var in vcf:
    ad = var.samples[0]['AD']  # Get the Allelic Depth (AD) field
    adf = [d/sum(ad) for d in ad]  # Calculate allele depth fractions
    
    if sum(ad)<args.min_dp:
        var.samples[0]['GT'] = (None, None)
    elif adf[1] > args.af and var.samples[0]['GT'] == (0, 1):
        var.samples[0]['GT'] = (1, 1)  # Modify genotype to homozygous variant
    vcf_lines.append(var)  # Print the modified variant


with pysam.VariantFile(args.output, 'w', header=vcf_lines[0].header) as out_vcf:
    sample = vcf_lines[0].header.samples[0]  # Assuming single sample VCF
    for rec in vcf_lines:
        out_vcf.write(rec)
