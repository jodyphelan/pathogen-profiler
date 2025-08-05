#!/usr/bin/env python
import argparse
import pysam

parser = argparse.ArgumentParser(description="Extract mixed calls from a VCF file and output them in BED format.")
parser.add_argument("--vcf", help="Path to the input VCF file.")
parser.add_argument("--output", help="Path to the output BED file.")
parser.add_argument("--lb", type=float, default=0.1, help="Lower bound for allele frequency.")
parser.add_argument("--ub", type=float, default=0.9, help="Upper bound for allele frequency.")
args = parser.parse_args()

def extract_mixed_calls(vcf_file):
    afs = []
    afs = []
    vcf = pysam.VariantFile(vcf_file)
    sample = list(vcf.header.samples)[0]  # Assuming we are interested in the first sample
    for var in vcf.fetch():
        AD = var.samples[sample]["AD"]
        alt_ad = sum(AD[1:]) if len(AD) > 1 else 0
        alt_af = alt_ad / sum(AD) if sum(AD) > 0 else 0
        afs.append({
            'chrom': var.chrom,
            'pos': var.pos,
            'ref': var.ref,
            'alt': var.alts[0] if var.alts else None,
            'alt_af': alt_af,
            'alt_ad': alt_ad,
            'total_ad': sum(AD)
        })
    return afs

def write_bed_file(afs, output_bed, lb_cut, ub_cut):
    with open(output_bed, "w") as bed_file:
        for r in afs:
            if lb_cut <= r['alt_af'] <= ub_cut:
                bed_file.write(f"{r['chrom']}\t{r['pos']-1}\t{r['pos']}\n")

if not args.vcf:
    args.vcf = '/dev/stdin'
if not args.output:
    args.output = '/dev/stdout'

afs = extract_mixed_calls(args.vcf)
write_bed_file(afs, args.output, args.lb, args.ub)