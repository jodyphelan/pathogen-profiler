#! /usr/bin/env python3
import argparse
import pysam

parser = argparse.ArgumentParser(description="Set GT to ./ for SNPs in indel areas.")
parser.add_argument("--vcf", help="Path to the input VCF file.")
parser.add_argument("--output", help="Path to the output BED file.")
parser.add_argument("--padding", type=int, default=50, help="Padding around indels.")
parser.add_argument("--remove", action="store_true", help="Remove SNPs in indel areas.")

args = parser.parse_args()

def is_indel(var):
    ref = var.ref
    alts = var.alts
    for alt in alts:
        if len(ref) != len(alt):
            return True
    return False

def load_vcf_lines(vcf_file):
    vcf = pysam.VariantFile(vcf_file, 'r')
    lines = []
    for rec in vcf.fetch():
        lines.append(rec)
    return lines

def extract_indel_area(vcf_lines, padding):
    indel_area = set()
    for rec in vcf_lines:
        if is_indel(rec):
            for pos in range(rec.start - padding, rec.stop + padding):
                indel_area.add((rec.chrom, pos))
    return indel_area

def set_gt_to_missing(vcf_lines, indel_area, output_file):
    with pysam.VariantFile(output_file, 'w', header=vcf_lines[0].header) as out_vcf:
        sample = vcf_lines[0].header.samples[0]  # Assuming single sample VCF
        for rec in vcf_lines:
            if not is_indel(rec):
                if (rec.chrom, rec.start) in indel_area:
                    if args.remove:
                        continue
                    # Set GT to ./.
                    rec.samples[sample]['GT'] = (None, None)
            out_vcf.write(rec)

if not args.vcf:
    args.vcf = '/dev/stdin'
if not args.output:
    args.output = '/dev/stdout'

vcf_lines = load_vcf_lines(args.vcf)
indel_area = extract_indel_area(vcf_lines, args.padding)
set_gt_to_missing(vcf_lines, indel_area, args.output)