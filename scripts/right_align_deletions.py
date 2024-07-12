#! /usr/bin/env python
import pysam
import argparse

parser = argparse.ArgumentParser(description='Right-align deletions in a VCF file')
parser.add_argument('input', help='VCF file with deletions to right-align')
parser.add_argument('ref', help='Reference genome in fasta format')
parser.add_argument('output', help='Output VCF file with right-aligned deletions')

args = parser.parse_args()

def check_coordinates(var: pysam.VariantRecord,refseq: pysam.FastaFile) -> tuple[int,int,str,tuple]:
    if var.alts[0] != '<DEL>' and len(var.ref)==len(var.alts[0]):
        return var.start, var.stop, var.ref, var.alts[0]
    start = var.start + 1
    if var.alts[0] == '<DEL>':
        end = var.stop
    else:
        end = var.start + len(var.ref) + 1
    svlen = end-start
    deleted_seq = refseq.fetch(var.chrom, start, end-1)
    nextseq = refseq.fetch(var.chrom, end-1, end+svlen-2)
    if deleted_seq == nextseq:
        last_del_nuc = deleted_seq[-1]
        ref = last_del_nuc + deleted_seq
        if var.alts[0] == '<DEL>':
            alt = var.alts[0]
            ref = ref[0]
        else:
            alt = last_del_nuc
        return end-2, end+svlen-1, ref, (alt,)
    else:
        return var.start, var.stop, var.ref, var.alts

vcf_in = pysam.VariantFile(args.input)
header = vcf_in.header
refseq = pysam.FastaFile(args.ref)

vcf_out = pysam.VariantFile(args.output, "w", header=header)
for var in vcf_in:
    start,end,ref,alts = check_coordinates(var,refseq)
    var.start = start
    var.stop = end
    var.ref = ref
    var.alts = alts
    vcf_out.write(var)