#! /usr/bin/env python
import pysam
import argparse
import logging
from pathogenprofiler.gff import load_gff
from collections import defaultdict
from tqdm import tqdm
import bisect

parser = argparse.ArgumentParser(description='Right-align deletions in a VCF file')
parser.add_argument('input', help='VCF file with deletions to right-align')
parser.add_argument('ref', help='Reference genome in fasta format')
parser.add_argument('gff', help='Gff annotation file')
parser.add_argument('output', help='Output VCF file with right-aligned deletions')
parser.add_argument('--debug',action='store_true',help='Print debug information')

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)

def right_align_deletion(var: pysam.VariantRecord, refseq: pysam.FastaFile) -> tuple[int,int,str,tuple]:
    """
    Right-align a deletion variant in a VCF record.
    """
    logging.debug(f"Right-aligned deletion: {var.chrom}:{var.start}-{var.stop} {var.ref} {var.alts}")

    next_base = refseq.fetch(var.chrom,var.stop, var.stop+1)  # Fetch the base at the stop position
    while next_base==var.ref[1]:
        var.ref = var.ref[1:] + next_base
        var.alts = (next_base,)
        var.start += 1
        next_base = refseq.fetch(var.chrom,var.stop, var.stop+1)  # Fetch the base at the stop position
    logging.debug(f"Right-aligned deletion: {var.chrom}:{var.start}-{var.stop} {var.ref} {var.alts}")

vcf_in = pysam.VariantFile(args.input)
header = vcf_in.header
refseq = pysam.FastaFile(args.ref)
genes = load_gff(args.gff)

gene_ends = defaultdict(list)
gene_names = defaultdict(list)

for chrom in refseq.references:
    for gene in genes:
            gene_ends[gene.chrom].append(gene.end)
            gene_names[gene.chrom].append(gene.name)

vcf_out = pysam.VariantFile(args.output, "w", header=header)
for var in vcf_in:
    logging.debug(f"Processing variant {var.chrom}:{var.start}-{var.stop}")
    # find if the variant overlaps with a gene end using the bisect module
    gene_ends_list = gene_ends[var.chrom]
    gene_end_index = bisect.bisect_left(gene_ends_list, var.start) 
    if gene_end_index<len(gene_ends_list) and gene_ends_list[gene_end_index] > var.start and gene_ends_list[gene_end_index] < var.stop:
        logging.debug(f"Variant overlaps with gene end {genes[gene_end_index].name}")
        direction = "right" if genes[gene_end_index].strand == "+" else "left"
        if direction=="right":
            right_align_deletion(var, refseq)

    vcf_out.write(var)