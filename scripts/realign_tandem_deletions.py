#! /usr/bin/env python
import pysam
import argparse
import logging
from pathogenprofiler.gff import load_gff
from collections import defaultdict
from tqdm import tqdm
import bisect

logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(description='Right-align deletions in a VCF file')
parser.add_argument('input', help='VCF file with deletions to right-align')
parser.add_argument('ref', help='Reference genome in fasta format')
parser.add_argument('gff', help='Gff annotation file')
parser.add_argument('output', help='Output VCF file with right-aligned deletions')

args = parser.parse_args()

def check_coordinates(var: pysam.VariantRecord,refseq: pysam.FastaFile, direction: str) -> tuple[int,int,str,tuple]:
    original_tuple = (var.start, var.stop, var.ref, var.alts)
    if var.alts[0] != '<DEL>' and len(var.ref)==len(var.alts[0]):
        return original_tuple
    start = var.start + 1
    if var.alts[0] == '<DEL>':
        end = var.stop
    else:
        end = var.start + len(var.ref) + 1
    svlen = end-start
    if svlen<=1:
        return original_tuple
    deleted_seq = refseq.fetch(var.chrom, start, end-1)
    if direction == "right":
        nextseq = refseq.fetch(var.chrom, end-1, end+svlen-2)
    else:
        nextseq = refseq.fetch(var.chrom, start-svlen, start-1)
    logging.info((svlen,deleted_seq,nextseq))
    if deleted_seq == nextseq:
        if direction == "right":
            last_del_nuc = deleted_seq[-1]
            ref = last_del_nuc + deleted_seq
            if var.alts[0] == '<DEL>':
                alt = var.alts[0]
                ref = ref[0]
            else:
                alt = last_del_nuc
            return end-2, end+svlen-1, ref, (alt,)
        else:
            first_del_nuc = deleted_seq[0]
            ref = deleted_seq + first_del_nuc
            if var.alts[0] == '<DEL>':
                alt = var.alts[0]
                ref = ref[-1]
            else:
                alt = first_del_nuc
            return start-svlen+1, start, ref, (alt,)
    else:
        return original_tuple

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
    # logging.info(f"Processing variant {var.chrom}:{var.start}-{var.stop}")
    # find if the variant overlaps with a gene end using the bisect module
    gene_ends_list = gene_ends[var.chrom]
    gene_end_index = bisect.bisect_left(gene_ends_list, var.start) 
    if gene_ends_list[gene_end_index] > var.start and gene_ends_list[gene_end_index] < var.stop:
        logging.info(f"Variant overlaps with gene end {genes[gene_end_index].name}")
        direction = "right" if genes[gene_end_index].strand == "+" else "left"
        start,end,ref,alts = check_coordinates(var,refseq,direction)
        var.start = start
        var.stop = end
        var.ref = ref
        var.alts = alts
    vcf_out.write(var)