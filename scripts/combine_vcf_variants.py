#! /usr/bin/env python3
import pysam
from typing import List
from collections import Counter, defaultdict
import argparse
from pathogenprofiler.models import GenomePosition
from pathogenprofiler.gff import load_gff, Exon
from itertools import product
import logging

parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcf',type=str,help='')
parser.add_argument('--gff',type=str,help='',required = True)
parser.add_argument('--ref',type=str,help='',required = True)
parser.add_argument('--out',type=str,help='')
parser.add_argument('--bam',type=str,help='')
parser.add_argument('--debug',action='store_true',help='Print debug information')
args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG)
else:
    logging.basicConfig(level=logging.INFO)


def get_alleles(
    read:pysam.libcalignedsegment.AlignedSegment,
    positions=List[GenomePosition],
    bystrand=False
):

    alleles = {}
    for read_pos, ref_pos, read_nt in read.get_aligned_pairs(with_seq=True):
        if read_nt is None:
            continue
        if read_pos is None:
            continue
        p = GenomePosition(chrom=read.reference_name,pos=ref_pos+1)
        if p not in positions:
            continue
        if read.query_qualities[read_pos] < 12:
            continue
        
        if read_nt.islower():
            read_nt = read.query_sequence[read_pos].upper()
        
        if bystrand:
            alleles[p] = read_nt.lower() if read.is_reverse else read_nt.upper()
        else:
            alleles[p] = read_nt

    return [alleles.get(p) for p in positions]


def get_haplotype_counts(
    bam:pysam.AlignmentFile,
    positions:List[GenomePosition],
    bystrand=False
):
    allele_combinations = []
    chrom = positions[0].chrom
    start = min(positions).pos
    end = max(positions).pos
    for read in bam.fetch(contig=chrom,start=start,end=end):
        logging.debug((read.qname, read.mapping_quality<10, read.is_secondary, read.is_supplementary, read.is_unmapped))
        if read.mapping_quality < 10:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.is_unmapped:
            continue
        alleles = get_alleles(read,positions=positions,bystrand=bystrand)
        logging.debug(alleles)
        if all(alleles):
            allele_combinations.append(''.join(alleles))
    counts = Counter(allele_combinations)
    logging.debug(counts)
    return counts



def get_overlapping_exons(chrom: str,pos: int,exons: List[Exon]):
    exons = [e for e in exons if e.start<=pos and e.end>=pos and chrom==e.chrom]
    if len(exons)==0:
        return None
    else:
        return exons[0]

def get_codon_pos(chrom: str,pos: int,exons: List[Exon]):
    e = get_overlapping_exons(chrom,pos,exons)
    if e==None:
        return (None,None)
    logging.debug(f"{chrom}:{pos} {vars(e)}")
    if e.strand=="+":
        codon_pos = (pos-e.start-e.phase)//3 + 1
    else:
        codon_pos = (e.end + e.phase - pos )//3 + 1
    return (e.id,codon_pos)
        

gff = load_gff(args.gff)
exons = []
for gene in gff:
    for transcript in gene.transcripts:
        for i,exon in enumerate(transcript.exons):
            exon.id = f"{gene.gene_id}_exon{i+1}"
            exons.append(exon)

def has_md_tag(bam):
    for read in bam.fetch():
        break
    if read.has_tag('MD'):
        return True
    else:
        return False
    
ref = pysam.FastaFile(args.ref)
coding_variants = defaultdict(list)
other_variants = []
vcf = pysam.VariantFile(args.vcf) if args.vcf else pysam.VariantFile('-')

if "AF" not in vcf.header.info.keys():
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
    vcf.header.add_line('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">')

for var in vcf:
    logging.debug(var)
    gene,cpos = get_codon_pos(var.chrom,var.pos,exons)
    if var.rlen!=1 or len(var.alts[0])!=1:
        other_variants.append(var)
    elif gene==None:
        other_variants.append(var)
    else:
        coding_variants[(gene,cpos)].append(var)




for key,variants in coding_variants.items():
    logging.debug((key, variants))
    if len(variants)==1:
        other_variants.append(variants[0])
        continue
    all_positions = [GenomePosition(chrom=v.chrom,pos=v.pos) for v in variants]
    positions = sorted(list(set(all_positions)))
    if len(positions)==2 and positions[0].pos == positions[1].pos-2:
        positions.insert(1,GenomePosition(chrom=positions[0].chrom,pos=positions[0].pos+1))

    ref_hap = ''.join([ref.fetch(p.chrom,p.pos-1,p.pos) for p in positions])
    if args.bam and has_md_tag(pysam.AlignmentFile(args.bam)):
        haplotypes_by_strand = get_haplotype_counts(pysam.AlignmentFile(args.bam),positions,bystrand=True)
        logging.debug(haplotypes_by_strand)
    else:
        # alt is just a combination of all the alt alleles
        alt_hap = ''.join([v.alts[0] for v in variants])
        ds = defaultdict(list)
        for v in variants:
            ds[v.pos].append(v.alts[0])
        v = variants[0]
        if 'DP4' in v.info:
            haplotypes_by_strand = {
                ref_hap: int(v.info['DP4'][0])+int(v.info['DP4'][1])
            }
        elif 'AD' in v.samples[0]:
            haplotypes_by_strand = {
                ref_hap: int(v.samples[0]['AD'][0])
            }
        if len(ds)==2 and list(ds)[0]+2 == list(ds)[1]:
            midpos = list(ds)[0]+1
            # insert the ref allele in between
            ds[midpos] = ref.fetch(v.chrom,midpos-1,midpos)
            # sort
            ds = dict(sorted(ds.items()))        

        for alts in product(*ds.values()):
            hap = ''.join(alts)
            if 'DP4' in v.info:
                haplotypes_by_strand[hap] = int(v.info['DP4'][2])+int(v.info['DP4'][3])
            elif 'SAF' in v.info:
                haplotypes_by_strand[hap] = int(v.info['SAF'][0])
                haplotypes_by_strand[hap.lower()] = int(v.info['SAR'][0])
            elif 'AD' in v.samples[0]:
                haplotypes_by_strand[hap] = int(v.samples[0]['AD'][1])


    haplotypes = Counter()
    for h,count in haplotypes_by_strand.items():
        haplotypes[h.upper()] += count
    dp = sum(haplotypes.values())

    if dp==0:
        # can happen if malformed VCF (e.g. https://github.com/jodyphelan/pathogen-profiler/issues/63)
        continue

    ref_fwd = haplotypes_by_strand.get(ref_hap.upper(),0)
    ref_rev = haplotypes_by_strand.get(ref_hap.lower(),0)
    for i,(hap,count) in enumerate(haplotypes.items()):
        hap_fwd = haplotypes_by_strand.get(hap.upper(),0)
        hap_rev = haplotypes_by_strand.get(hap.lower(),0)
        if hap.upper()==ref_hap.upper(): 
            continue
        variant = variants[0].copy()
        variant.alts = (hap,)
        variant.ref = ref_hap
        variant.info.update({'AF':count/dp})
        if 'DP4' in variant.info:
            variant.info['DP4'] = [ref_fwd,ref_rev,hap_fwd,hap_rev]
        if 'SAF' in variant.info:
            variant.info['SAF'] = hap_fwd
            variant.info['SAR'] = hap_rev
        if 'ADF' in var.samples[0]:
            variant.samples[0]['ADF'] = (ref_fwd,hap_fwd)
            variant.samples[0]['ADR'] = (ref_rev,hap_rev)
        if 'AD' in variant.samples[0]:
            variant.samples[0]['AD'] = [dp - count, hap_fwd + hap_rev]
        
        other_variants.append(variant)


new_vcf = pysam.VariantFile(args.out,'w',header = vcf.header) if args.out else pysam.VariantFile('-','w',header = vcf.header)
for var in other_variants:
    new_vcf.write(var)
new_vcf.close()