#! /usr/bin/env python
import pathogenprofiler as pp
import argparse

def main(args):
	if args.bed:
		pp.split_bed(args.bed,args.size,reformat=args.reformat)
	else:
		fasta = pp.fasta(args.fasta)
		fasta.splitchr(args.size,reformat=args.reformat)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('fasta', help='Fasta file')
parser.add_argument('size', type=int,help='Chunk Size')
parser.add_argument('--bed', default=None, help='Bed file')
parser.add_argument('--reformat', action='store_true', help='Bed file')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
