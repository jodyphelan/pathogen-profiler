#! /usr/bin/env python
import pathogenprofiler as pp
import argparse
import sys

parser = argparse.ArgumentParser(description='Get genome chunks')
parser.add_argument('--fasta',help='Genome accession',required=True)
parser.add_argument('--num',help='Number of chunks',type=int,required = True)

args = parser.parse_args()

# Get the genome chunks
for region in pp.get_genome_chunks(args.fasta,args.num):
    sys.stdout.write(region + "\t" + region.replace(":","_").replace("-","_") + "\n") 