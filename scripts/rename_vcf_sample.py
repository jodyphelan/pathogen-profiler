#! /usr/bin/env python3
import argparse
import sys

parser = argparse.ArgumentParser(description='add required annotations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--sample-name',type=str,help='Sample name')
args = parser.parse_args()

for l in sys.stdin:
    row = l.strip().split()
    if row[0]=="#CHROM":
        row[9] = args.sample_name
        l = "\t".join(row)+"\n"
    sys.stdout.write(l)