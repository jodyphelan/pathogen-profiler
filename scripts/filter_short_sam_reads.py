#! /usr/bin/env python3
import sys
import re
import argparse

def main(args):
    for l in sys.stdin:
        row = l.strip().split()
        if l[0]=="@":
            sys.stdout.write(l)
            continue
        num_matches = sum([int(x) for x in re.findall("(\d+)M",row[5])])
        if num_matches<args.min_match:
            continue

        sys.stdout.write(l)

parser = argparse.ArgumentParser(description='add required annotations',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--min-match',type=int,default=100,help='Reference file (lofreq required)')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
