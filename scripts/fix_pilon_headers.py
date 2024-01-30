#! /usr/bin/env python3
import sys
import re
import argparse

def main(args):
    for l in sys.stdin:
        if l[0]=="#":
            # temp fix for issue with pilon vcf header
            if "ID=DP,Number=1,Type=String" in l:
                l = re.sub(r'ID=DP,Number=1,Type=String', r'ID=DP,Number=1,Type=Integer',l)
            if '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' in l:
                l = re.sub(r'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">', r'##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',l)
            if '##INFO=<ID=SVLEN,Number=.,Type=String,Description="Difference in length between REF and ALT alleles">' in l:
                l = re.sub(r'##INFO=<ID=SVLEN,Number=.,Type=String,Description="Difference in length between REF and ALT alleles">', r'##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Insertion length for SVTYPE=INS.">',l)
            if l.startswith('#CHROM'):
                l = re.sub(r'SAMPLE', args.sample,l)
        sys.stdout.write(l)

parser = argparse.ArgumentParser(description='Fix pilon headers',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.set_defaults(func=main)
parser.add_argument('--sample',help='Sample name',required=True)
args = parser.parse_args()
args.func(args)