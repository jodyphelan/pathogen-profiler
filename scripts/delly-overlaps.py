import sys
import subprocess
import argparse
import pathogenprofiler as pp

def main(args):
	pp.filecheck(args.bcf)
	pp.filecheck(args.bed)
	bcf = pp.delly_bcf(args.bcf)
	overlaps = bcf.overlap_bed(args.bed)
	for overlap in overlaps:
		overlap["len"] = int(overlap["end"]) - int(overlap["start"])
		overlap["sample_name"] = args.sample_name
		print("%(sample_name)s\t%(chr)s\t%(start)s\t%(end)s\t%(len)s\t%(region)s" % overlap)
parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('bcf',type=str,help='BCF or VCF file')
parser.add_argument('bed',type=str,help='BED file containing regions of interest')
parser.add_argument('sample_name',type=str,help='BED file containing regions of interest')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
