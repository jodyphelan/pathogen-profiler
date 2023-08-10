#! /usr/bin/env python3
import sys
from tqdm import tqdm
import argparse
import re

def main(args):
    ad_cutoff = args.fraction
    for l in tqdm(sys.stdin):
        line = l.strip().replace("|","/")
        something_changed = False
        row = line.split()
        if l[0]=="#":
            sys.stdout.write(l)
            continue
        alleles = [row[3]]+row[4].split(",")
        if len(alleles)>9: continue

        if "*" in row[4]:
            something_changed = True
            for i in range(9,len(row)):
                if row[i][0]=="." or row[i][2]==".": continue
                if alleles[int(row[i][0])]=="*" or alleles[int(row[i][2])]=="*":
                    tmp = list(row[i])
                    tmp[0] = "."
                    tmp[2] = "."
                    row[i]="".join(tmp)


        idx = row[8].split(":").index("AD")

        for i in range(9,len(row)):

            fmt = row[i].split(":")
            ad = [int(x) for x in fmt[idx].split(",")]
            total_ad = sum(ad)
            if total_ad==0:continue
            adf = [ad[j]/total_ad for j in range(len(ad))]

            if max(adf)>=ad_cutoff:
                new_gt = adf.index(max(adf))
                if alleles[new_gt]=="*": new_gt = "."
                fmt[0] = f"{new_gt}/{new_gt}"
                something_changed = True
                row[i] = ":".join(fmt)
            else:
                something_changed = True
                fmt[0] = "./."
                row[i] = ":".join(fmt)
        if something_changed:
            sys.stdout.write("\t".join(row)+"\n")
        else:
            sys.stdout.write(l)

parser = argparse.ArgumentParser(description='TBProfiler pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fraction',default=0.7,type=float,help='Fraction of coverage to assign major')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
