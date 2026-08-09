#! /usr/bin/env python3
import pysam
import re
import sys

infile = sys.argv[1] if len(sys.argv)>1 else "-"
vcf_in = pysam.VariantFile(infile)
vcf_out = pysam.VariantFile("-", "w", header=vcf_in.header)

regex = re.compile(r"^(p\.Val1|p\.Leu1|p\.Ile1)([A-Z?*]+)")

with vcf_out as out:
    for rec in vcf_in:
        for i,ann in enumerate(rec.info['ANN']):
            ann_fields = ann.split("|")
            if r:=regex.match(ann_fields[10]):
                ann_fields[10] = "p.Met1" + r.group(2)
                if ann_fields[1]=='initiator_codon_variant':
                    ann_fields[10] = "p.Met1="
                ann = "|".join(ann_fields)
                new_info = list(rec.info['ANN'])
                new_info[i] = ann
                rec.info['ANN'] = new_info
        out.write(rec)