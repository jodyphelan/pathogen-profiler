#! /usr/bin/env python
import argparse
import pysam

parser = argparse.ArgumentParser(description='Convert fasta to vcf')
parser.add_argument('aln', help='Alignment file')
parser.add_argument('out', help='Output vcf file')

args = parser.parse_args()



aln = pysam.FastaFile(args.aln)
chrom = aln.references[0]
sample_name = aln.references[1]
chrom_len = aln.get_reference_length(chrom)

variants = []
for i,(ref,alt) in enumerate(zip(aln.fetch(chrom),aln.fetch(aln.references[1]))):
    if ref!=alt:
        variants.append(
            {
                "chrom": aln.references[0],
                "pos": i,
                "ref": ref,
                "alt": alt if alt!="N" else "*"
            }
        )

# create vcf file

header = pysam.VariantHeader()
header.add_sample(sample_name)
header.contigs.add(chrom, length=chrom_len)
header.add_line('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">')
header.add_line('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">')
header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

vcf = pysam.VariantFile(args.out, "w", header=header)

# add a variant
for variant in variants:
    record = vcf.new_record(
        contig=variant["chrom"],
        start=variant["pos"],
        alleles=[variant["ref"], variant["alt"]],
        id=variant['ref'] + str(variant['pos']+1) + (variant['alt'] if variant["alt"]!="*" else variant['ref'] )
    )
    if variant["alt"]=="*":
        record.info["AC"] = [0]
        record.info["AN"] = 0
        record.samples["sample1"]["GT"] = (None)

    else:
        record.info["AC"] = [1]
        record.info["AN"] = 1
        record.samples["sample1"]["GT"] = (1)

    # write GT field
    vcf.write(record)

vcf.close()




