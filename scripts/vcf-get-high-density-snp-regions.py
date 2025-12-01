#! /usr/bin/env python3
import argparse
import pysam

parser = argparse.ArgumentParser(description="Set GT to ./ for SNPs in indel areas.")
parser.add_argument("--vcf", help="Path to the input VCF file.")
parser.add_argument("--ref", help="Path to the reference genome file.",required=True)
parser.add_argument("--output", help="Path to the output BED file.")
parser.add_argument("--window-size", type=int, default=50, help="Window size for high-density SNP regions.")
parser.add_argument("--max-snps", type=int, default=5, help="Maximum number of SNPs in a window.")

args = parser.parse_args()



if not args.vcf:
    args.vcf = '/dev/stdin'
if not args.output:
    args.output = '/dev/stdout'

ref = pysam.FastaFile(args.ref)  # Open the reference genome file
ref_seqname = ref.references[0]  # Assuming single contig/reference
chrom_len = ref.get_reference_length(ref_seqname)

vcf = pysam.VariantFile(args.vcf)  # Open the input VCF file

snp_positions = []
snp_vector = [0] * chrom_len
for record in vcf:
    if record.alts[0] == "*":
        continue  # Skip deletions
    if record.alts and len(record.alts[0]) == 1 and len(record.ref) == 1:  # Check for SNP
        snp_vector[record.pos - 1] = 1  # Mark SNP position
        snp_positions.append(record.pos - 1)

windows = []
for i,pos in enumerate(snp_positions):
    start = max(0, pos - args.window_size // 2)
    end = min(chrom_len, pos + args.window_size // 2)
    num_snps_in_window = sum(snp_vector[start:end])
    if num_snps_in_window > args.max_snps:
        windows.append({
            'chrom': ref_seqname,
            'start': start,
            'end': end,
            'num_snps': num_snps_in_window
        })

def merge_intervals(intervals):
    if not intervals:
        return []
    # Sort intervals by start position
    intervals.sort(key=lambda x: x['start'])
    merged = [intervals[0]]
    for current in intervals[1:]:
        last = merged[-1]
        if current['start'] <= last['end']:  # Overlap
            merged[-1] = {
                'chrom': last['chrom'],
                'start': last['start'],
                'end': max(last['end'], current['end']),
                'num_snps': max(last['num_snps'], current['num_snps'])
            }  # Merge
        else:
            merged.append(current)
    return merged

merged_windows = merge_intervals(windows)

with open(args.output, "w") as bed_file:
    for window in merged_windows:
        bed_file.write("{chrom}\t{start}\t{end}\t{num_snps}\n".format(**window))
