import pysam
import argparse
import math

parser = argparse.ArgumentParser(description="Add variant annotations to a VCF file.")


parser.add_argument(
    "bam", type=str, help="Input BAM file."
)



args = parser.parse_args()



import pysam
from scipy.stats import ranksums



rc = {"A": "T", "T": "A", "C": "G", "G": "C"}

def read_pos_rank_sum(bam: str, chrom: str, pos: int, ref: str, alt: str):
    bam = pysam.AlignmentFile(bam, "rb")
    
    ref_positions = []
    alt_positions = []
    
    for read in bam.fetch(chrom, pos, pos + 1):
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue
        
        # Get the read position that aligns to the genomic position
        read_pos = read.get_reference_positions(full_length=True)
        try:
            read_idx = read_pos.index(pos-1)
        except ValueError:
            continue  # This read doesn't actually cover the position (e.g., deletion)

        base = read.query_sequence[read_idx]
        
        # Store read-relative position (0 to read length)
        rel_pos = read_idx
        

        if base.upper() == ref.upper():
            ref_positions.append(rel_pos)
        elif base.upper() == alt.upper():
            alt_positions.append(rel_pos)
        else:
            continue  # Ignore mismatches that don't match ref or alt
    

    if len(ref_positions) == 0 or len(alt_positions) == 0:
        print("Insufficient data: need both ref and alt reads.")
        return None

    # Perform Wilcoxon rank-sum test
    stat, p_value = ranksums(alt_positions, ref_positions)
    
    return {
        "statistic": stat,
        "p_value": p_value,
        "alt_positions": alt_positions,
        "ref_positions": ref_positions
    }

def alt_near_read_ends_pct(bam: str, chrom: str, pos: int, alt: str, margin: int = 10):
    bam = pysam.AlignmentFile(bam, "rb")
    
    alt_read_positions = []
    near_end_count = 0
    
    for read in bam.fetch(chrom, pos, pos + 1):
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue
        
        # Get aligned reference-to-read positions
        read_pos = read.get_reference_positions(full_length=True)
        try:
            read_idx = read_pos.index(pos-1)
        except ValueError:
            continue  # skip reads that do not align to the position (e.g. deletions)

        # Skip if out of read bounds
        if read_idx >= len(read.query_sequence):
            continue

        base = read.query_sequence[read_idx]

        
        if base.upper() == alt.upper():
            alt_read_positions.append(read_idx)
            read_len = read.query_length
            # Check if variant is within 'margin' bases from either end
            if read_idx < margin or (read_len - read_idx - 1) < margin:
                near_end_count += 1

    total_alt_reads = len(alt_read_positions)
    if total_alt_reads == 0:
        return {
            "percent_near_ends": None,
            "total_alt_reads": 0,
            "near_end_reads": 0
        }

    percent = 100.0 * near_end_count / total_alt_reads
    return {
        "percent_near_ends": percent,
        "total_alt_reads": total_alt_reads,
        "near_end_reads": near_end_count
    }



def calculate_epp(bam: str, chrom: str, pos: int, alt: str, margin: int = 10):
    bam = pysam.AlignmentFile(bam, "rb")
    EL = 0
    ER = 0

    for read in bam.fetch(chrom, pos, pos + 1):
        # check read name
        if read.is_unmapped or read.is_duplicate or read.is_secondary:
            continue

        ref_pos = read.get_reference_positions(full_length=True)
        try:
            read_idx = ref_pos.index(pos-1)
        except ValueError:
            continue  # read doesn't align over this position

        if read_idx >= len(read.query_sequence):
            continue

        
        

        base = read.query_sequence[read_idx]

        if base.upper() != alt.upper():
            continue

        read_len = read.query_length
        if read_idx < margin:
            EL += 1
        elif (read_len - read_idx - 1) < margin:
            ER += 1


    n = EL + ER
    if n == 0:
        return {
            "EPP": None,
            "EL": EL,
            "ER": ER
        }

    delta = abs((EL / n) - 0.5)
    p = 2 * math.exp(-2 * n * delta ** 2)
    p = min(max(p, 1e-300), 1.0)  # Clamp for numerical safety

    epp = -10 * math.log10(p)
    return {
        "EPP": epp,
        "EL": EL,
        "ER": ER
    }



vcf_in = pysam.VariantFile("-")
header = vcf_in.header

# add float INFO field
tags = [
    {'ID':'AEP','Number':'1','Type':'Float','Description':'Allele End Proximity'},
    {'ID':'AEPT','Number':'1','Type':'Integer','Description':'Total Alt Reads'},
    {'ID':'AEPE','Number':'1','Type':'Integer','Description':'Alt Reads Near Ends'},
]
for tag in tags:
    header.add_meta('INFO',items=tag.items())


vcf_out = pysam.VariantFile("-", "w", header=header)

for var in vcf_in:
    x = alt_near_read_ends_pct(
        bam=args.bam,
        chrom=var.chrom,
        pos=var.pos,
        alt=var.alts[0],
        margin=10
    )
    var.info["AEP"] = x["percent_near_ends"]
    var.info["AEPT"] = x["total_alt_reads"]
    var.info["AEPE"] = x["near_end_reads"]
    vcf_out.write(var)
    

