import logging
from .utils import run_cmd, cmd_out, TempFilePrefix
import pysam
import argparse
from uuid import uuid4
import numpy as np
import tempfile
from .globals import g
from .fastq import Fastq
from .cli import get_bam_file
from .bam import Bam

def robust_bounds(data, k=5):
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    mad_scaled = 1.4826 * mad
    lower = median - k * mad_scaled
    upper = median + k * mad_scaled
    logging.debug(f"Robust bounds for depth: {median} [{lower}, {upper}]")
    return lower, upper

def generate_low_dp_mask(bam: str,ref: str,outfile: str,min_dp: int = 10) -> None:
    refseq = pysam.FastaFile(ref)
    seqname = refseq.references[0]
    dp = np.zeros(refseq.get_reference_length(seqname))
    for l in cmd_out(f"samtools depth {bam}"):
        row = l.strip().split("\t")
        dp[int(row[1])-1] = int(row[2])

    lower, upper = robust_bounds(dp)
    if min_dp>lower:
        lower = min_dp


    masked_positions = []
    for i, d in enumerate(dp):
        if d < lower or d > upper:
            masked_positions.append(i)

    with open(outfile,"w") as O:
        for i in masked_positions:
            O.write(f"{seqname}\t{i}\t{i+1}\n")



def generate_low_dp_mask_vcf(vcf: str,outfile: str,min_dp: int = 10) -> None:
    missing_positions = []
    vcf_obj = pysam.VariantFile(vcf)
    for rec in vcf_obj:
        # use AD field if available
        if 'AD' in rec.samples[0]:
            dp = sum(rec.samples[0]['AD'])
        else:
            dp = rec.samples[0]['DP']
        if dp<min_dp:
            missing_positions.append((rec.chrom,rec.pos))

    # write missing positions to bed file
    with open(outfile,"w") as O:
        for x in missing_positions:
            O.write(f"{x[0]}\t{x[1]}\t{x[1]+1}\n")

def prepare_sample_consensus(
        sample_name: str,
        ref: str,
        input_vcf: str, 
        output_file: str,
        excluded_regions: str,
        low_dp_regions: str = None,

    ) -> str:
    with TempFilePrefix() as tmp:
        tmp_vcf = f"{tmp}.{sample_name}.vcf.gz"
        run_cmd(f"""
            bcftools view {input_vcf} \
                | annotate_maaf.py \
                | bcftools filter -S . -e 'GT="alt" && MAAF<0.7' \
                | snp-gap.py \
                | rename_vcf_sample.py --sample-name {sample_name} \
                | vcf-process-for-consensus.py \
                | bcftools view -v snps -Oz -o {tmp_vcf}
        """)
        run_cmd(f"bcftools index {tmp_vcf}")
        
        run_cmd(f"vcf-extract-mixed-pos-bed.py --vcf {tmp_vcf} --lb 0.2 --ub 0.8 > {tmp_vcf}.mixed_positions.bed ")
        if low_dp_regions:
            mask_cmd = f"-m {low_dp_regions} -m {tmp_vcf}.mixed_positions.bed -m {excluded_regions}"
        else:
            mask_cmd = f"-m {excluded_regions}"


        
        run_cmd(f"bcftools consensus --sample {sample_name} {mask_cmd} -f {ref} {tmp_vcf} | sed 's/>/>{sample_name} /' > {output_file}")
        return output_file

def cli_prepare_sample_consensus(sample: str,input_vcf: str,args: argparse.Namespace, strain: str = None) -> str:
    if (
        strain 
        and ("strain-specific-references" in args.conf) 
        and (strain in g.strain_specific_references)
    ):
        logging.debug(f"Strain specific reference exists for {strain}")
        strain_specific_ref = g.strain_specific_references[strain]
        ref_file = strain_specific_ref.fasta
        excluded_regions = strain_specific_ref.mask
        args.bam = None
        bam_for_consensus = Bam(
            bam_file=get_bam_file(args,strain_specific_ref),
            prefix=args.prefix,
            platform="illumina",
            threads=args.threads
        )
        # bam.call_variants(conf["ref"], caller=args.caller, filters = conf['variant_filters'], threads=args.threads, calling_params=args.calling_params, samclip = args.samclip, cli_args=vars(args))
        input_vcf = bam_for_consensus.call_variants(
            ref_file=ref_file,
            caller=args.caller,
            filters=args.conf['variant_filters'],
            threads=args.threads
        ).filename
    else:
        ref_file = args.conf['ref']
        excluded_regions = args.conf['bedmask']
        bam_for_consensus = args.bam
        


    mask_bed = f"{args.files_prefix}.{sample}.mask.bed"
    
    if hasattr(args,'supplementary_bam') and args.supplementary_bam:
        bam_for_consensus = args.supplementary_bam
    if bam_for_consensus:
        generate_low_dp_mask(f"{bam_for_consensus}",ref_file,mask_bed)
    elif args.low_dp_mask:
        mask_bed = args.low_dp_mask
    elif args.vcf:
        generate_low_dp_mask_vcf(args.vcf,mask_bed)
    else:
        mask_bed = None

    output_file = f'{args.files_prefix}.consensus.fa'
    prepare_sample_consensus(
        sample_name=sample,
        ref=ref_file,
        input_vcf=input_vcf,
        output_file=output_file,
        excluded_regions=excluded_regions,
        low_dp_regions=mask_bed,
    )
    return output_file

def consensus_fasta_to_vcf(consensus_fasta: str, ref: str, outfile: str) -> None:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_aln = f"{tmpdir}/temp_aln.fa"
        run_cmd(f"cat {ref} {consensus_fasta}> {tmp_aln}")
        run_cmd(f"fa2vcf.py {tmp_aln} {outfile}")

def cli_get_consensus_vcf(sample: str,input_vcf: str,args: argparse.Namespace) -> str:
    consensus_file = cli_prepare_sample_consensus(sample,input_vcf,args)
    outfile = f"{args.files_prefix}.masked.vcf"
    consensus_fasta_to_vcf(consensus_file, args.conf['ref'], outfile)
    return outfile

