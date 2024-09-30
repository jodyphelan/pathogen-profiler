import pysam

def check_bam_for_rg(bam_file) -> None:
    """
    Check if the BAM file has a read group.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    if 'RG' not in bam.header:
        raise Exception(f"No read group found in {bam_file}. Please create create your bam with a read group or add using `samtools addreplacerg`.")
    
def check_vcf_chrom_match(vcf_file, ref_file) -> None:
    """
    Check if the chromosomes in the VCF file match the reference file.
    """
    vcf = pysam.VariantFile(vcf_file)
    ref = pysam.FastaFile(ref_file)
    for chrom in vcf.header.contigs:
        if chrom not in ref.references:
            raise Exception(f"Chromosome {chrom} in VCF file {vcf_file} not found in reference file {ref_file}.")

def check_bam_chrom_match(bam_file, ref_file) -> None:
    """
    Check if the chromosomes in the BAM file match the reference file.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    ref = pysam.FastaFile(ref_file)
    for chrom in bam.references:
        if chrom not in ref.references:
            raise Exception(f"Chromosome {chrom} in BAM file {bam_file} not found in reference file {ref_file}.")