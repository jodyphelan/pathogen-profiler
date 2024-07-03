import pysam

def check_bam_for_rg(bam_file) -> None:
    """
    Check if the BAM file has a read group.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    if 'RG' not in bam.header:
        raise Exception(f"No read group found in {bam_file}. Please create create your bam with a read group or add using `samtools addreplacerg`.")
    
