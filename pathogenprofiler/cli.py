from .fastq import fastq
from .utils import run_cmd
from .bam import bam
from .db import get_species_db

def speciate(args,bam_region=None):
    if args.external_species_db:
        conf = get_species_db(args.software_name,args.external_species_db)
    else:
        conf = get_species_db(args.software_name,args.species_db)
    
    if args.bam==None:
        fastq_class = fastq(args.read1,args.read2)
        kmer_dump = fastq_class.get_kmer_counts(args.files_prefix,threads=args.threads)

    else:
        if bam_region:
            region_arg = bam_region if bam_region else ""
            run_cmd(f"samtools view -b {args.bam} {region_arg} | samtools fastq > {args.files_prefix}.tmp.fq")
            kmer_dump = fastq(f"{args.files_prefix}.tmp.fq").get_kmer_counts(args.files_prefix,threads=args.threads)
        else:
            bam_class = bam(args.bam,args.files_prefix,"illumina")
            if bam_class.filetype=="cram":
                run_cmd(f"samtools fastq {bam_class.bam_file} > {args.files_prefix}.tmp.fq")
                kmer_dump = fastq(f"{args.files_prefix}.tmp.fq").get_kmer_counts(args.files_prefix,threads=args.threads)
            else:
                kmer_dump = bam_class.get_kmer_counts(args.files_prefix,threads=args.threads)
    
    species = kmer_dump.get_taxonomic_support(conf['kmers'])
    return species

def get_bam_file(args):
    ### Create bam file if fastq has been supplied ###
    if args.bam==None:
        if args.read1 and args.read2 and args.no_trim:
            # Paired + no trimming
            fastq_obj = fastq(args.read1,args.read2)
        elif args.read1 and args.read2 and not args.no_trim:
            # Paired + trimming
            untrimmed_fastq_obj = fastq(args.read1,args.read2)
            fastq_obj = untrimmed_fastq_obj.trim(args.files_prefix,threads=args.threads)
        elif args.read1 and not args.read2 and args.no_trim:
            # Unpaired + trimming
            fastq_obj = fastq(args.read1,args.read2)
        elif args.read1 and not args.read2 and not args.no_trim:
            # Unpaired + trimming
            untrimmed_fastq_obj = fastq(args.read1)
            fastq_obj = untrimmed_fastq_obj.trim(args.files_prefix,threads=args.threads)
        else:
            exit("\nPlease provide a bam file or a fastq file(s)...Exiting!\n")
        bam_obj = fastq_obj.map_to_ref(
            ref_file=args.conf["ref"], prefix=args.files_prefix,sample_name=args.prefix,
            aligner=args.mapper, platform=args.platform, threads=args.threads
        )
        bam_file = bam_obj.bam_file
    else:
        bam_file = args.bam

    return bam_file