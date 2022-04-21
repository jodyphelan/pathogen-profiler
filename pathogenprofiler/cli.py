from .fastq import fastq
from .utils import infolog, run_cmd
from .bam import bam
from .db import get_species_db, get_resistance_db
from .fasta import fasta
from .profiler import bam_profiler, fasta_profiler, vcf_profiler
import json

def get_resistance_db_from_species_prediction(args,species_prediction):
    if args.resistance_db:
        return get_resistance_db(args.software_name,args.resistance_db)

    if len(species_prediction['prediction'])>1:
        infolog(f"Multiple species found.\n")
        return None
    if len(species_prediction['prediction'])==0:
        infolog(f"Species classification failed.\n")
        return None
    if len(species_prediction['prediction'])==1:
        infolog("No resistance database was specified. Attempting to use database based on species prediction...\n")
        db_name = species_prediction['prediction'][0]["species"].replace(" ","_")
        conf = get_resistance_db(args.software_name,db_name)
        if not conf:
            infolog(f"No resistance db found for {db_name}.\n")
        return conf
    


def run_profiler(args):
    if args.read1 or args.bam:
        args.bam_file = get_bam_file(args)
        results = bam_profiler(
            conf=args.conf, bam_file=args.bam_file, prefix=args.files_prefix, platform=args.platform,
            caller=args.caller, threads=args.threads, no_flagstat=args.no_flagstat,
            run_delly = args.run_delly, calling_params=args.calling_params,
            coverage_fraction_threshold=args.coverage_fraction_threshold,
            missing_cov_threshold=args.missing_cov_threshold, samclip=args.no_clip,
            min_depth=args.min_depth,delly_vcf_file=args.delly_vcf,call_wg=args.call_whole_genome,
            variant_annotations=args.add_variant_annotations
        )
        results["input_data_source"] = "fastq" if args.read1 else "bam"
    elif args.fasta:
        results = fasta_profiler(conf=args.conf,prefix=args.files_prefix,filename=args.fasta)
        results["input_data_source"] = "fasta"
    elif args.vcf:
        results = vcf_profiler(conf=args.conf,prefix=args.files_prefix,sample_name=args.prefix,vcf_file=args.vcf,delly_vcf_file=args.delly_vcf)
        results["input_data_source"] = "vcf"
    return results

def speciate(args,bam_region=None):
    conf = get_species_db(args.software_name,args.species_db)
    
    if "read1" in vars(args) and args.read1:
        fastq_class = fastq(args.read1,args.read2)
        kmer_dump = fastq_class.get_kmer_counts(args.files_prefix,threads=args.threads)
    elif "bam" in vars(args) and args.bam:
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
    elif "fasta" in vars(args) and args.fasta:
        kmer_dump = fasta(args.fasta).get_kmer_counts(args.files_prefix,threads=args.threads)
    if "output_kmer_counts" not in vars(args):
        args.output_kmer_counts = None
    else:
        args.output_kmer_counts = f"{args.prefix}.kmers.txt" 
    species = kmer_dump.get_taxonomic_support(conf['kmers'],args.output_kmer_counts)
    return {"prediction":species,"species_db_version":json.load(open(conf['version']))}

def get_bam_file(args):
    ### Create bam file if fastq has been supplied ###
    if args.bam is None:
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


def set_platform_params(args):
    if args.platform=="nanopore":
        args.mapper = "minimap2"
        if args.caller=="gatk":
            args.caller = "freebayes"
        args.no_trim=True
        args.run_delly = True
    else:
        if "no_delly" in vars(args):
            args.run_delly = False if args.no_delly else True

    return args

