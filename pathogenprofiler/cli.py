from .fastq import Fastq
from .utils import run_cmd, cmd_out
from .bam import Bam
from .db import get_db
from .fasta import Fasta
from .profiler import bam_profiler, fasta_profiler, vcf_profiler
import json
import logging

def get_variant_filters(args):
    filters = {}
    for f in [
        'depth','af','strand','sv_depth','sv_af','sv_len'
    ]:
        if not hasattr(args,f):
            continue
        vals = getattr(args,f).split(",")
        if vals[0]=="-":
            filters[f+"_hard"] = None
            filters[f+"_soft"] = None    
        else:
            filters[f+"_hard"] = float(vals[0]) if "." in vals[0] else int(vals[0])
            filters[f+"_soft"] = float(vals[1]) if "." in vals[1] else int(vals[1])
    return filters



def get_resistance_db_from_species_prediction(args ,species_prediction):
    if args.resistance_db:
        return get_db(args.software_name,args.resistance_db)
    if species_prediction==None:
        logging.info(f"Species classification failed.\n")
        return None
    species_prediction = species_prediction.replace(" ","_")
    conf = get_db(args.software_name,species_prediction)
    if conf is None:
        logging.info(f"No resistance db found for {species_prediction}.\n")
    return conf
    # if len(species_prediction['prediction'])>1:
    #     logging.info(f"Multiple species found.\n")
    #     return None
    # if len(species_prediction['prediction'])==1:
    #     logging.info("No resistance database was specified. Attempting to use database based on species prediction...\n")
    #     db_name = species_prediction['prediction'][0]["species"].replace(" ","_")
    #     conf = get_db(args.software_name,db_name)
        
    
def test_vcf_for_lofreq(vcf_file):
    lofreq = False
    for l in cmd_out(f"bcftools view -h {vcf_file}"):
        if "source=lofreq call" in l:
            lofreq = True
    return lofreq

def run_profiler(args):
    if args.read1 or args.bam:
        args.bam_file = get_bam_file(args)
        results = bam_profiler(
            conf=args.conf, bam_file=args.bam_file, prefix=args.files_prefix, platform=args.platform,
            caller=args.caller, threads=args.threads, no_flagstat=args.no_flagstat,
            run_delly = args.run_delly, calling_params=args.calling_params,
            samclip=args.no_clip,delly_vcf_file=args.delly_vcf,
            call_wg=args.call_whole_genome,variant_annotations=args.add_variant_annotations
        )
        results["input_data_source"] = "fastq" if args.read1 else "bam"
    elif args.fasta:
        results = fasta_profiler(conf=args.conf,prefix=args.files_prefix,filename=args.fasta)
        results["input_data_source"] = "fasta"
    elif args.vcf:
        if test_vcf_for_lofreq(args.vcf):
            tmp_vcf_file = f"{args.files_prefix}.tmp.vcf.gz"
            run_cmd(f"bcftools view {args.vcf} | modify_lofreq_vcf.py | bcftools view -Oz -o {tmp_vcf_file}")
            args.vcf = tmp_vcf_file
        results = vcf_profiler(conf=args.conf,prefix=args.files_prefix,sample_name=args.prefix,vcf_file=args.vcf,delly_vcf_file=args.delly_vcf)
        results["input_data_source"] = "vcf"
    return results

def kmer_speciate(args,bam_region=None):
    conf = get_db(args.software_name,args.species_db)
    if conf==None:
        logging.error(
            f"\nDatabase '{args.species_db}' not found. You may need to load this database first... Exiting!\n"
        )
    
    if "read1" in vars(args) and args.read1:
        fastq = Fastq(args.read1,args.read2)
        kmer_dump = fastq.get_kmer_counts(args.files_prefix,threads=args.threads,max_mem=args.ram,counter = args.kmer_counter)
    elif "bam" in vars(args) and args.bam:
        if bam_region:
            region_arg = bam_region if bam_region else ""
            run_cmd(f"samtools view -b {args.bam} {region_arg} | samtools fastq > {args.files_prefix}.tmp.fq")
            kmer_dump = Fastq(f"{args.files_prefix}.tmp.fq").get_kmer_counts(args.files_prefix,threads=args.threads,max_mem=args.ram,counter = args.kmer_counter)
        else:
            bam = Bam(args.bam,args.files_prefix,"illumina")
            run_cmd(f"samtools fastq {bam.bam_file} > {args.files_prefix}.tmp.fq")
            kmer_dump = Fastq(f"{args.files_prefix}.tmp.fq").get_kmer_counts(args.files_prefix,threads=args.threads,max_mem=args.ram,counter = args.kmer_counter)
    elif "fasta" in vars(args) and args.fasta:
        kmer_dump = Fasta(args.fasta).get_kmer_counts(args.files_prefix,threads=args.threads,max_mem=args.ram,counter = args.kmer_counter)
    if "output_kmer_counts" not in vars(args):
        args.output_kmer_counts = None
    else:
        args.output_kmer_counts = f"{args.prefix}.kmers.txt"  if args.output_kmer_counts else False
    species = kmer_dump.get_taxonomic_support(conf['kmers'],args.output_kmer_counts)
    return species

def get_bam_file(args):
    ### Create bam file if fastq has been supplied ###
    if args.bam is None:
        if args.read1 and args.read2 and args.no_trim:
            # Paired + no trimming
            fastq_obj = Fastq(args.read1,args.read2)
        elif args.read1 and args.read2 and not args.no_trim:
            # Paired + trimming
            untrimmed_fastq_obj = Fastq(args.read1,args.read2)
            fastq_obj = untrimmed_fastq_obj.trim(args.files_prefix,threads=args.threads)
        elif args.read1 and not args.read2 and args.no_trim:
            # Unpaired + no trimming
            fastq_obj = Fastq(args.read1,args.read2)
        elif args.read1 and not args.read2 and not args.no_trim:
            # Unpaired + trimming
            untrimmed_fastq_obj = Fastq(args.read1)
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
    if args.platform in ("nanopore","pacbio"):
        args.mapper = "minimap2"
        if args.caller=="gatk":
            args.caller = "freebayes"
        args.no_trim=True
        args.run_delly = True
    else:
        if "no_delly" in vars(args):
            args.run_delly = False if args.no_delly else True

    return args

