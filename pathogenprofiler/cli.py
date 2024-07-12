from .fastq import Fastq
from .utils import run_cmd, cmd_out
from .bam import Bam
from .db import get_db
from .fasta import Fasta, Paf
from .profiler import vcf_profiler, bam_barcoder, vcf_barcoder
import logging
from typing import List, Union
import argparse
from .models import Variant, DrVariant, Gene, DrGene, SpeciesPrediction, Species, BarcodeResult
from .mutation_db import MutationDB
from .vcf import Vcf
from .sanity import check_bam_for_rg, check_vcf_chrom_match

def get_variant_filters(args):
    filters = {}
    for f in [
        'depth','af','strand','sv_depth','sv_af','sv_len'
    ]:
        if not hasattr(args,f):
            continue
        vals = getattr(args,f).split(",")
        if len(vals)==2:
            filters[f+"_hard"] = float(vals[0]) if "." in vals[0] else int(vals[0])
            filters[f+"_soft"] = float(vals[1]) if "." in vals[1] else int(vals[1])
        elif len(vals)==1:
            if vals[0]=="-":
                filters[f+"_hard"] = None
                filters[f+"_soft"] = None   
            else:
                filters[f+"_hard"] = 0.0
                filters[f+"_soft"] = float(vals[0]) if "." in vals[0] else int(vals[0])
        else:
            logging.error(f"Error parsing {f} ({getattr(args,f)}) filter. Please provide 1 or 2 values separated by a comma")
            exit()
    return filters

def get_fastq_qc(args: argparse.Namespace):
    fastq = Fastq(args.read1,args.read2)
    return fastq.get_qc()

def run_bam_qc(args: argparse.Namespace):
    bam = Bam(args.bam, args.files_prefix, platform=args.platform, threads=args.threads)
    qc = bam.get_bam_qc(
        bed_file=args.conf["bed"],
        ref_file=args.conf["ref"],
        depth_cutoff=args.conf['variant_filters']['depth_soft'],
    )
    mutation_db = MutationDB(args.conf["json_db"])
    mutation_db.annotate_missing_positions(qc.missing_positions)
    return qc

def run_fasta_qc(args: argparse.Namespace):
    fasta = Fasta(args.fasta)
    qc = fasta.get_fasta_qc()
    if hasattr(args,'paf') and args.paf:
        paf = Paf(args.paf)
        qc.target_qc = paf.get_target_qc(
            bed_file=args.conf["bed"]
        )
    return qc

def run_vcf_qc(args: argparse.Namespace):
    vcf = Vcf(args.vcf)
    qc = vcf.get_vcf_qc()
    
    return qc

def get_resistance_db_from_species_prediction(args: argparse.Namespace,species_prediction:SpeciesPrediction):
    logging.debug("Attempting to load db with species prediction")
    number_of_species = len(set([s.species for s in species_prediction.species]))
    if number_of_species==1:
        return get_db(args.software_name,species_prediction.species[0].species.replace(" ","_")) 
    else:
        return None
    

    
def test_vcf_for_lofreq(vcf_file):
    lofreq = False
    for l in cmd_out(f"bcftools view -h {vcf_file}"):
        if "source=lofreq call" in l:
            lofreq = True
    return lofreq

def get_input_data_source(args):
    if args.read1:
        return "fastq"
    elif args.bam:
        return "bam"
    elif args.fasta:
        return "fasta"
    elif args.vcf:
        return "vcf"
    else:
        return None
    
# def get_qc(args: argparse.Namespace):
#     if args.qc:
#         return True
#     else:
#         return False
    


def process_args(args: argparse.Namespace) -> None:
    """
    Process the input arguments and return a Namespace object
    
    Arguments
    ---------
    args : argparse.Namespace
        The input arguments

    Returns
    -------
    argparse.Namespace
        The processed arguments
    
    Examples
    --------
    >>> from pathogenprofiler import cli
    >>> import argparse
    >>> args = argparse.Namespace()
    >>> args.software_name = "ntm-profiler"
    >>> args.resistance_db = "Mycobacterium_abscessus"
    >>> args.platform = "illumina"
    >>> args.no_samclip = False
    >>> args.no_coverage_qc = False
    >>> args.read1 = "test.fastq.gz"
    >>> args = cli.process_args(args)
    >>> args.software_name
    'ntm-profiler'
    """
    args = set_platform_params(args)
    args.samclip = True if not args.no_samclip else False
    args.coverage_qc = True if not args.no_coverage_qc else False
    args.data_source = get_input_data_source(args)
    if 'conf' in vars(args) and args.conf:
        args.conf['variant_filters'] = get_variant_filters(args)

def get_vcf_from_bam(args: argparse.Namespace):
    conf = args.conf
    ### Create bam object and call variants ###
    bam = Bam(args.bam, args.files_prefix, platform=args.platform, threads=args.threads)
    if args.call_whole_genome:
        wg_vcf_obj = bam.call_variants(conf["ref"], caller=args.caller, filters = conf['variant_filters'], threads=args.threads, calling_params=args.calling_params, samclip = args.samclip, cli_args=vars(args))
        vcf_obj = wg_vcf_obj
        # TODO optional?
        # vcf_obj = wg_vcf_obj.view_regions(conf["bed"])
    else:
        vcf_obj = bam.call_variants(conf["ref"], caller=args.caller, filters = conf['variant_filters'], bed_file=conf["bed"], threads=args.threads, calling_params=args.calling_params, samclip = args.samclip, cli_args=vars(args))

    ### Run delly if specified ###
    final_target_vcf_file = args.files_prefix+".targets.vcf.gz"
    if not args.no_delly:
        delly_vcf_obj = bam.run_delly(conf['ref'],conf['bed'])
        if delly_vcf_obj is not None:
            run_cmd("bcftools index %s" % delly_vcf_obj.filename)
            run_cmd("bcftools concat %s %s | bcftools sort -Oz -o %s" % (vcf_obj.filename,delly_vcf_obj.filename,final_target_vcf_file))
        else:
            return vcf_obj.filename
            #run_cmd("mv %s %s" % (vcf_obj.filename, final_target_vcf_file))
    else:
        return vcf_obj.filename
        #run_cmd("mv %s %s" % (vcf_obj.filename, final_target_vcf_file))
    
    return final_target_vcf_file

def get_vcf_file(args: argparse.Namespace):
    if args.vcf:
        args.vcf = args.vcf
    elif args.fasta:
        args.paf = Fasta(args.fasta).align_to_ref(args.conf["ref"],args.files_prefix)
        args.vcf = Paf(args.paf).get_ref_variants(args.conf["ref"], args.prefix, args.files_prefix)
    elif args.bam:
        args.vcf = get_vcf_from_bam(args)

def run_barcoder(args: argparse.Namespace) -> List[BarcodeResult]:
    if args.data_source in ('fastq', 'bam'):
        if not args.bam:
            quit()
        else:
            barcode_result = bam_barcoder(args)
    elif args.data_source == 'fasta':
        barcode_result = vcf_barcoder(args)
    elif args.data_source == 'vcf':
        barcode_result = vcf_barcoder(args)

    return barcode_result

def run_profiler(args) -> List[Union[Variant,DrVariant,Gene,DrGene]]:
    if args.read1 or args.bam:
        args.bam = get_bam_file(args)
        get_vcf_file(args)
        args.supplementary_bam = args.bam
        
    elif args.fasta:
        get_vcf_file(args)

    elif args.vcf:
        if test_vcf_for_lofreq(args.vcf):
            tmp_vcf_file = f"{args.files_prefix}.tmp.vcf.gz"
            run_cmd(f"bcftools view {args.vcf} | modify_lofreq_vcf.py --sample {args.prefix} | bcftools view -Oz -o {tmp_vcf_file}")
            args.vcf = tmp_vcf_file
    # check_vcf_chrom_match(args.vcf,args.conf["ref"])
    annotated_variants = vcf_profiler(args)

    return annotated_variants

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
        check_bam_for_rg(args.bam)
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

    if args.fasta:
        args.depth = "0,0"
        args.strand = "0,0"
    return args

def get_sourmash_hit(args):
    args.species_conf = get_db(args.software_name,args.species_db)
    if args.read1:
        if args.read2:
            fastq = Fastq(args.read1,args.read2)
        else:
            fastq = Fastq(args.read1)
        sourmash_sig = fastq.sourmash_sketch(args.files_prefix)
    elif args.fasta:
        fasta = Fasta(args.fasta)
        sourmash_sig = fasta.sourmash_sketch(args.files_prefix)
    elif args.bam:
        run_cmd(f"samtools fastq {args.bam} > {args.files_prefix}.tmp.fastq")
        fq_file = f"{args.files_prefix}.tmp.fastq"
        fastq = Fastq(fq_file)
        sourmash_sig = fastq.sourmash_sketch(args.files_prefix)

    sourmash_sig = sourmash_sig.gather(args.species_conf["sourmash_db"],args.species_conf["sourmash_db_info"],intersect_bp=500000,f_match_threshold=0.1)
    result =  []

    if len(sourmash_sig)>0:
        result = sourmash_sig
    
    return result



def set_species(args: argparse.Namespace) -> SpeciesPrediction:
    """
    Get a SpeciesPrediction object based on the input arguments
    
    Arguments
    ---------
    args : argparse.Namespace
        The input arguments

    Returns
    -------
    SpeciesPrediction
        A SpeciesPrediction object
    
    Examples
    --------
    >>> from pathogenprofiler import cli
    >>> import argparse
    >>> args = argparse.Namespace()
    >>> args.software_name = "ntm-profiler"
    >>> args.resistance_db = "Mycobacterium_abscessus"
    >>> species_prediction = cli.set_species(args)
    >>> species_prediction.species
    [Species(species='Mycobacterium abscessus', prediction_info=None)]
    """
    conf = get_db(args.software_name,args.resistance_db)
    species = Species(
        species=conf["species"]
    )
    data = {
        "species":[species],
        "prediction_method":"user_defined",
    }
    return SpeciesPrediction(**data)
        
        


def get_sourmash_species_prediction(args: argparse.Namespace) -> SpeciesPrediction:
    """
    Get a SpeciesPrediction object based on prediction using sourmash

    Arguments
    ---------
    args : argparse.Namespace
        The input arguments

    Returns
    -------
    SpeciesPrediction
        A SpeciesPrediction object

    """
    conf = get_db(args.software_name,args.species_db)
    sourmash_species_prediction = get_sourmash_hit(args)
    species = []
    for obj in sourmash_species_prediction:
        species.append(
            Species(
                species=obj['species'],
                prediction_method='sourmash',
                prediction_info=obj,
            )
        )
    return SpeciesPrediction(
        species=species,
        prediction_method='sourmash',
        species_db = conf['version']
    )



