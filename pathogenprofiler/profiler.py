from .utils import run_cmd, bed2gene_lookup
from .bam import Bam
from .barcode import barcode
from .mutation_db import db_compare
from .models import BarcodeResult, DrVariant, Variant, Gene, DrGene
from .vcf import Vcf,DellyVcf
from .fasta import Fasta
import statistics as stats
import os
import logging
from .models import Variant
from typing import List, Union
import argparse
from .rules import apply_rules

def bam_barcoder(args: argparse.Namespace) -> List[BarcodeResult]:
    conf = args.conf
    bam = Bam(args.bam, args.files_prefix, platform=args.platform, threads=args.threads)
    barcode_mutations = bam.get_bed_gt(conf["barcode"],conf["ref"], caller=args.caller,platform=args.platform)        
    barcode_assignment = barcode(barcode_mutations,conf["barcode"])
    return barcode_assignment

def vcf_barcoder(args: argparse.Namespace) -> List[BarcodeResult]:
    conf = args.conf
    vcf = Vcf(args.vcf)
    barcode_mutations = vcf.get_bed_gt(conf["barcode"],conf["ref"])        
    barcode_assignment = barcode(barcode_mutations,conf["barcode"])
    return barcode_assignment

def bam_profiler(args: argparse.Namespace) -> List[Union[Variant,DrVariant,Gene,DrGene]]:
    logging.warning("Please ensure that this BAM was made using the same reference as in the database. If you are not sure what reference was used it is best to remap the reads.")
    conf = args.conf
    ### Create bam object and call variants ###
    bam = Bam(args.bam, args.files_prefix, platform=args.platform, threads=args.threads)
    if args.call_whole_genome:
        wg_vcf_obj = bam.call_variants(conf["ref"], caller=args.caller, filters = conf['variant_filters'], threads=args.threads, calling_params=args.calling_params, samclip = args.samclip)
        vcf_obj = wg_vcf_obj.view_regions(conf["bed"])
    else:
        vcf_obj = bam.call_variants(conf["ref"], caller=args.caller, filters = conf['variant_filters'], bed_file=conf["bed"], threads=args.threads, calling_params=args.calling_params, samclip = args.samclip)

    ### Run delly if specified ###
    if not args.no_delly:
        final_target_vcf_file = args.files_prefix+".targets.vcf.gz"
        delly_vcf_obj = bam.run_delly(conf['bed'])
        if delly_vcf_obj is not None:
            run_cmd("bcftools index %s" % delly_vcf_obj.filename)
            run_cmd("bcftools concat %s %s | bcftools sort -Oz -o %s" % (vcf_obj.filename,delly_vcf_obj.filename,final_target_vcf_file))
        else:
            run_cmd("mv %s %s" % (vcf_obj.filename, final_target_vcf_file))
    else:
        run_cmd("mv %s %s" % (vcf_obj.filename, final_target_vcf_file))
    
            


    ### Annotate variants ###
    annoted_variants = vcf_variant_profiler(conf, args.files_prefix, final_target_vcf_file)
    return annoted_variants
    
    

# def fasta_profiler(args: argparse.Namespace) -> List[Union[Variant,DrVariant,Gene,DrGene]]:
#     conf = args.conf
#     fasta = Fasta(args.fasta)
#     wg_vcf_file = fasta.get_ref_variants(conf["ref"], args.prefix)
#     quit()
#     annotated_vaiants = vcf_variant_profiler(conf, args.files_prefix, wg_vcf_file)

#     return annotated_vaiants    

def vcf_is_indexed(vcf_file: str) -> bool:
    if os.path.exists(vcf_file+".tbi"):
        return True
    elif os.path.exists(vcf_file+".csi"):
        return True
    else:
        return False

def vcf_profiler(args: argparse.Namespace) -> List[Union[Variant,DrVariant,Gene,DrGene]]:
    conf = args.conf
    vcf_targets_file = "%s.targets.vcf.gz" % args.files_prefix
    if not vcf_is_indexed(args.vcf):
        run_cmd("bcftools index %s" % args.vcf)
    # run_cmd("bcftools view -R %s %s -Oz -o %s" % (conf["bed"],args.vcf,vcf_targets_file))
    annotated_variants = vcf_variant_profiler(conf, args.files_prefix, args.vcf)
    return annotated_variants
    
 
    



def vcf_variant_profiler(conf: dict, prefix: str, vcf_file: str) -> List[Union[Variant,Gene]]:
    vcf_targets_file = "%s.targets_for_profile.vcf.gz" % prefix
    if not vcf_is_indexed(vcf_file):
        run_cmd("bcftools index %s" % vcf_file)
    run_cmd("bcftools view -R %s %s -Oz -o %s" % (conf["bed"],vcf_file,vcf_targets_file))
    vcf_obj = Vcf(vcf_targets_file)
    vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None))
    variants = vcf_obj.load_ann(conf['variant_filters'],bed_file=conf["bed"],keep_variant_types = ["ablation","upstream","synonymous","noncoding"])

    # compare against database of variants
    annotated_variants = db_compare(db=conf["json_db"], variants=variants)
    
    

    # set the default consequence
    for obj in annotated_variants:
        # check if it is a variant or a gene class
        if isinstance(obj, Variant):
            obj.set_default_csq()

    # Add gene names based on the bed file
    gene_id2name = bed2gene_lookup(conf["bed"])
    for obj in annotated_variants:
        obj.set_gene_name(gene_id2name)
    
    # if 'rules' in conf:
    #     rules_applied = apply_rules(conf['rules'], annotated_variants)

    # # Convert variant objects to DrVariant if they cause resistance
    # for var in annotated_variants:
    #     var.convert_to_dr_element()

    return annotated_variants
