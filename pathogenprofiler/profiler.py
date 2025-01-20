from .utils import run_cmd, bed2gene_lookup
from .bam import Bam
from .barcode import barcode
from .mutation_db import db_compare
from .models import BarcodeResult, DrVariant, Variant, Gene, DrGene
from .vcf import Vcf
import os
from .models import Variant
from typing import List, Union
import argparse
import logging

def bam_barcoder(args: argparse.Namespace) -> List[BarcodeResult]:
    conf = args.conf
    if 'barcode' not in conf:
        return []
    bam = Bam(args.bam, args.files_prefix, platform=args.platform, threads=args.threads)
    if not hasattr(args,'barcode_snps'):
        args.barcode_snps = None
    barcode_mutations = bam.get_bed_gt(conf["barcode"],conf["ref"], caller=args.caller,platform=args.platform)  
    stdev_cutoff = args.barcode_stdev if hasattr(args,'barcode_stdev') else None      
    barcode_assignment = barcode(barcode_mutations,conf["barcode"],args.barcode_snps,stdev_cutoff=stdev_cutoff)
    return barcode_assignment

def vcf_barcoder(args: argparse.Namespace) -> List[BarcodeResult]:
    conf = args.conf
    if 'barcode' not in conf:
        return []
    vcf = Vcf(args.vcf)
    barcode_mutations = vcf.get_bed_gt(conf["barcode"],conf["ref"])        
    if not hasattr(args,'barcode_snps'):
        args.barcode_snps = None
    stdev_cutoff = args.barcode_stdev if hasattr(args,'barcode_stdev') else None      
    barcode_assignment = barcode(barcode_mutations,conf["barcode"],args.barcode_snps,stdev_cutoff=stdev_cutoff)
    return barcode_assignment




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
    if args.caller == 'freebayes-haplotype':
        args.supplementary_bam = None
    annotated_variants = vcf_variant_profiler(
        conf=conf, 
        prefix=args.files_prefix, 
        vcf_file=args.vcf, 
        bam_for_phasing=args.supplementary_bam,
        db_dir=args.db_dir if hasattr(args,'db_dir') else None
    )
    return annotated_variants
    
 
    



def vcf_variant_profiler(conf: dict, prefix: str, vcf_file: str, bam_for_phasing: str = None, db_dir: str = None) -> List[Union[Variant,Gene]]:
    vcf_targets_file = "%s.targets_for_profile.vcf.gz" % prefix
    if not vcf_is_indexed(vcf_file):
        run_cmd("bcftools index %s" % vcf_file)
    run_cmd("bcftools view -R %s %s -Oz -o %s" % (conf["bed"],vcf_file,vcf_targets_file))
    vcf_obj = Vcf(vcf_targets_file)
    vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms=conf.get("chromosome_conversion",None),bam_for_phasing=bam_for_phasing, db_dir=db_dir)
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
