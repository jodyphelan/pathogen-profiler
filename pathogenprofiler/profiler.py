from .utils import run_cmd
from .bam import Bam
from .barcode import barcode, db_compare
from .vcf import Vcf,DellyVcf
from .fasta import Fasta
import statistics as stats
import os
import logging



def bam_profiler(conf, bam_file, prefix, platform, caller, threads=1, no_flagstat=False, run_delly=True, calling_params=None, delly_vcf_file=None, min_depth = 10, samclip=False, variant_annotations = False, call_wg=False,coverage_tool="bedtools"):
    logging.warning("Please ensure that this BAM was made using the same reference as in the database. If you are not sure what reference was used it is best to remap the reads.")

    ### Put user specified arguments to lower case ###
    platform = platform.lower()
    caller = caller.lower()

    ### Set caller to freebayes if platform is nanopre and wrong caller has been used ###
    if platform in ("nanopore","pacbio"):
        run_delly = False
        if caller=="gatk":
            caller = "freebayes"

    ### Create bam object and call variants ###
    bam = Bam(bam_file, prefix, platform=platform, threads=threads)
    if call_wg:
        wg_vcf_obj = bam.call_variants(conf["ref"], caller=caller, filters = conf['variant_filters'], threads=threads, calling_params=calling_params, samclip = samclip, min_dp=min_depth)
        vcf_obj = wg_vcf_obj.view_regions(conf["bed"])
    else:
        vcf_obj = bam.call_variants(conf["ref"], caller=caller, filters = conf['variant_filters'], bed_file=conf["bed"], threads=threads, calling_params=calling_params, samclip = samclip, min_dp=min_depth)
    if variant_annotations:
        vcf_obj = vcf_obj.add_annotations(conf["ref"],bam.bam_file)
    else:
        ann_vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None))
    ann = ann_vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types = ["ablation","upstream","synonymous","noncoding"])

    # bam.get_region_qc(conf["bed"],conf["ref"],min_dp=min_depth)

    results = {}
    
    ### Get % and num reads mapping ###
    if not no_flagstat:
        results['qc'] = {}
        bam.calculate_bamstats()
        results['qc']['pct_reads_mapped'] = bam.pct_reads_mapped
        results['qc']['num_reads_mapped'] = bam.mapped_reads
        results['qc']['region_qc'] = bam.get_region_qc(bed_file=conf['bed'],cutoff=min_depth)
        results['qc']['region_median_depth'] = stats.median([x['median_depth'] for x in results['qc']['region_qc']])
        results["qc"]["missing_positions"] = bam.get_missing_genomic_positions(cutoff=min_depth)
        if 'amplicon' not in conf or conf['amplicon']==False:
            results['qc']['genome_median_depth'] = bam.get_median_depth(ref_file=conf['ref'],software=coverage_tool)
        
    results["variants"]  = ann

    if "barcode" in conf:
        if platform in ("nanopore","pacbio"):
            mutations = bam.get_bed_gt(conf["barcode"],conf["ref"], caller="bcftools",platform=platform)
        else:
            mutations = bam.get_bed_gt(conf["barcode"],conf["ref"], caller=caller,platform=platform)

            
        results["barcode"] = barcode(mutations,conf["barcode"])
    

    ### Run delly if specified ###
    if run_delly:
        if delly_vcf_file:
            delly_vcf_obj = Vcf(delly_vcf_file)
        else:
            delly_vcf_obj = bam.run_delly(conf['bed'])
        if delly_vcf_obj is not None:
            results["delly"] = "success"
            # delly_vcf_obj = delly_vcf_obj.get_robust_calls(prefix,conf["bed"])
            ann_vcf_obj = delly_vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None),split_indels=False)
            results["variants"].extend(ann_vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types=["ablation"] ))
        else:
            results["delly"] = "fail"

    ### Compare variants to database ###
    results = db_compare(db=conf["json_db"], mutations=results)
    return results


def fasta_profiler(conf, prefix, filename):
    fasta = Fasta(filename)
    wg_vcf_file = fasta.get_ref_variants(conf["ref"], prefix)
    wg_vcf_obj = Vcf(wg_vcf_file)
    vcf_file = prefix+".targets.vcf.gz"
    run_cmd("bcftools view -c 1 %s -Oz -o %s -T %s" % (wg_vcf_file,vcf_file,conf['bed']))
    vcf_obj = Vcf(vcf_file)
    vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None))
    ann = vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types = ["ablation","upstream","synonymous","noncoding"])

    results = {
        "variants":ann,
        "qc":{
            "gene_coverage":fasta.get_aln_coverage(conf['bed']),
            "num_reads_mapped":"NA"}
    }

    if "barcode" in conf:
        mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])
        results["barcode"] = barcode(mutations,conf["barcode"])
    results = db_compare(db=conf["json_db"], mutations=results)
    return results

def vcf_is_indexed(vcf_file):
    if os.path.exists(vcf_file+".tbi"):
        return True
    elif os.path.exists(vcf_file+".csi"):
        return True
    else:
        return False

def vcf_profiler(conf, prefix, sample_name, vcf_file,delly_vcf_file):
    vcf_targets_file = "%s.targets.vcf.gz" % prefix
    if not vcf_is_indexed(vcf_file):
        run_cmd("bcftools index %s" % vcf_file)
    run_cmd("bcftools view -R %s %s -Oz -o %s" % (conf["bed"],vcf_file,vcf_targets_file))
    vcf_obj = Vcf(vcf_targets_file)
    vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None))
    ann = vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types = ["ablation","upstream","synonymous","noncoding"])
    results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
    results["variants"]  = ann

    if delly_vcf_file:
        delly_vcf_obj = DellyVcf(delly_vcf_file)
        # delly_vcf_obj = delly_vcf_obj.get_robust_calls(prefix,conf["bed"])
        ann_vcf_obj = delly_vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None),split_indels=False)
        results["variants"].extend(ann_vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types=["ablation"]))
 
    if "barcode" in conf:
        mutations = Vcf(vcf_file).get_bed_gt(conf["barcode"], conf["ref"])
        results["barcode"] = barcode(mutations,conf["barcode"])
    results = db_compare(db=conf["json_db"], mutations=results)

    return results
