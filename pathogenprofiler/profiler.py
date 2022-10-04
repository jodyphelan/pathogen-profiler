from .utils import infolog, run_cmd, debug
from .bam import bam
from .barcode import barcode, db_compare
from .vcf import vcf,delly_bcf
from .fasta import fasta
from .variant_set import variant_set
import os




def bam_profiler(conf, bam_file, prefix, platform, caller, dir, threads=1, no_flagstat=False, run_delly=True, calling_params=None, delly_vcf_file=None, run_coverage=True, coverage_fraction_threshold=0, min_depth = 10, missing_cov_threshold=10, min_af=0.1, samclip=False, variant_annotations = False, call_wg=False,coverage_tool="bedtools"):
    infolog("Using %s\n\nPlease ensure that this BAM was made using the same reference as in the database.\nIf you are not sure what reference was used it is best to remap the reads." % bam_file)

    ### Put user specified arguments to lower case ###
    platform = platform.lower()
    caller = caller.lower()

    ### Set caller to freebayes if platform is nanopre and wrong caller has been used ###
    if platform in ("nanopore","pacbio"):
        run_delly = False
        if caller=="gatk":
            caller = "freebayes"

    ### Create bam object and call variants ###
    bam_obj = bam(bam_file, prefix, platform=platform, threads=threads)
    if call_wg:
        wg_vcf_obj = bam_obj.call_variants(conf["ref"], caller=caller, threads=threads, calling_params=calling_params, samclip = samclip, min_dp=min_depth)
        vcf_obj = wg_vcf_obj.view_regions(conf["bed"])
    else:
        vcf_obj = bam_obj.call_variants(conf["ref"], caller=caller, bed_file=conf["bed"], threads=threads, calling_params=calling_params, samclip = samclip, min_dp=min_depth)
    if variant_annotations:
        vcf_obj = vcf_obj.add_annotations(conf["ref"],bam_obj.bam_file)
    else:
        ann_vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None))
    ann = ann_vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types = ["upstream","synonymous","noncoding"],min_af=min_af)


    ### Get % and num reads mapping ###
    if no_flagstat:
        bam_obj.pct_reads_mapped = "NA"
        bam_obj.num_reads_mapped = "NA"
        bam_obj.median_coverage = "NA"
    else:
        bam_obj.flagstat()
        bam_obj.get_median_coverage(ref_file=conf["ref"],software=coverage_tool)

    ### Put results into a dictionary ###
    results = {
        "variants":[],
        "qc":{
            "pct_reads_mapped":bam_obj.pct_reads_mapped,
            "num_reads_mapped":bam_obj.num_reads_mapped,
            "median_coverage":bam_obj.median_coverage
        }
    }

    if run_coverage:
        results["qc"]["gene_coverage"] = bam_obj.get_region_coverage(conf["bed"], fraction_threshold= coverage_fraction_threshold)
        results["qc"]["missing_positions"] = bam_obj.get_missing_genomic_positions(cutoff=missing_cov_threshold)
    results["variants"]  = ann

    if "barcode" in conf:
        if platform in ("nanopore","pacbio"):
            mutations = bam_obj.get_bed_gt(conf["barcode"],conf["ref"], caller="bcftools",platform=platform)
        else:
            mutations = bam_obj.get_bed_gt(conf["barcode"],conf["ref"], caller=caller,platform=platform)

            
        results["barcode"] = barcode(mutations,conf["barcode"])
    

    ### Run delly if specified ###
    if run_delly:
        if delly_vcf_file:
            delly_vcf_obj = delly_bcf(delly_vcf_file)
        else:
            delly_vcf_obj = bam_obj.run_delly()
        if delly_vcf_obj is not None:
            results["delly"] = "success"
            delly_vcf_obj = delly_vcf_obj.get_robust_calls(prefix,conf["bed"])
            ann_vcf_obj = delly_vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None),split_indels=False)
            results["variants"].extend(ann_vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types=["ablation"] ))
        else:
            results["delly"] = "fail"

    if call_wg:
        var_set = variant_set(wg_vcf_obj.filename,exclude_bed="/home/jody/refgenome/excluded_loci.bed")
        results["close_samples"] = var_set.get_close_samples(os.path.join(dir,"results"))
    ### Compare variants to database ###
    results = db_compare(db=conf["json_db"], mutations=results)
    return results


def fasta_profiler(conf, prefix, filename):
    fasta_obj = fasta(filename)
    wg_vcf_file = fasta_obj.get_ref_variants(conf["ref"], prefix)
    wg_vcf_obj = vcf(wg_vcf_file)
    vcf_file = prefix+".targets.vcf.gz"
    run_cmd("bcftools view -c 1 %s -Oz -o %s -T %s" % (wg_vcf_file,vcf_file,conf['bed']))
    vcf_obj = vcf(vcf_file)
    vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None))
    ann = vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types = ["upstream","synonymous","noncoding"])

    results = {
        "variants":ann,
        "qc":{
            "gene_coverage":fasta_obj.get_aln_coverage(conf['bed']),
            "num_reads_mapped":"NA"}
    }

    if "barcode" in conf:
        mutations = wg_vcf_obj.get_bed_gt(conf["barcode"], conf["ref"])
        results["barcode"] = barcode(mutations,conf["barcode"])
    results = db_compare(db=conf["json_db"], mutations=results)
    return results

def vcf_profiler(conf, prefix, sample_name, vcf_file,delly_vcf_file):
    vcf_targets_file = "%s.targets.vcf.gz" % prefix
    run_cmd("bcftools view -T %s %s -Oz -o %s" % (conf["bed"],vcf_file,vcf_targets_file))
    vcf_obj = vcf(vcf_targets_file)
    vcf_obj = vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None))
    ann = vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types = ["upstream","synonymous","noncoding"])
    results = {"variants":[],"missing_pos":[],"qc":{"pct_reads_mapped":"NA","num_reads_mapped":"NA"}}
    results["variants"]  = ann

    if delly_vcf_file:
        delly_vcf_obj = delly_bcf(delly_vcf_file)
        delly_vcf_obj = delly_vcf_obj.get_robust_calls(prefix,conf["bed"])
        ann_vcf_obj = delly_vcf_obj.run_snpeff(conf["snpEff_db"],conf["ref"],conf["gff"],rename_chroms= conf.get("chromosome_conversion",None),split_indels=False)
        results["variants"].extend(ann_vcf_obj.load_ann(bed_file=conf["bed"],keep_variant_types=["ablation"]))
 
    if "barcode" in conf:
        mutations = vcf(vcf_file).get_bed_gt(conf["barcode"], conf["ref"])
        results["barcode"] = barcode(mutations,conf["barcode"])
    results = db_compare(db=conf["json_db"], mutations=results)

    return results
