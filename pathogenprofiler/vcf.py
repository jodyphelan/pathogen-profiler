from .utils import TempFilePrefix, run_cmd, cmd_out,add_arguments_to_self, index_bcf,tabix, load_bed
import logging
from .fasta import Fasta
from collections import defaultdict
import re
import sys
import os.path
import json
import pysam
from .models import Variant,Consequence, VcfQC, GenomePosition
from typing import Optional, List, Tuple, Dict

def get_default_snpeff_config():
    share_dir = os.path.join(sys.base_prefix, 'share')
    # check folder exists
    if os.path.isdir(share_dir):
        snp_eff_dirs = [d for d in os.listdir(share_dir) if d.startswith('snpeff-')]
        if len(snp_eff_dirs) == 1:
            snp_eff_dir = snp_eff_dirs[0]
            return os.path.join(share_dir, snp_eff_dir, 'snpEff.config')
    
def get_custom_snpeff_config(software_name: str):
    # check folder exists
    share_dir = os.path.join(sys.base_prefix, 'share')
    software_dir = os.path.join(share_dir, software_name)
    if os.path.isdir(software_dir):
        config = os.path.join(software_dir, 'snpEff.config')
        if os.path.isfile(config):
            return config

def get_stand_support(var: pysam.VariantRecord,alt: str) -> Tuple[int,int]:
    alt_index = list(var.alts).index(alt)
    forward_support = None
    reverse_support = None
    if 'SAF' in var.info:
        forward_support = var.info['SAF'][alt_index]
        reverse_support = var.info['SAR'][alt_index]
    elif 'ADF' in var.samples[0]:
        forward_support = var.samples[0]['ADF'][alt_index+1]
        reverse_support = var.samples[0]['ADR'][alt_index+1]
    elif 'SB' in var.samples[0]:
        logging.debug(var.samples[0]['SB'])
        forward_support = var.samples[0]['SB'][1*2]
        reverse_support = var.samples[0]['SB'][1*2+1]
    elif 'DP4' in var.info:
        forward_support = var.info['DP4'][1*2]
        reverse_support = var.info['DP4'][1*2+1]
    else:
        forward_support = None
        reverse_support = None
    
    logging.debug(f"Forward support: {forward_support}, Reverse support: {reverse_support}")
    return forward_support, reverse_support

def get_sv_ad(var):
    return [
        var.samples[0]['DR']+var.samples[0]['RR'],
        var.samples[0]['DV']+var.samples[0]['RV']
    ]
class Vcf:
    def __init__(self,filename,prefix=None,threads=1):
        self.samples = []
        add_arguments_to_self(self,locals())
        if prefix==None:
            if filename[-4:] == ".bcf":
                self.prefix = filename[:-4]
            elif filename[-7:] == ".vcf.gz":
                self.prefix = filename[:-7]
            elif filename[-8:] == ".gvcf.gz":
                self.prefix = filename[:-8]
            elif filename[-4:] == ".vcf":
                self.prefix = filename[:-4]
            else:
                self.prefix = filename
        else:
            self.prefix = prefix
        if filename[-4:] == ".bcf":
            index_bcf(filename,threads)
        else:
            tabix(filename,threads)
        self.vcf_dir = "/".join(os.path.abspath(filename).split("/")[:-1])
        for l in cmd_out("bcftools query -l %(filename)s" % vars(self)):
            self.samples.append(l.rstrip())
        self.nsamples = len(self.samples)
        header = "\n".join([l.strip() for l in cmd_out("bcftools view -h %(filename)s" % vars(self))])
        if "bcftools_callCommand" in header:
            self.caller = "bcftools"
        elif "source=lofreq call" in header:
            self.caller = "lofreq"
        elif "GATKCommandLine=<ID=HaplotypeCaller" in header:
            self.caller = "gatk"
        elif 'source="Pilon' in header:
            self.caller = 'pilon'
        elif 'source=freeBayes' in header:
            self.caller = 'freebayes'
        else:
            self.caller = 'Unknown'
        

    def view_regions(self,bed_file):
        self.bed_file = bed_file
        self.newfile = "%(prefix)s.region_subset.vcf.gz" % vars(self)
        run_cmd("bcftools view -R %(bed_file)s %(filename)s -Oz -o %(newfile)s" % vars(self))
        return Vcf(self.newfile)


    def run_snpeff(self, db: str,ref_file: str,gff_file: str,rename_chroms: dict = None, split_indels: bool=True, bam_for_phasing: str=None, db_dir: str=None):
        logging.info("Running snpEff")
        add_arguments_to_self(self,locals())
        self.vcf_csq_file = self.prefix+".csq.vcf.gz"
        self.rename_cmd = f"rename_vcf_chrom.py --source {' '.join(rename_chroms['source'])} --target {' '.join(rename_chroms['target'])} |" if rename_chroms else ""
        self.re_rename_cmd = f"| rename_vcf_chrom.py --source {' '.join(rename_chroms['target'])} --target {' '.join(rename_chroms['source'])}" if rename_chroms else ""
        
        if db_dir:
            self.snpeff_config_opt = f'-config {db_dir}/snpeff/snpEff.config'
        else:
            self.snpeff_config_opt = f'-config {get_default_snpeff_config()}'

        if split_indels:
            with TempFilePrefix() as tmp:
                self.tmp_file1 = f"{tmp}.1.vcf.gz"
                self.tmp_file2 = f"{tmp}.2.vcf.gz"
                self.tmp_file3 = f"{tmp}.3.vcf.gz"
                if bam_for_phasing:
                    self.phasing_bam = f"--bam {bam_for_phasing}"
                else:
                    self.phasing_bam = ""
                run_cmd("bcftools view -a %(filename)s | bcftools view -c 1 | realign_tandem_deletions.py - %(ref_file)s %(gff_file)s - |  combine_vcf_variants.py --ref %(ref_file)s --gff %(gff_file)s %(phasing_bam)s | %(rename_cmd)s snpEff ann %(snpeff_config_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools sort -Oz -o %(vcf_csq_file)s" % vars(self))
                
        else :
            run_cmd("bcftools view %(filename)s | %(rename_cmd)s snpEff ann %(snpeff_config_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools view -Oz -o %(vcf_csq_file)s" % vars(self))
        return Vcf(self.vcf_csq_file,self.prefix)



    def load_ann(
        self,
        filter_params: dict,
        max_promoter_length: int = 1000, 
        bed_file: Optional[str] = None,
        exclude_variant_types: Optional[List[str]] = None,
        keep_variant_types: Optional[List[str]] = None
    ) -> List[Variant]:
        logging.info("Loading snpEff annotations")
        filter_out = []
        filter_types = {
                "intergenic":["intergenic_region"],
                "intragenic":["intragenic_variant"],
                "noncoding":["non_coding_transcript_variant","non_coding_transcript_exon_variant"],
                "downstream":["downstream_gene_variant"],
                "upstream":["upstream_gene_variant"],
                "intron":["intron_variant"],
                "synonymous":["synonymous_variant"],
                "splice":["splice_region_variant&intron_variant","splice_region_variant&synonymous_variant"],
                "ablation":["transcript_ablation"],
                "utr":["5_prime_UTR_variant","3_prime_UTR_variant"],
                "other":["feature_ablation"]
            }
        if exclude_variant_types:
            for t in exclude_variant_types:
                filter_out += filter_types[t]
        if keep_variant_types:
            for t in filter_types:
                if t in keep_variant_types: continue
                filter_out += filter_types[t]

        if bed_file:
            genes_to_keep = set()
            for l in open(bed_file):
                row = l.strip().split()
                genes_to_keep.add(row[3])
                genes_to_keep.add(row[4])

        variants = []
        vcf = pysam.VariantFile(self.filename)
        for var in vcf:
            logging.debug(str(var)[:5000])
            chrom = var.chrom
            pos = var.pos
            ref = var.ref
            alleles = var.alleles
            if "SVTYPE" in var.info:
                if var.info['SVTYPE']!="DEL":
                    continue
                ad = get_sv_ad(var)
                varlen = var.stop - var.pos
                sv = True
            else:
                ad = [int(x) for x in var.samples[0]['AD']]
                varlen = None
                sv = False
            if sum(ad)==0:
                continue
            af_dict = {alleles[i]:ad[i]/sum(ad) for i in range(len(alleles))}
            if "ANN" not in var.info:
                ## TODO - add a way to get the annotations from the VCF
                continue
            ann_strs = var.info['ANN']
            ann_list = [x.split("|") for x in ann_strs]
            for alt in alleles[1:]:
                strand_support = get_stand_support(var,alt)
                if strand_support[0] != None:
                    freq = sum(strand_support)/sum(ad)
                else:
                    if 'QNAME' in var.info:
                        freq = 1.0
                    else:
                        freq = af_dict[alt]

                    
                tmp_var = Variant(
                    chrom = chrom,
                    pos = int(pos),
                    ref = ref,
                    alt = alt,
                    depth = sum(ad),
                    freq = freq,
                    forward_reads = strand_support[0],
                    reverse_reads = strand_support[1],
                    sv = sv,
                    sv_len = varlen,
                )

                tmp_var.filter = filter_variant(tmp_var,filter_params)
                if tmp_var.filter=="hard_fail":
                    continue
                

                for ann in ann_list:
                    ann[3] = ann[3].replace("gene:","")
                    ann[4] = ann[4].replace("gene:","")
                    if ann[0]!=alt:
                        continue
                    if ann[1] in filter_out:
                        continue
                    if bed_file:
                        if ann[3] in genes_to_keep or ann[4] in genes_to_keep:
                            pass
                        else:
                            continue
                    if ann[1]=="upstream_gene_variant":
                        r = re.search("[cn].-([0-9]+)",ann[9])
                        if int(r.group(1))>max_promoter_length:
                            continue

                    tmp_var.consequences.append(
                        Consequence(
                            gene_name = ann[3],
                            gene_id = ann[4],
                            feature_id = ann[6],
                            type = ann[1],
                            nucleotide_change = ann[9],
                            protein_change = ann[10],
                            sequence_hgvs= f'{chrom}:g.{pos}{ref}>{alt}'
                        )
                    )

                # genes_in_csq = set([x.gene_id for x in tmp_var.consequences])
                # for gene in genes_in_csq:
                #     tmp_var.consequences.append(
                #         Consequence(
                #             gene_name = None,
                #             gene_id = gene,
                #             feature_id = None,
                #             type = "genomic_change",
                #             nucleotide_change = f'{chrom}:g.{pos}{ref}>{alt}',
                #             protein_change = None
                #         )
                #     )
                
                
                if len(tmp_var.consequences)==0:
                    continue
                variants.append(tmp_var)
        # variants = uniqify_dict_list(variants)            
        return variants


    
    def add_annotations(self,ref_file,bam_file):
        add_arguments_to_self(self,locals())
        self.new_file = self.prefix + ".ann.vcf.gz"

        run_cmd("gatk VariantAnnotator -R %(ref_file)s -I %(bam_file)s -V %(filename)s -O %(new_file)s -A MappingQualityRankSumTest -A ReadPosRankSumTest -A QualByDepth -A BaseQualityRankSumTest -A TandemRepeat -A StrandOddsRatio -OVI false" % vars(self))
        return Vcf(self.new_file,self.prefix)


    def get_positions(self) -> List[GenomePosition]:
        results = []
        for l in cmd_out("bcftools query -f '%%CHROM\\t%%POS\\n' %s" % self.filename):
            row = l.split()
            results.append(GenomePosition(chrom=row[0],posposition=int(row[1])))
        return results

    def get_bed_gt(self,bed_file: str,ref_file: str) -> Dict[GenomePosition,Dict[str,int]]:
        self.bed_file = bed_file
        self.ref_file = ref_file
        cmd = f"bcftools view -T {bed_file}  {self.filename}" + r" | bcftools query -u -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n'"
        results = defaultdict(dict)
        ref_seq = Fasta(ref_file).fa_dict
        for l in cmd_out(cmd):
            #Chromosome    4348079    0/0    51
            chrom,pos,ref,alt,gt,ad = l.rstrip().split()
            p = GenomePosition(chrom=chrom,pos=int(pos))
            d = {}
            alts = alt.split(",")
            ad = [int(x) for x in ad.split(",")] if ad!="." else [0,100]
            if gt=="0/0":
                d[ref] = ad[0]
            elif gt=="./.":
                d[ref] = 0
            else:
                for i,a in enumerate([ref]+alts):
                    d[a] = ad[i]
            results[p] = d
        bed = load_bed(bed_file)
        for r in bed:
            p = GenomePosition(chrom=r.chrom,pos=r.end)
            if p not in results:
                results[p] = {ref_seq[p.chrom][p.pos-1]:50}
        return results

    def get_gatk_annotations(self):
        possible_annotations = ["ReadPosRankSum","QD","MQRankSum","ClippingRankSum","BaseQRankSum","SOR","TandemRepeat"]
        lines = []
        for l in cmd_out("bcftools view -h %s | grep INFO" % self.filename):
            lines.append(l)
        lines = "".join(lines)
        found_annotations = []
        for x in possible_annotations:
            if x in lines:
                found_annotations.append(x)
        return found_annotations
    
    def get_vcf_qc(self):
        for l in cmd_out(f'bcftools stats {self.filename}'):
            row = l.strip().split("\t")
            if row[0]=='SN' and 'number of records' in l:
                num_variants = int(l.split("\t")[3].strip())
        return VcfQC(total_variants=num_variants)

# class DellyVcf(Vcf):
#     def __init__(self,filename):
#         Vcf.__init__(self,filename)
#     def get_robust_calls(self,prefix,bed_file = None):
#         self.tmpfile = f"{prefix}.tmp.delly.vcf.gz"
#         self.outfile = f"{prefix}.delly.vcf.gz"
#         run_cmd("bcftools view -c 2  %(filename)s | bcftools view -e '(INFO/END-POS)>=100000' -Oz -o %(tmpfile)s" % vars(self))
#         if bed_file:
#             run_cmd(f"bcftools index {self.tmpfile}")
#             run_cmd(f"bcftools view -R {bed_file} {self.tmpfile} -Oz -o {self.outfile}")
#         else:
#             self.outfile = self.tmpfile
#         return DellyVcf(self.outfile)

def uniqify_dict_list(data):
    s = []
    for obj in data:
        t = json.dumps(obj)
        if t not in s:
            s.append(t)
    
    return [json.loads(d) for d in s]


def var_qc_test(var: Variant,min_depth: int,min_af: float,strand_support: int) -> bool:
    """Test if a variant passes QC"""
    fail = False
    if min_depth!=None and var.depth<min_depth:
        fail = True
    if min_af!=None and var.freq<min_af:
        fail = True
    if strand_support!=None and var.forward_reads!=None and var.forward_reads<strand_support:
        fail = True
    if strand_support!=None and var.reverse_reads!=None and var.reverse_reads<strand_support:
        fail = True
    return fail

def sv_var_qc_test(var: Variant,min_depth: int,min_af: float, sv_len: int) -> bool:
    
    fail = False
    if min_depth!=None and var.depth<min_depth:
        fail = True
    if min_af!=None and var.freq<min_af:
        fail = True
    if sv_len!=None and var.sv_len>sv_len:
        fail = True
    return fail

def filter_variant(var,filter_params):
    """
    Filter a variant based on the filter parameters
    
    Parameters
    ----------
    var : Variant
        The variant to filter
    filter_params : dict
        The filter parameters
    
    Returns
    -------
    str
        The QC status of the variant
    
    Examples
    --------
    >>> from pathogenprofiler import Variant, generate_example_variant
    >>> filter_params = {
    ...     "depth_hard":5,
    ...     "depth_soft":10,
    ...     "af_hard":0,
    ...     "af_soft":0.1,
    ...     "strand_hard":0,
    ...     "strand_soft":3
    ... }
    >>> var = generate_example_variant(forward_reads=1,reverse_reads=1)
    >>> filter_variant(var,filter_params)
    'hard_fail'
    >>> var = generate_example_variant(forward_reads=3,reverse_reads=3)
    >>> filter_variant(var,filter_params)
    'soft_fail'
    >>> var = generate_example_variant(forward_reads=10,reverse_reads=10)
    >>> filter_variant(var,filter_params)
    'pass'
    """
    qc = "pass"
    if var.sv==True:
        if sv_var_qc_test(var,filter_params["sv_depth_hard"],filter_params["sv_af_hard"],filter_params["sv_len_hard"]):
            qc = "hard_fail"
        elif sv_var_qc_test(var,filter_params["sv_depth_soft"],filter_params["sv_af_soft"],filter_params["sv_len_soft"]):
            qc = "soft_fail"
    else:
        if var_qc_test(var,filter_params["depth_hard"],filter_params["af_hard"],filter_params["strand_hard"]):
            qc = "hard_fail"
        elif var_qc_test(var,filter_params["depth_soft"],filter_params["af_soft"],filter_params["strand_soft"]):
            qc = "soft_fail"
    logging.debug(qc)
    return qc
