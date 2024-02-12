from .utils import TempFilePrefix, run_cmd, cmd_out,add_arguments_to_self,rm_files, index_bcf,tabix, load_bed
import logging
from .fasta import Fasta
from collections import defaultdict
import re
from uuid import uuid4
import sys
import os.path
import json
import pysam
from .models import Variant,Consequence, VcfQC, GenomePosition
from typing import Optional, List, Tuple, Dict

# def get_stand_support(var,alt,caller):
#     alt_index = list(var.alts).index(alt)
#     forward_support = None
#     reverse_support = None
#     if caller == "freebayes" or caller=="lofreq":
#         forward_support = var.info['SAF'][alt_index]
#         reverse_support = var.info['SAR'][alt_index]
#     elif caller=="bcftools":
#         forward_support = var.samples[0]['ADF'][alt_index+1]
#         reverse_support = var.samples[0]['ADR'][alt_index+1]
#     else:
#         forward_support = None
#         reverse_support = None
#     return forward_support, reverse_support

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
    else:
        forward_support = None
        reverse_support = None
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


    def set_snpeff_datadir(self):
        """
            Look for snpEff database directory and if not found, set it to current working directory.
            Store the directory found in self.snpeff_data_dir
        """
        snpeff_basedir = None
        for path_el in os.environ.get('PATH','').split(os.pathsep):
            path_el = path_el.strip('"')

            fpath = os.path.join(path_el, 'snpEff')

            if os.path.isfile(fpath) and os.access(fpath, os.X_OK):
                snpeff_basedir = path_el.rstrip('/').rstrip('bin')
                break
        
        snpeff_data_dir = None
        if snpeff_basedir is not None:
            snpeff_shared_dir = os.path.join(snpeff_basedir, 'share')
            for path_el in os.listdir(snpeff_shared_dir):
                if path_el.startswith('snpeff-'):
                    snpeff_data_dir = os.path.join(snpeff_shared_dir, path_el, 'data')
                    snpeff_db_dir = os.path.join(snpeff_data_dir, self.db)
                    if not os.path.isdir(snpeff_data_dir) and os.access(os.path.join(snpeff_shared_dir, path_el), os.W_OK | os.X_OK | os.R_OK):
                        os.mkdir(snpeff_data_dir)
                    if (os.path.isdir(snpeff_db_dir) or 
                        (os.path.isdir(snpeff_data_dir) and os.access(snpeff_data_dir, os.W_OK | os.X_OK | os.R_OK))
                        ):
                        self.snpeff_data_dir = snpeff_data_dir
                        return snpeff_data_dir

        # if we have got here, we need to try an store the snpEff DB in the current working directory
        snpeff_data_dir = os.getcwd()
        if os.access(snpeff_data_dir, os.W_OK | os.R_OK):
            self.snpeff_data_dir = snpeff_data_dir
            return snpeff_data_dir
        
        return None

    def run_snpeff(self,db,ref_file,gff_file,rename_chroms = None, split_indels=True):
        logging.info("Running snpEff")
        add_arguments_to_self(self,locals())
        self.vcf_csq_file = self.prefix+".csq.vcf.gz"
        self.rename_cmd = f"rename_vcf_chrom.py --source {' '.join(rename_chroms['source'])} --target {' '.join(rename_chroms['target'])} |" if rename_chroms else ""
        self.re_rename_cmd = f"| rename_vcf_chrom.py --source {' '.join(rename_chroms['target'])} --target {' '.join(rename_chroms['source'])}" if rename_chroms else ""
        if self.set_snpeff_datadir() is None:
            logging.warning("snpEff database not found and no writeable directory to store database in, analysis might fail", file=sys.stderr)
            self.snpeff_data_dir_opt = ''
        else:
            self.snpeff_data_dir_opt = '-dataDir %(snpeff_data_dir)s' % vars(self)
        if split_indels:
            with TempFilePrefix() as tmp:
                self.tmp_file1 = f"{tmp}.1.vcf.gz"
                self.tmp_file2 = f"{tmp}.2.vcf.gz"
                self.tmp_file3 = f"{tmp}.3.vcf.gz"

                run_cmd("bcftools view -c 1 -a %(filename)s | bcftools view -v snps | combine_vcf_variants.py --ref %(ref_file)s --gff %(gff_file)s | %(rename_cmd)s snpEff ann %(snpeff_data_dir_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools sort -Oz -o %(tmp_file1)s && bcftools index %(tmp_file1)s" % vars(self))
                run_cmd("bcftools view -c 1 -a %(filename)s | bcftools view -v indels | %(rename_cmd)s snpEff ann %(snpeff_data_dir_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools sort -Oz -o %(tmp_file2)s && bcftools index %(tmp_file2)s" % vars(self))
                run_cmd("bcftools view -c 1 -a %(filename)s | bcftools view -v other | %(rename_cmd)s snpEff ann %(snpeff_data_dir_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools sort -Oz -o %(tmp_file3)s && bcftools index %(tmp_file3)s" % vars(self))
                run_cmd("bcftools concat -a %(tmp_file1)s %(tmp_file2)s %(tmp_file3)s | bcftools sort -Oz -o %(vcf_csq_file)s" % vars(self))
        else :
            run_cmd("bcftools view %(filename)s | %(rename_cmd)s snpEff ann %(snpeff_data_dir_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools view -Oz -o %(vcf_csq_file)s" % vars(self))
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
            logging.debug(var)
            chrom = var.chrom
            pos = var.pos
            ref = var.ref
            alleles = var.alleles
            alt_str = list(var.alts)[0]
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
            ann_strs = var.info['ANN']
            ann_list = [x.split("|") for x in ann_strs]
            for alt in alleles[1:]:
                strand_support = get_stand_support(var,alt)

                tmp_var = Variant(
                    chrom = chrom,
                    pos = int(pos),
                    ref = ref,
                    alt = alt,
                    depth = sum(ad),
                    freq = af_dict[alt],
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
                            protein_change = ann[10]
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
            print(l)
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

class DellyVcf(Vcf):
    def __init__(self,filename):
        Vcf.__init__(self,filename)
    def get_robust_calls(self,prefix,bed_file = None):
        self.tmpfile = f"{prefix}.tmp.delly.vcf.gz"
        self.outfile = f"{prefix}.delly.vcf.gz"
        run_cmd("bcftools view -c 2  %(filename)s | bcftools view -e '(INFO/END-POS)>=100000' -Oz -o %(tmpfile)s" % vars(self))
        if bed_file:
            run_cmd(f"bcftools index {self.tmpfile}")
            run_cmd(f"bcftools view -R {bed_file} {self.tmpfile} -Oz -o {self.outfile}")
        else:
            self.outfile = self.tmpfile
        return DellyVcf(self.outfile)

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
    return qc
