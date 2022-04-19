from .utils import run_cmd, cmd_out,add_arguments_to_self,rm_files, index_bcf,tabix, log, load_bed, debug, warninglog
from .fasta import fasta
from collections import defaultdict
import re
from uuid import uuid4
import sys
import os.path

# re_seq = re.compile("([0-9\-]*)([A-Z\*]+)")
# re_I = re.compile("([A-Z\*]+)")
# number_re = re.compile("[0-9\-]+")

# def parse_mutation(x):
#     tmp = x.split(">")
#     aa_changed = True if len(tmp)>1 else False
#     re_obj = re_seq.search(tmp[0])
#     change_num = re_obj.group(1)
#     ref_aa = re_obj.group(2)
#     alt_aa = re_seq.search(tmp[1]).group(2) if aa_changed else None
#     return change_num,ref_aa,alt_aa

class vcf:
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
        for l in cmd_out("bcftools query -l %(filename)s" % vars(self)):
            self.samples.append(l.rstrip())
    def view_regions(self,bed_file):
        self.bed_file = bed_file
        self.newfile = "%(prefix)s.targets.vcf.gz" % vars(self)
        run_cmd("bcftools view -R %(bed_file)s %(filename)s -Oz -o %(newfile)s" % vars(self))
        return vcf(self.newfile)


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
        add_arguments_to_self(self,locals())
        self.vcf_csq_file = self.prefix+".csq.vcf.gz"
        self.rename_cmd = f"rename_vcf_chrom.py --source {' '.join(rename_chroms['source'])} --target {' '.join(rename_chroms['target'])} |" if rename_chroms else ""
        self.re_rename_cmd = f"| rename_vcf_chrom.py --source {' '.join(rename_chroms['target'])} --target {' '.join(rename_chroms['source'])}" if rename_chroms else ""
        if self.set_snpeff_datadir() is None:
            print("WARNING: snpEff database not found and no writeable directory to store database in, analysis might fail", file=sys.stderr)
            self.snpeff_data_dir_opt = ''
        else:
            self.snpeff_data_dir_opt = '-dataDir %(snpeff_data_dir)s' % vars(self)
        if split_indels:
            self.tmp_file1 = "%s.vcf.gz" % uuid4()
            self.tmp_file2 = "%s.vcf.gz" % uuid4()

            run_cmd("bcftools norm -m - %(filename)s | bcftools view -v snps | combine_vcf_variants.py --ref %(ref_file)s --gff %(gff_file)s | %(rename_cmd)s snpEff ann %(snpeff_data_dir_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools sort -Oz -o %(tmp_file1)s && bcftools index %(tmp_file1)s" % vars(self))
            run_cmd("bcftools norm -m - %(filename)s | bcftools view -v indels | %(rename_cmd)s snpEff ann %(snpeff_data_dir_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools sort -Oz -o %(tmp_file2)s && bcftools index %(tmp_file2)s" % vars(self))
            run_cmd("bcftools concat -a %(tmp_file1)s %(tmp_file2)s | bcftools sort -Oz -o %(vcf_csq_file)s" % vars(self))
            rm_files([self.tmp_file1, self.tmp_file2, self.tmp_file1+".csi", self.tmp_file2+".csi"])
        else :
            run_cmd("bcftools view %(filename)s | %(rename_cmd)s snpEff ann %(snpeff_data_dir_opt)s -noLog -noStats %(db)s - %(re_rename_cmd)s | bcftools view -Oz -o %(vcf_csq_file)s" % vars(self))
        return vcf(self.vcf_csq_file,self.prefix)



    def load_ann(self,max_promoter_length=1000, bed_file=None,exclude_variant_types = None,keep_variant_types=None):
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
                "utr":["5_prime_UTR_variant","3_prime_UTR_variant"]
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
        for l in cmd_out(f"bcftools query -u -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%ANN\\t[%AD]\\n' {self.filename}"):
            chrom,pos,ref,alt_str,ann_str,ad_str = l.strip().split()
            
            alleles = [ref] + alt_str.split(",")
            if alt_str=="<DEL>":
                af_dict = {"<DEL>":1.0}
            else:
                ad = [int(x) for x in ad_str.split(",")]
                af_dict = {alleles[i]:ad[i]/sum(ad) for i in range(len(alleles))}
            ann_list = [x.split("|") for x in ann_str.split(",")]
            for alt in alleles[1:]:
                tmp_var = {
                    "chrom": chrom,
                    "genome_pos": int(pos),
                    "ref": ref,
                    "alt":alt,
                    "freq":af_dict[alt],
                    "consequences":[]
                }
                # if pos=="1673425":
                    # import pdb; pdb.set_trace()
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

                    tmp = {
                        "gene_name":ann[3],
                        "gene_id":ann[4],
                        "feature_id":ann[6],
                        "type":ann[1],
                        "nucleotide_change":ann[9],
                        "protein_change":ann[10],
                    }
                    tmp_var["consequences"].append(tmp)
                variants.append(tmp_var)
                    
        return variants

    
    def add_annotations(self,ref_file,bam_file):
        add_arguments_to_self(self,locals())
        self.new_file = self.prefix + ".ann.vcf.gz"

        run_cmd("gatk VariantAnnotator -R %(ref_file)s -I %(bam_file)s -V %(filename)s -O %(new_file)s -A MappingQualityRankSumTest -A ReadPosRankSumTest -A QualByDepth -A BaseQualityRankSumTest -A TandemRepeat -A StrandOddsRatio -OVI false" % vars(self))
        return vcf(self.new_file,self.prefix)


    def get_positions(self):
        results = []
        for l in cmd_out("bcftools query -f '%%CHROM\\t%%POS\\n' %s" % self.filename):
            row = l.split()
            results.append((row[0],int(row[1])))
        return results

    def get_bed_gt(self,bed_file,ref_file):
        add_arguments_to_self(self,locals())
        cmd = "bcftools convert --gvcf2vcf -f %(ref_file)s %(filename)s  | bcftools view -T %(bed_file)s  | bcftools query -u -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        results = defaultdict(lambda : defaultdict(dict))
        ref_seq = fasta(ref_file).fa_dict
        for l in cmd_out(cmd):
            #Chromosome    4348079    0/0    51
            chrom,pos,ref,alt,gt,ad = l.rstrip().split()
            pos =int(pos)
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
            results[chrom][pos] = d
        bed = load_bed(bed_file,[1,3,5],1,3)
        for chrom in bed:
            for pos in bed[chrom]:
                if int(pos) not in results[chrom]:
                    results[chrom][int(pos)] = {ref_seq[chrom][pos-1]:50}
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

class delly_bcf(vcf):
    def __init__(self,filename):
         vcf.__init__(self,filename)
    def get_robust_calls(self,prefix,bed_file = None):
        self.tmpfile = f"{prefix}.tmp.delly.vcf.gz"
        self.outfile = f"{prefix}.delly.vcf.gz"
        run_cmd("bcftools view -c 2  %(filename)s | bcftools view -e '(INFO/END-POS)>=100000' -Oz -o %(tmpfile)s" % vars(self))
        if bed_file:
            run_cmd(f"bcftools index {self.tmpfile}")
            run_cmd(f"bcftools view -R {bed_file} {self.tmpfile} -Oz -o {self.outfile}")
        else:
            self.outfile = self.tmpfile
        return delly_bcf(self.outfile)

