from glob import glob
from .kmer import KmerDump
from .utils import add_arguments_to_self, run_cmd, cmd_out, filecheck, index_bam, run_cmd_parallel_on_genome
from .vcf import Vcf
from collections import defaultdict
import json
from uuid import uuid4
import os
import platform 
import statistics as stats
import logging 

class Bam:
    """
    A class to perform operations on BAM files such as SNP calling
    """
    def __init__(self,bam_file,prefix,platform,threads=1):
        add_arguments_to_self(self, locals())
        logging.debug("Creating Bam object with %s" % (bam_file))
        filecheck(self.bam_file)
        index_bam(bam_file,threads=threads)
        self.filetype = "cram" if bam_file[-5:]==".cram" else "bam"
        

    def run_delly(self,bed_file):
        logging.info("Running delly")
        self.bed_file = bed_file
        if self.platform=="illumina":
            run_cmd("delly call -t DEL -g %(ref_file)s %(bam_file)s -o %(prefix)s.delly.bcf" % vars(self))
            run_cmd("bcftools view -c 2 %(prefix)s.delly.bcf | bcftools view -e '(INFO/END-POS)>=100000' -Oz -o %(prefix)s.delly.vcf.gz" % vars(self))
            run_cmd("bcftools index %(prefix)s.delly.vcf.gz" % vars(self))
            run_cmd("bcftools view -R %(bed_file)s %(prefix)s.delly.vcf.gz -Oz -o %(prefix)s.delly.targets.vcf.gz" % vars(self))
            return Vcf("%(prefix)s.delly.targets.vcf.gz" % vars(self))
        else:
            run_cmd("delly lr -t DEL -g %(ref_file)s %(bam_file)s -o %(prefix)s.delly.bcf" % vars(self))
            return Vcf("%(prefix)s.delly.bcf" % vars(self))
                
    def call_variants(self,ref_file,caller,filters,bed_file=None,threads=1,calling_params=None, samclip=False):
        add_arguments_to_self(self, locals())
        filecheck(ref_file)
        self.caller = caller.lower()
        self.af_hard = filters['af_hard'] if 'af_hard' in filters else 0
        # Set up final vcf file name
        # Make the windows for parallel calling based on chunking the whole
        # genome or by providing a bed file
        if bed_file:
            self.vcf_file = "%s.targets.vcf.gz" % (self.prefix) if bed_file else "%s.vcf.gz" % (self.prefix)
            self.windows_cmd = "cat %(bed_file)s | awk '{print $1\":\"$2\"-\"$3\" \"$1\"_\"$2\"_\"$3}'" % vars(self)
        else:
            self.vcf_file = "%s.vcf.gz" % (self.prefix) if bed_file else "%s.vcf.gz" % (self.prefix)
            self.windows_cmd = "get_genome_regions.py --fasta %(ref_file)s --num %(threads)s" % vars(self)

        self.samclip_cmd = "| samclip --ref %(ref_file)s" % vars(self) if samclip else ""
        # Run through different options. 
        if self.platform == "nanopore" and self.caller=="bcftools":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD,ADF,ADR -r {region} %(bam_file)s | bcftools call -mv | annotate_maaf.py | bcftools +fill-tags | bcftools view -c 1 | bcftools norm -f %(ref_file)s | bcftools filter -e 'IMF < 0.7' -S 0 -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform == "nanopore" and self.caller=="freebayes":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "freebayes -f %(ref_file)s -F %(af_hard)s -r {region} --haplotype-length -1 %(calling_params)s %(bam_file)s | annotate_maaf.py | bcftools +fill-tags | bcftools view -c 1 | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="nanopore" and self.caller == "pilon":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = """pilon --genome %(ref_file)s --targets {region} %(calling_params)s --nanopore %(bam_file)s --variant --output %(prefix)s.{region_safe} && bcftools view -i 'AF>0' -c 1 %(prefix)s.{region_safe}.vcf |add_dummy_AD.py | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz""" % vars(self)
        elif self.platform == "pacbio" and self.caller=="freebayes":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "freebayes -f %(ref_file)s -r {region} --haplotype-length -1 %(calling_params)s %(bam_file)s | annotate_maaf.py | bcftools +fill-tags | bcftools view -c 1 | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="pacbio" and self.caller == "pilon":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = """pilon --genome %(ref_file)s --targets {region} %(calling_params)s --pacbio %(bam_file)s --variant --output %(prefix)s.{region_safe} && bcftools view -i 'AF>0' -c 1 %(prefix)s.{region_safe}.vcf |add_dummy_AD.py | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz""" % vars(self)
        elif self.platform=="illumina" and self.caller == "bcftools":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(prefix)s.{region_safe}.tmp.bam && samtools index %(prefix)s.{region_safe}.tmp.bam && bcftools mpileup -f %(ref_file)s %(calling_params)s -a DP,AD,ADF,ADR -r {region} %(prefix)s.{region_safe}.tmp.bam | bcftools call -mv | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="illumina" and self.caller == "gatk":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(prefix)s.{region_safe}.tmp.bam && samtools index %(prefix)s.{region_safe}.tmp.bam && gatk HaplotypeCaller -R %(ref_file)s -I %(prefix)s.{region_safe}.tmp.bam -O /dev/stdout -L {region} %(calling_params)s -OVI false | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="illumina" and self.caller == "freebayes":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "samtools view -T %(ref_file)s -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(prefix)s.{region_safe}.tmp.bam && samtools index %(prefix)s.{region_safe}.tmp.bam && freebayes -f %(ref_file)s -r {region} --haplotype-length -1 %(calling_params)s %(prefix)s.{region_safe}.tmp.bam | bcftools view -c 1 | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="illumina" and self.caller == "pilon":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "samtools view -T %(ref_file)s -f 0x1 -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(prefix)s.{region_safe}.tmp.bam && samtools index %(prefix)s.{region_safe}.tmp.bam && pilon --genome %(ref_file)s --targets {region} --diploid %(calling_params)s --frags %(prefix)s.{region_safe}.tmp.bam --variant --output %(prefix)s.{region_safe} && bcftools view -e " % vars(self) +  r"""'ALT="."'""" +  " -c 1 %(prefix)s.{region_safe}.vcf | add_dummy_AD.py | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)
        elif self.platform=="illumina" and self.caller == "lofreq":
            self.calling_params = calling_params if calling_params else ""
            self.calling_cmd = "samtools view -T %(ref_file)s  -h %(bam_file)s {region} %(samclip_cmd)s | samtools view -b > %(prefix)s.{region_safe}.tmp.bam && samtools index %(prefix)s.{region_safe}.tmp.bam && lofreq call --call-indels -f %(ref_file)s -r {region} %(calling_params)s  %(prefix)s.{region_safe}.tmp.bam  | modify_lofreq_vcf.py | add_dummy_AD.py --ref %(ref_file)s --sample-name %(prefix)s --add-dp | bcftools view -c 1 | bcftools norm -f %(ref_file)s -Oz -o %(prefix)s.{region_safe}.vcf.gz" % vars(self)
        else:
            logging.debug("Unknown combination %(platform)s + %(caller)s" % vars(self))
        logging.info("Running variant calling")
        run_cmd_parallel_on_genome(self.calling_cmd,ref_file,bed_file = bed_file,threads=threads,desc="Calling variants")
        cmd = "bcftools index  %(prefix)s.{region_safe}.vcf.gz" % vars(self) 
        run_cmd_parallel_on_genome(cmd,ref_file,bed_file = bed_file,threads=threads,desc="Indexing variants")
        run_cmd("bcftools concat -aD `%(windows_cmd)s | awk '{print \"%(prefix)s.\"$2\".vcf.gz\"}'` | bcftools view -Oz -o %(vcf_file)s" % vars(self))
        run_cmd("rm `%(windows_cmd)s | awk '{print \"%(prefix)s.\"$2\".vcf.gz*\"}'`" % vars(self))
        if self.platform=="illumina":
            run_cmd("rm `%(windows_cmd)s | awk '{print \"%(prefix)s.\"$2\".tmp.bam*\"}'`" % vars(self))

        return Vcf(self.vcf_file)
    
    def get_median_depth(self,ref_file,software="bedtools"):
        logging.info("Calculating median depth")
        if software=="bedtools":
            lines = []
            for l in cmd_out("bedtools genomecov -ibam %s" % (self.bam_file)):
                arr = l.split()
                if arr[0]=="genome":
                    lines.append(arr)
            midpoint =  int(lines[0][3])/2
            x = 0
            for row in lines:
                x = x + int(row[2])
                if x>midpoint:
                    break
            self.median_coverage = int(row[1])
            return int(row[1])
        elif software=="mosdepth":
            self.median_coverage = None
            tmp = str(uuid4())
            run_cmd(f"mosdepth {tmp} {self.bam_file} -f {ref_file}")
            for l in open(f"{tmp}.mosdepth.summary.txt"):
                row = l.strip().split()
                if row[0]=="total":
                    self.median_coverage = int(row[1])
            for f in glob(f"{tmp}*"):
                os.remove(f)
            return int(float(row[3]))

    def calculate_median_coverage(self,ref_file,software="bedtools"):
        logging.info("Calculating median depth")
        if software=="bedtools":
            lines = []
            for l in cmd_out("bedtools genomecov -ibam %s" % (self.bam_file)):
                arr = l.split()

                if arr[0]=="genome":
                    lines.append(arr)
            midpoint =  int(lines[0][3])/2
            x = 0
            for row in lines:
                x = x + int(row[2])
                if x>midpoint:
                    break
            self.median_coverage = int(row[1])
            return int(row[1])
        elif software=="mosdepth":
            self.median_coverage = None
            tmp = str(uuid4())
            run_cmd(f"mosdepth {tmp} {self.bam_file} -f {ref_file}")
            for l in open(f"{tmp}.mosdepth.summary.txt"):
                row = l.strip().split()
                if row[0]=="total":
                    self.median_coverage = int(row[1])
            for f in glob(f"{tmp}*"):
                os.remove(f)
            return int(float(row[3]))


    def get_bed_gt(self,bed_file,ref_file,caller,platform):
        logging.info("Getting genotypes for positions in bed file")
        add_arguments_to_self(self, locals())
        results = defaultdict(lambda : defaultdict(dict))
        run_cmd("samtools view -Mb -L %(bed_file)s %(bam_file)s -T %(ref_file)s > %(prefix)s.tmp.bam" % vars(self))
        run_cmd("samtools index %(prefix)s.tmp.bam" % vars(self))
        if platform in ("nanopore","pacbio"):
            caller="bcftools"
        if caller == "gatk":
            cmd = "gatk HaplotypeCaller -I %(prefix)s.tmp.bam -R %(ref_file)s -L %(bed_file)s -OVI false -O /dev/stdout | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        elif caller == "freebayes":
            cmd = "freebayes -f %(ref_file)s -t %(bed_file)s %(prefix)s.tmp.bam --haplotype-length -1 | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        elif caller == "bcftools":
            cmd = "bcftools mpileup -f %(ref_file)s -T %(bed_file)s %(prefix)s.tmp.bam -BI -a AD | bcftools call -mv | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)
        else:
            cmd = "freebayes -f %(ref_file)s -t %(bed_file)s %(prefix)s.tmp.bam --haplotype-length -1 | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%GT\\t%%AD]\\n'" % vars(self)

        for l in cmd_out(cmd):
            # Chromosome    4348079    0/0    51
            chrom, pos, ref, alt, gt, ad = l.rstrip().split()
            pos = int(pos)
            d = {}
            alts = alt.split(",")
            ad = [int(x) for x in ad.split(",")]
            # if gt == "0/0":
            #     d[ref] = ad[0]
            # elif gt == "./.":
            if gt == "./.":
                d[ref] = 0
            else:
                genotypes = list([ref]+alts)
                if platform in ("nanopore","pacbio"):
                    idx = ad.index(max(ad))
                    d[genotypes[idx]] = ad[idx]
                else:
                    for i, a in enumerate(genotypes):
                        d[a] = ad[i]
            results[chrom][pos] = d

        ref_nt = {}
        for l in cmd_out("bedtools getfasta -fi %s -bed %s" % (ref_file,bed_file)):
            if l[0]==">":
                tmp = l.strip().replace(">","").split(":")
                tmp_chrom = tmp[0]
                tmp_pos = int(tmp[1].split("-")[1])
            else:
                ref_nt[(tmp_chrom,tmp_pos)] = l.strip().upper()

        for l in cmd_out(f"samtools view -b -L {bed_file} {self.prefix}.tmp.bam | bedtools coverage -a {bed_file} -b - -d -sorted"):
            row = l.strip().split()
            chrom = row[0]
            pos = int(row[2])
            cov = int(row[-1])
            if chrom not in results or pos not in results[chrom]:
                results[chrom][pos] = {ref_nt[(chrom,pos)]:cov}

        os.remove(f"{self.prefix}.tmp.bam")
        os.remove(f"{self.prefix}.tmp.bam.bai")

        return results

    def calculate_region_coverage(self,bed_file,depth_threshold=0,region_column=None):
        logging.info("Calculating coverage for regions in bed file")

        add_arguments_to_self(self, locals())
        numrows =len(open(bed_file).readline().split())
        if region_column is None:
            if numrows==7:
                region_column = 7
            else:
                region_column = 4
        self.region_cov = defaultdict(list)
        self.region_qc = []
        self.genome_coverage = []

        for l in cmd_out(f"samtools view -Mb -L {bed_file} {self.bam_file} | bedtools coverage -a {bed_file} -b - -d -sorted"):
            row = l.split()
            region = row[region_column-1]
            depth = int(row[-1])
            genomic_position = int(row[1]) + int(row[-2]) -1
            self.genome_coverage.append((genomic_position, depth))
            self.region_cov[region].append(depth)

        for region in self.region_cov:
            region_len = len(self.region_cov[region])
            pos_pass_thresh = len([d for d in self.region_cov[region] if d>=depth_threshold])
            self.region_qc.append({
                "region":region, 
                "pct_depth_pass":round(pos_pass_thresh/region_len*100,2), 
                "median_depth":stats.median(self.region_cov[region]),
            })

    def calculate_bamstats(self):
        logging.info("Calculating bamstats")
        temp_file = str(uuid4())
        run_cmd(f"samtools stats {self.bam_file} > {temp_file}")
        for l in open(temp_file):
            row = l.strip().split("\t")
            if row[0]!="SN": continue
            if row[1]=='raw total sequences:': self.total_reads = int(row[2])
            if row[1]=='reads mapped:': self.mapped_reads = int(row[2])
        self.pct_reads_mapped = round(self.mapped_reads/self.total_reads*100,2)
        os.remove(temp_file)
    
    def get_missing_genomic_positions(self,bed_file=None,cutoff=10):
        logging.info("Getting missing genomic positions")
        if not hasattr(self,"genome_coverage"):
            self.calculate_region_coverage(bed_file)
        return [x[0] for x in self.genome_coverage if x[1]<cutoff]

    def get_region_qc(self,bed_file=None,cutoff=10):
        logging.info("Getting qc metrics for regions")
        if not hasattr(self,"region_qc"):
            self.calculate_region_coverage(bed_file,depth_threshold=cutoff)
        return self.region_qc
    
    def get_kmer_counts(self,prefix,klen = 31,threads=1,max_mem=2,counter = "kmc"):
        logging.info("Getting kmer counts")
        if counter=="kmc":
            if threads>32:
                threads = 32
            tmp_prefix = str(uuid4())
            os.mkdir(tmp_prefix)
            bins = "-n128" if platform.system()=="Darwin" else ""
            run_cmd(f"kmc {bins} -fbam -m{max_mem} -t{threads} -k{klen} {self.bam_file} {tmp_prefix} {tmp_prefix}")
            run_cmd(f"kmc_dump {tmp_prefix} {prefix}.kmers.txt")
            run_cmd(f"rm -r {tmp_prefix}*")
            return KmerDump(f"{prefix}.kmers.txt",counter)
        else:
            logging.error("Can't use dsk for bam files, please use kmc instead")